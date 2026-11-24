import os
import numpy as np
import petsc4py
from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import getfacevelocity
    from gospl._fortran import advecupwind

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()


class WaveMesh(object):
    """
    Wave energy and wave shadows (opt-in ``waves:`` block,
    docs/DESIGN_CARBONATES.md 2.10, phase C7).

    For each direction of the wave climate, a relative wave energy ``E``
    (0 to 1) is carried from the open ocean across the sea and lost on land
    and in shallow water:

    .. math::
        (\\hat{\\mathbf{d}} \\cdot \\nabla + \\kappa)\\, E = 0, \\qquad
        \\kappa = \\frac{1}{L_d} + \\frac{1}{L_s}
                  \\max\\left(0, 1 - \\frac{w}{w_s}\\right),

    with ``E = 1`` on the edges where the waves enter the domain, ``E = 0`` on
    land (identity rows: the shadow behind islands and headlands), ``w`` the
    water depth, ``L_d`` the attenuation length in deep water (``decay``),
    ``w_s`` the depth below which the waves feel the bottom (``shoal_depth``)
    and ``L_s`` the attenuation length at the shoreline (``shoal_length``).
    The climate is a list of directions with weights; ``E`` is their weighted
    sum. The operator is the first-order upwind finite-volume scheme of the
    orographic rain (``getfacevelocity`` + ``advecupwind``), whose
    coefficients are computed once per direction; the matrix is rebuilt every
    step because the coastline moves. Flat meshes only.

    The field is a response control (``wave``) of the carbonate factories and
    the energy of the coastal erosion (phase C4).
    """

    def __init__(self):
        """Cached coefficients per direction and the solver (registered in
        ``destroy_DMPlex``)."""
        self.waveE = None
        self._waveCoef = None
        self._waveMat = None
        self._ksp_wave = None
        if not getattr(self, "waveOn", False):
            return
        if not self.flatModel:
            raise ValueError(
                "waves: wave energy is only available on flat meshes for now "
                "(a global mesh needs a direction field)."
            )
        self.waveE = np.zeros(self.lpoints, dtype=np.float64)
        return

    def _waveCoefficients(self):
        """Upwind advection coefficients and inflow edges of every direction
        of the climate (unit speed, so the operator is ``d . grad`` in 1/m).
        Computed once (the mesh does not change)."""
        self._waveCoef = []
        for k, (frm, _) in enumerate(self.waveClimate):
            a = np.deg2rad(frm)
            u0, v0 = -np.sin(a), -np.cos(a)      # waves travel away from `frm`
            vel = np.zeros((self.lpoints, 3))
            vel[:, 0] = u0
            vel[:, 1] = v0
            getfacevelocity(self.lpoints, vel)
            lcoeff = advecupwind(self.lpoints, 1.0)
            inflow = self._oroInflowBorders(u0, v0)
            self._waveCoef.append((lcoeff[:, 0] - 1.0, lcoeff[:, 1:].copy(), inflow))

    def _waveSolve(self, diag, off, fixed, rhs):
        """Assemble ``diag`` + ``off`` with identity rows at ``fixed`` and
        solve with the right-hand side ``rhs`` (the imposed values on the
        fixed rows, 0 elsewhere); returns the local solution. Collective."""
        d = diag.copy()
        o = off.copy()
        d[fixed] = 1.0
        o[fixed, :] = 0.0
        coeffs = np.empty((self.lpoints, 1 + self.maxnb), dtype=np.float64)
        coeffs[:, 0] = d
        coeffs[:, 1:] = o[:, : self.maxnb]
        if self._waveMat is not None:
            self._waveMat.destroy()
        self._waveMat = self._assembleDiffMatCSR(coeffs)
        if self._ksp_wave is None:
            ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
            ksp.setType("gmres")
            ksp.getPC().setType("bjacobi")
            ksp.setTolerances(rtol=1.0e-10, atol=1.0e-14, max_it=2000)
            ksp.setOptionsPrefix("wave_")
            petsc4py.PETSc.Options()["wave_sub_pc_factor_shift_type"] = "nonzero"
            ksp.setFromOptions()
            self._ksp_wave = ksp
        self._ksp_wave.setOperators(self._waveMat, self._waveMat)
        self.tmpL.setArray(np.where(fixed, rhs, 0.0))
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.tmp1.set(0.0)
        self._ksp_wave.solve(self.tmp, self.tmp1)
        if self._ksp_wave.getConvergedReason() < 0 and MPIrank == 0:
            print("[waves] wave energy solve did not converge (%d)"
                  % self._ksp_wave.getConvergedReason(), flush=True)
        self.dm.globalToLocal(self.tmp1, self.tmpL)
        return self.tmpL.getArray().copy()

    def waveEnergy(self):
        """
        Relative wave energy ``waveE`` (0 to 1) of the current surface: the
        weighted sum over the wave climate of the solutions of
        ``(d . grad + kappa) E = 0`` (see the class). Collective. Scratch:
        ``tmp``, ``tmpL``, ``tmp1``.
        """
        t0 = process_time()
        if self._waveCoef is None:
            self._waveCoefficients()
        h = self.hLocal.getArray()
        w = self.sealevel - h
        land = w <= 0.0
        kappa = 1.0 / self.waveDecay + np.maximum(0.0, 1.0 - w / self.waveShoalDepth) \
            / self.waveShoalLength
        E = np.zeros(self.lpoints, dtype=np.float64)
        for (frm, wgt), (adiag, off, inflow) in zip(self.waveClimate, self._waveCoef):
            fixed = land.copy()
            val = np.zeros(self.lpoints)
            # waves enter through the inflow edges (open sea there)
            inn = np.zeros(self.lpoints, dtype=bool)
            inn[inflow] = True
            inn &= ~land
            fixed |= inn
            val[inn] = 1.0
            Ek = self._waveSolve(adiag + kappa, off, fixed, val)
            Ek[land] = 0.0
            E += wgt * np.clip(Ek, 0.0, 1.0)
        self.waveE = E
        if MPIrank == 0 and self.verbose:
            print("[waves] wave energy (%0.02f seconds)" % (process_time() - t0),
                  flush=True)
        return E
