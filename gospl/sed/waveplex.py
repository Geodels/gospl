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
    r"""
    Wave energy and wave shadows (opt-in ``waves:`` block,
    docs/DESIGN_CARBONATES.md 2.10, phase C7), on flat and global meshes.

    For each direction of the wave climate, a relative wave energy ``E``
    (0 to 1) travels across the sea along the local wave direction
    :math:`\hat{\mathbf{d}}`, grows over a fetch and is lost on land and in
    shallow water:

    .. math::
        \left(\hat{\mathbf{d}} \cdot \nabla + \frac{1}{L_f} + \kappa\right) E
        = \frac{1}{L_f}, \qquad
        \kappa = \frac{1}{L_d} + \frac{1}{L_s}
                  \max\left(0, 1 - \frac{w}{w_s}\right),

    with ``E = 0`` on land (identity rows: the shadow behind islands and
    headlands, from which the waves grow again over the fetch), ``w`` the
    water depth, ``L_f`` the fetch length over which the waves reach their
    full energy (``fetch``; without it there is no growth term), ``L_d`` the
    attenuation length in deep water (``decay``), ``w_s`` the depth below which
    the waves feel the bottom (``shoal_depth``) and ``L_s`` the attenuation
    length at the shoreline (``shoal_length``). On a flat mesh the domain
    edges facing the waves bring the open-ocean energy in (``E = 1``); a
    global mesh has no edges and needs the fetch.

    The local direction is a bearing (the direction the waves come from)
    turned into a tangent vector with the local east and north: one bearing
    everywhere, a map of bearings, or the ``zonal`` belts (trade winds from
    the east within 30 degrees of the equator, westerlies to 60 degrees,
    polar easterlies beyond). The climate is a list of directions with
    weights; ``E`` is their weighted sum. The operator is the first-order
    upwind finite-volume scheme of the orographic rain and of the tectonic
    advection (``getfacevelocity`` with a velocity per node +
    ``advecupwind``), whose coefficients are computed once per direction;
    the matrix is rebuilt every step because the coastline moves.

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
        if not self.flatModel and self.waveFetch is None:
            raise ValueError(
                "waves: a global mesh has no domain edges to bring the waves in; "
                "set `fetch` (m), the length over which the waves grow to full "
                "energy over open water."
            )
        if self.flatModel and any(f == "zonal" for f, _ in self.waveClimate):
            raise ValueError(
                "waves: the `zonal` belts need latitudes, i.e. a global mesh; on a "
                "flat mesh give a bearing or a from_map.")
        self.waveE = np.zeros(self.lpoints, dtype=np.float64)
        return

    def _waveBearings(self, frm):
        """Bearing (degrees, the direction the waves come from) at every
        local node for one climate entry: a number, ``zonal`` belts, or a
        ``(file, key)`` map."""
        if isinstance(frm, tuple):
            d = np.load(frm[0] + ".npz")
            arr = np.asarray(d[frm[1]], dtype=np.float64)
            if arr.shape[0] != self.mpoints:
                raise ValueError("waves: bearing map %s has %d values, the mesh %d."
                                 % (list(frm), arr.shape[0], self.mpoints))
            return arr[self.locIDs]
        if frm == "zonal":
            lat = self._waveLatitude()
            b = np.where(np.abs(lat) < 30.0, 90.0, 270.0)       # trades, westerlies
            return np.where(np.abs(lat) >= 60.0, 90.0, b)       # polar easterlies
        return np.full(self.lpoints, float(frm))

    def _waveLatitude(self):
        """Latitude (degrees) of the local nodes (global mesh); on a flat mesh
        the ``y`` coordinate is not a latitude, so the zonal belts are refused
        there by the parser."""
        xyz = self.lcoords
        r = np.linalg.norm(xyz, axis=1)
        return np.degrees(np.arcsin(np.clip(xyz[:, 2] / np.maximum(r, 1.0), -1.0, 1.0)))

    def _waveVelocity(self, bearing):
        """Unit wave-travel vectors per node from bearings (the direction the
        waves come from, 0 = north, 90 = east): along ``(x, y)`` on a flat
        mesh, in the local tangent plane (east, north) on a sphere."""
        a = np.deg2rad(bearing)
        ue, un = -np.sin(a), -np.cos(a)                     # travel away from `from`
        vel = np.zeros((self.lpoints, 3))
        if self.flatModel:
            vel[:, 0] = ue
            vel[:, 1] = un
            return vel
        xyz = self.lcoords
        lon = np.arctan2(xyz[:, 1], xyz[:, 0])
        lat = np.arcsin(np.clip(xyz[:, 2] / np.maximum(np.linalg.norm(xyz, axis=1), 1.0),
                                -1.0, 1.0))
        east = np.column_stack([-np.sin(lon), np.cos(lon), np.zeros_like(lon)])
        north = np.column_stack([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon),
                                 np.cos(lat)])
        return ue[:, None] * east + un[:, None] * north

    def _waveInflow(self, vel, grazing=0.05):
        """Domain-edge nodes (flat mesh) where the waves enter: the local
        travel direction points into the domain, ``v . n < -grazing`` with
        ``n`` the outward edge normal estimated from the neighbours (the test
        of ``_oroInflowBorders``, with a direction per node). Rank-local."""
        b = np.asarray(self.advectBorders, dtype=np.int64)
        if b.size == 0:
            return b
        ngb = self.FVmesh_ngbID[b]
        valid = ngb >= 0
        xy = self.lcoords[:, :2]
        cnt = np.maximum(valid.sum(axis=1), 1)
        cx = np.where(valid, xy[np.maximum(ngb, 0), 0], 0.0).sum(axis=1) / cnt
        cy = np.where(valid, xy[np.maximum(ngb, 0), 1], 0.0).sum(axis=1) / cnt
        nx, ny = xy[b, 0] - cx, xy[b, 1] - cy
        nrm = np.hypot(nx, ny)
        cosv = np.divide(vel[b, 0] * nx + vel[b, 1] * ny, nrm, out=np.zeros_like(nrm),
                         where=nrm > 0.0)
        return b[cosv < -grazing]

    def _waveCoefficients(self):
        """Upwind advection coefficients (unit speed, so the operator is
        ``d . grad`` in 1/m) and the inflow edges (flat meshes) of every
        direction of the climate. Computed once (the mesh does not change);
        ``getfacevelocity`` sets module state of the Fortran kernels, so the
        coefficients are copied at once."""
        self._waveCoef = []
        for frm, _ in self.waveClimate:
            bearing = self._waveBearings(frm)
            vel = self._waveVelocity(bearing)
            getfacevelocity(self.lpoints, vel)
            lcoeff = advecupwind(self.lpoints, 1.0)
            inflow = self._waveInflow(vel) if self.flatModel else np.zeros(0, dtype=np.int64)
            self._waveCoef.append((lcoeff[:, 0] - 1.0, lcoeff[:, 1:].copy(), inflow))

    def _waveSolve(self, diag, off, fixed, rhs):
        """Assemble ``diag`` + ``off`` with identity rows at ``fixed`` and
        solve with the right-hand side ``rhs`` (the imposed values on the
        fixed rows, the source elsewhere); returns the local solution.
        Collective."""
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
        self.tmpL.setArray(rhs)
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
        ``(d . grad + 1/L_f + kappa) E = 1/L_f`` (see the class). Collective. Scratch:
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
        grow = 0.0 if self.waveFetch is None else 1.0 / self.waveFetch
        E = np.zeros(self.lpoints, dtype=np.float64)
        for (frm, wgt), (adiag, off, inflow) in zip(self.waveClimate, self._waveCoef):
            fixed = land.copy()
            # the growth source over open water, and the open-ocean energy on
            # the inflow edges (flat meshes)
            val = np.where(land, 0.0, grow)
            inn = np.zeros(self.lpoints, dtype=bool)
            inn[inflow] = True
            inn &= ~land
            fixed |= inn
            val[inn] = 1.0
            Ek = self._waveSolve(adiag + kappa + grow, off, fixed, val)
            Ek[land] = 0.0
            E += wgt * np.clip(Ek, 0.0, 1.0)
        self.waveE = E
        if MPIrank == 0 and self.verbose:
            print("[waves] wave energy (%0.02f seconds)" % (process_time() - t0),
                  flush=True)
        return E
