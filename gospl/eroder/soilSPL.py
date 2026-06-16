import os
import gc
import sys
import petsc4py
import numpy as np
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time

from gospl.tools.constants import BEDROCK_EXPOSED

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import fctcoeff
    from gospl._fortran import local_spl
    from gospl._fortran import jacobiancoeff

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()


class soilSPL(object):
    """
    The class computes river incision expressed using a **stream power formulation** function of river discharge and slope also **accounting for soil production**.

    A non-linear diffusion of soil based on soil thickness is also implemented in this class.

    If the user has turned-on the sedimentation capability, this class will solve implicitly the **stream power formulation** accounting for a sediment transport/deposition term (`Yuan et al, 2019 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JF004867>`_).
    """

    def __init__(self, *args, **kwargs):
        """
        Initialisation of `soilSPL` class.
        """

        # Soil-SPL SNES controls. These are normally set from the YAML soil
        # block by the input parser; fall back to robust defaults when soilSPL
        # is built standalone (e.g. unit tests).
        self.soil_rtol = getattr(self, "soil_rtol", 1.0e-6)
        self.soil_atol = getattr(self, "soil_atol", 1.0e-6)
        self.soil_maxit = getattr(self, "soil_maxit", 500)
        self.soil_pc = getattr(self, "soil_pc", "hypre")

        # Cached SNES + helper vectors for the soil-aware SPL solver; created
        # lazily on first call and reused across timesteps.
        self._snes_soil = None
        self._snes_soil_f = None
        self._snes_soil_x = None

        # Robustness fallback SNES (L-BFGS quasi-Newton); created lazily and
        # only used on timesteps where the primary solver fails to converge.
        self._snes_soil_fb = None
        self._snes_soil_fb_f = None

        # Cached TS for the non-linear soil diffusion (rosw scheme); created
        # lazily on first call and reused across timesteps.
        self._ts_soil = None
        self._ts_soil_x = None
        self._ts_soil_f = None

        if self.Sperc > 0:
            self.soil_transition = -np.log(self.Sperc) * self.Hs
        else:
            self.soil_transition = 100.0
        self.Gsoil = self.hGlobal.duplicate()
        self.Lsoil = self.hLocal.duplicate()
        self.lHbed = self.hLocal.duplicate()
        self.gHbed = self.hGlobal.duplicate()

        # Allocate initial soil thicknesses
        if self.soilFile is not None:
            loadData = np.load(self.soilFile)
            soilH = loadData[self.soilData]
            self.Lsoil.setArray(soilH[self.locIDs])
            self.dm.localToGlobal(self.Lsoil, self.Gsoil)
        else:
            # Create an uniform soil thickness distribution
            self.Gsoil.set(self.cstSoilH)
            self.Lsoil.set(self.cstSoilH)

        # If temperatures dataset is provided then compute the corresponding soil production rate
        if self.tempFile is not None:
            Tref = self.tempRef + 273.15
            loadData = np.load(self.tempFile)
            # Subset the full-mesh temperature map to this rank's local nodes
            # (mirrors the soilFile branch above). Without [self.locIDs] the
            # derived prodSoil stays global (mpoints) and broadcasts against
            # the local hSoil/rainVal arrays only when MPIsize==1; in parallel
            # (lpoints < mpoints) it raises a ValueError shape mismatch.
            T = loadData[self.tempData][self.locIDs] + 273.15  # Conversion to Kelvin
            # Compute Arrhenius term, including Ea / R T0 term
            R = 8.314  # Gas constant (J/mol/K)
            Arr_terms = self.energyAct * (1./Tref - 1./T) / R
            self.prodSoil = self.P0 * np.exp(Arr_terms)
        else:
            self.prodSoil = self.P0 * np.ones(len(self.locIDs))

        return

    def _form_residual_soil(self, snes, h, F):
        """
        The nonlinear system (SNES) at each time step is solved iteratively by assessing the residual of the SPL equation accounting for erosion, deposition (transport-limited) and soil production.

        Parameters:
        -----------
        snes : PETSc.SNES: The snes object.
        h : PETSc.Vec: The current solution vector (h^{n+1}) at the new time step.
        F : PETSc.Vec: The residual vector to be filled.
        """

        # Current state
        self.dm.globalToLocal(h, self.hl)
        h_array = self.hl.getArray()

        # Compute slope
        S = local_spl(self.flowDir, h_array, self.rcvIDi,
                      self.distRcvi, self.wghtVali)
        S[S < 0.] = 0.

        # Compute upstream sediment flux
        if self.fDepa > 0:
            self.tmp.waxpy(-1.0, h, self.hOld)
            self.tmp1.pointwiseMult(self.tmp, self.areaGlobal)
            self.tmp1.scale(1. / self.dt)
            self._solve_KSP(True, self.fMati, self.tmp1, self.tmp)
            self.dm.globalToLocal(self.tmp, self.tmpL)
            Qt = self.tmpL.getArray()
            Qt[Qt < 0.] = 0.
            Qt[self.seaID] = 0.

        # Compute soil thickness based on changes in elevation and soil production rates
        hSoil = self.soilH + h_array - self.hOldArray
        hSoil[hSoil < 0] = 0.
        if self.tempFile is not None:
            # In this case, the soil production needs to be scaled with local rainfall as the soil production parameter refers here to the precipitation factor a0 in Norton et al. (2013) EQ. 8
            # http://dx.doi.org/10.1016/j.geomorph.2013.08.030
            hSoil += self.dt * self.prodSoil * self.rainVal * np.exp(-hSoil / self.Hs)
        else:
            hSoil += self.dt * self.prodSoil * np.exp(-hSoil / self.Hs)
        self.nsoilH = hSoil.copy()

        # Residuals based on the equation
        # h(t+dt) (1-G) - h(t) (1-G) + dt * K * A^m * S^n - dt * G * Qt / Area = 0
        res = (h_array - self.hOldArray) * (1.0 - self.fDep)
        res += self.Kbr * np.exp(-hSoil / self.h_star) * S**self.spl_n
        res += self.K_soil * (1.0 - np.exp(-hSoil / self.h_star)) * S**self.spl_n
        if self.fDepa > 0:
            res -= self.fDep * self.dt * Qt / self.larea

        # Residual vector
        F.setArray(res[self.glIDs])

        return

    def _monitorsoil(self, snes, its, norm):
        """
        Non-linear SPL with soil production solver convergence evaluation.
        """

        if MPIrank == 0 and its % 10 == 0:
            print(f"  ---  Non-linear soil SPL solver iteration {its}, Residual norm: {norm}", flush=True)

    def _build_soil_snes(self, primary=True):
        """
        Construct (and return) a soil-SPL SNES solver and its residual vector.

        Two configurations are available:

        * **primary** (``primary=True``) -- a Nonlinear GMRES accelerator
          (``ngmres``) right-preconditioned by Nonlinear Richardson
          (``nrichardson``). The Richardson sweep is what actually applies the
          Krylov solve (``cg``) and the multigrid preconditioner (HYPRE
          BoomerAMG by default, or ``self.soil_pc``); a *bare* ``ngmres`` ignores
          the KSP/PC entirely, which is why the previous setup stalled to
          ``SNES_DIVERGED_MAX_IT`` on the stiff soil-production residual.
        * **fallback** (``primary=False``) -- a limited-memory quasi-Newton
          solver (``qn``, L-BFGS) with a matrix-free critical-point line search.
          It builds an approximate Jacobian from secant updates (no analytic
          Jacobian required) and is markedly more robust for stiff problems; it
          runs only when the primary solver fails to converge.

        Each SNES gets its own PETSc options prefix so per-solver options (e.g.
        the line-search type) do not leak into the model's other SNES/KSP
        objects.

        :arg primary: select the primary (True) or fallback (False) solver.

        :return: the configured ``(SNES, residual Vec)`` pair.
        """

        opts = petsc4py.PETSc.Options()
        snes = petsc4py.PETSc.SNES().create(comm=petsc4py.PETSc.COMM_WORLD)
        f = self.hGlobal.duplicate()
        snes.setFunction(self._form_residual_soil, f)
        if self.verbose:
            snes.setMonitor(self._monitorsoil)

        if primary:
            snes.setOptionsPrefix("soilspl_")
            snes.setTolerances(rtol=self.soil_rtol, atol=self.soil_atol,
                               max_it=self.soil_maxit)
            snes.setType("ngmres")
            # Nonlinear right-preconditioner: one Richardson sweep per outer
            # iteration is what engages the Krylov solve + multigrid PC.
            npc = snes.getNPC()
            npc.setType("nrichardson")
            npc.setTolerances(max_it=1)
            ksp = npc.getKSP()
            ksp.setType("cg")
            pc = ksp.getPC()
            pc.setType(self.soil_pc)
            if self.soil_pc == "hypre":
                opts["pc_hypre_type"] = "boomeramg"
                pc.setFromOptions()
            ksp.setTolerances(rtol=1.0e-6)
        else:
            snes.setOptionsPrefix("soilsplfb_")
            # Relax the relative tolerance and grant a larger iteration budget:
            # the fallback only runs when the primary has already stalled.
            snes.setTolerances(rtol=max(self.soil_rtol * 100.0, 1.0e-4),
                               atol=self.soil_atol,
                               max_it=max(2 * self.soil_maxit, 200))
            snes.setType("qn")
            opts["soilsplfb_snes_qn_type"] = "lbfgs"
            # Critical-point line search: robust and matrix-free (the 'bt'
            # backtracking search requires a Jacobian, which qn does not form).
            opts["soilsplfb_snes_linesearch_type"] = "cp"

        snes.setFromOptions()
        return snes, f

    def _solveSoil(self):
        """
        Solves the non-linear stream power law for the transport limited and soil case. This calls the following *private function*:

        - _form_residual_soil

        .. note::

            PETSc SNES approach is used to solve the nonlinear equation without forming an analytic Jacobian. The primary solver is a Nonlinear GMRES accelerator (``ngmres``) right-preconditioned by Nonlinear Richardson (``nrichardson``), whose Krylov solve uses a Preconditioned Conjugate Gradient (``cg``) method with a multi-grid preconditioner (HYPRE BoomerAMG by default). If the primary solver stalls on the stiff soil-production residual, a limited-memory quasi-Newton fallback (``qn``, L-BFGS) with a critical-point line search is used. The iteration budget, tolerances and preconditioner are configurable through the YAML ``soil`` block (``maxIter``, ``rtol``, ``atol``, ``pcType``).
        """

        self.hOldArray = self.hLocal.getArray().copy()

        # Get soil thickness from previous time step
        self.soilH = self.Lsoil.getArray().copy()
        # Consider bedrock exposed when soil thickness is below 10 cm
        self.soilH[self.soilH < BEDROCK_EXPOSED] = 0.0

        # Upstream-averaged mean annual precipitation rate based on drainage area
        PA = self.FAL.getArray()

        # Define erosion limiter to prevent formation of flat
        dh = (self.hOldArray[:, None] - self.hOldArray[self.rcvIDi]).max(axis=1)
        elimiter = np.divide(dh, dh + 1.0e-2, out=np.zeros_like(dh),
                             where=dh != 0)

        # Per-node erodibility multiplier from the top of the local
        # stratigraphic column (1.0 = use self.K as-is). Only scales the
        # *bedrock* SPL coefficient; the soil-layer K is governed by
        # `self.Ksoil` and is left unchanged.
        # Fold in the dual-lithology erodibility blend (1.0 everywhere when
        # single-fraction, so behaviour is unchanged): K_eff = K*surfK*litK.
        surfK = self._surfaceK() * self._surfaceLithoK()

        # Incorporate the effect of local mean annual precipitation rate on erodibility (for soil and bedrock)
        if self.sedfacVal is not None:
            self.Kbr = self.K * surfK * self.sedfacVal * (self.rainVal ** self.coeffd)
        else:
            self.Kbr = self.K * surfK * (self.rainVal ** self.coeffd)
        self.Kbr *= self.dt * (PA ** self.spl_m) * elimiter
        self.Kbr[self.seaID] = 0.0

        self.K_soil = self.Ksoil * self.dt * (PA ** self.spl_m) * elimiter
        self.K_soil[self.seaID] = 0.0

        # Dimensionless depositional coefficient
        self.fDep = np.divide(self.fDepa * self.larea, PA, out=np.zeros_like(PA), where=PA != 0)
        self.fDep[self.seaID] = 0.
        self.fDep[self.fDep > 0.99] = 0.99
        if self.flatModel:
            self.fDep[self.idBorders] = 0.

        if self._snes_soil is None:
            self._snes_soil, self._snes_soil_f = self._build_soil_snes(primary=True)
            self._snes_soil_x = self.hGlobal.duplicate()

        snes = self._snes_soil
        x = self._snes_soil_x
        self.hGlobal.copy(result=x)
        snes.solve(None, x)
        r = snes.getConvergedReason()

        # Robustness net (mirrors the flow KSP's primary -> fallback path): if
        # the accelerated fixed-point solver stalls (typically
        # SNES_DIVERGED_MAX_IT on the stiff soil-production residual), retry
        # from the same initial guess with the limited-memory quasi-Newton
        # fallback. It only runs on the timesteps where the primary failed.
        if r < 0:
            if self._snes_soil_fb is None:
                self._snes_soil_fb, self._snes_soil_fb_f = self._build_soil_snes(
                    primary=False
                )
            fb = self._snes_soil_fb
            r0, it0 = r, snes.getIterationNumber()
            self.hGlobal.copy(result=x)
            fb.solve(None, x)
            r = fb.getConvergedReason()
            if MPIrank == 0:
                if r >= 0:
                    print(
                        "Soil SPL: primary SNES stalled (reason %d after %d its); "
                        "quasi-Newton fallback converged (reason %d, %d its)."
                        % (r0, it0, r, fb.getIterationNumber()),
                        flush=True,
                    )
                else:
                    print(
                        "Soil SPL SNES failed to converge: primary (reason %d) and "
                        "quasi-Newton fallback (reason %d) both diverged; continuing "
                        "with the best available iterate." % (r0, r),
                        flush=True,
                    )

        # Get eroded sediment thicknesses
        self.tmp.waxpy(-1.0, self.hOld, x)

        # Update soil thicknesses
        nHsoil = self.nsoilH.copy()
        # No subaerial soil production underwater: the soil-production term
        # scales with rainfall (Norton et al. 2013), so without this mask the
        # submarine nodes accumulate a spurious rainfall-scaled soil cover.
        nHsoil[self.seaID] = 0.0
        nHsoil[nHsoil < BEDROCK_EXPOSED] = 0.
        # Limit soil thickness
        nHsoil[nHsoil > self.soil_transition] = self.soil_transition
        self.Lsoil.setArray(nHsoil)
        self.dm.localToGlobal(self.Lsoil, self.Gsoil)

        petsc4py.PETSc.garbage_cleanup()

        return

    def _getEroDepRateSoil(self):
        """
        This function computes erosion deposition rates in metres per year and associated soil evolution. This is done on the filled elevation.

        The approach is based on **BasicHySa** governing equations from Terrainbento (as described in Appendix B20 from `Barnhart et al. (2019) <https://gmd.copernicus.org/articles/12/1267/2019/gmd-12-1267-2019.pdf>`_).

        .. note::

           The approach uses a continuous layer of soil-alluvium, which influences both hillslope and river-induced erosion. It relies on the SPACE algorithm of `Shobe et al. (2017) <https://gmd.copernicus.org/articles/10/4577/2017/>`_.
        """

        t0 = process_time()

        # Build the SPL erosion arrays
        if self.flexOn:
            self.hLocal.copy(result=self.hOldFlex)

        self._solveSoil()

        # Update erosion/deposition rate (thickness convention: positive
        # for deposition, negative for incision; same sign as cumED and
        # the on-disk EDrate field). See SPL.py for the convention note.
        E = self.tmp.getArray().copy()
        E = np.divide(E, self.dt)
        self.Eb.setArray(E)
        self.dm.globalToLocal(self.Eb, self.EbLocal)
        E = self.EbLocal.getArray().copy()
        if self.flatModel:
            E[self.idBorders] = 0.0
        E[self.lsink] = 0.0
        self.EbLocal.setArray(E)
        self.dm.localToGlobal(self.EbLocal, self.Eb)

        if MPIrank == 0 and self.verbose:
            print(
                "Finalise erosion deposition rates (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def updateSoilThickness(self):
        """
        Updates soil thickness through time.
        """

        self.dm.globalToLocal(self.tmp, self.tmpL)
        nHsoil = self.Lsoil.getArray() + self.tmpL.getArray()

        # Limit soil thickness
        nHsoil[nHsoil < 0.] = 0.
        nHsoil[nHsoil > self.soil_transition] = self.soil_transition

        self.Lsoil.setArray(nHsoil)
        self.dm.localToGlobal(self.Lsoil, self.Gsoil)

        return

    def erodepSPLsoil(self):
        """
        Modified **stream power law** model used to represent erosion by rivers also taking into account the role played by sediments in modulating erosion and deposition rate, considering **non-linear slope dependency** and accounting for soil production.

        It calls the private function `_getEroDepRateSoil` described above. Once erosion/deposition rates have been calculated, the function computes local thicknesses and soil evolution for the considered time step and update local elevation and cumulative erosion, deposition values.
        """

        t0 = process_time()

        # Computes the erosion deposition rates based on flow accumulation
        self.Eb.set(0.0)
        self.hGlobal.copy(result=self.hOld)
        self.dm.globalToLocal(self.hOld, self.hOldLocal)
        self._getEroDepRateSoil()
        self._glacialAbrasion()

        # Get erosion / deposition thicknesses (Eb is in thickness rate
        # convention: positive deposition, negative incision). See SPL.py.
        Eb = self.Eb.getArray().copy()
        self.tmp.setArray(Eb * self.dt)
        self.cumED.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)
        self.tmp1.pointwiseMult(self.tmp, self.areaGlobal)

        # Update stratigraphic layers
        if self.stratNb > 0:
            self.erodeStrat()
            self.deposeStrat()

        # Update erosion/deposition rates
        self.dm.globalToLocal(self.tmp, self.tmpL)
        add_rate = self.tmpL.getArray() / self.dt
        self.EbLocal.setArray(add_rate)

        # Destroy flow matrices
        self.fMati.destroy()
        self.fMat.destroy()

        if MPIrank == 0 and self.verbose:
            print(
                "Get Erosion Deposition values (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        if self.memclear:
            del Eb
            gc.collect()

        return

    def _evalFunctionSoil(self, ts, t, x, xdot, f):
        """
        The non-linear system for soil diffusion is solved iteratively using PETSc time stepping and SNES solution and is based on Rosenbrock W-scheme (``rosw``).

        Here again, we evaluate the residual function on a DMPlex for an implicit time-stepping method.

        Parameters:
        -----------
        ts : PETSc.TS: The time-stepper object.
        t : float: The current time.
        x : PETSc.Vec: The current solution vector (h^{n+1}) at the new time step.
        xdot : PETSc.Vec: The time derivative approximation (h^{n+1} - h^n) / dt.
        f : PETSc.Vec: The residual vector to be filled.
        """

        self.dm.globalToLocal(x, self.hl)
        with self.hl as hl, self.lHbed as zb, xdot as hdot:
            dh = hl - zb
            dh[dh < 0.1] = 0.0
            Cd = self.minDiff + np.multiply(self.Cd, (1.0 - np.exp(-dh / self.H0)))
            nlvec = fctcoeff(hl, Cd)
            f.setArray(hdot + nlvec[self.glIDs])

        return

    def _evalJacobianSoil(self, ts, t, x, xdot, a, A, B):
        """
        The non-linear system for soil diffusion is solved iteratively using PETSc time stepping and SNES solution and is based on Rosenbrock W-scheme (``rosw``).

        Here, we define the Jacobian matrix A and the preconditioner matrix B on a DMPlex.

        Parameters:
        -----------
        ts : PETSc.TS: The time-stepper object.
        t : float: The current time.
        x : PETSc.Vec: The current solution vector (h^{n+1}) at the new time step.
        xdot : PETSc.Vec: The time derivative approximation (h^{n+1} - h^n) / dt.
        a : float:  The shift factor for implicit methods.
        A : PETSc.Mat: The Jacobian matrix to be filled.
        B : PETSc.Mat: The preconditioner matrix to be filled.
        """

        self.dm.globalToLocal(x, self.hl)

        with self.hl as hl, self.lHbed as zb:
            dh = hl - zb
            dh[dh < 0.1] = 0.0
            Cd = self.minDiff + np.multiply(self.Cd, (1.0 - np.exp(-dh / self.H0)))

            # Coefficient derivatives
            Cp = np.multiply(self.Cd, np.exp(-dh / self.H0) / self.H0)
            nlC = jacobiancoeff(hl, Cd, Cp)

            # Combine the diagonal and off-diagonal columns into one
            # setValuesLocal call per row (was two), halving the Python -> PETSc
            # transitions. CSR-style batch is not safe here because the mesh
            # lgmap is keyed by mesh-local indices (length lpoints), not by
            # matrix-local rows (length n_owned).
            diag_col = np.arange(self.lpoints, dtype=petsc4py.PETSc.IntType)
            cols_2d = np.column_stack([diag_col[:, None], self.FVmesh_ngbID]).astype(
                petsc4py.PETSc.IntType
            )
            vals_2d = np.column_stack([(a + nlC[:, 0])[:, None], nlC[:, 1:]])
            for row in range(self.lpoints):
                B.setValuesLocal(row, cols_2d[row], vals_2d[row])
            B.assemble()

            if A != B:
                A.assemble()

        return True

    def _evalSolutionSoil(self, t, x):
        """
        Evaluate the initial solution of the SNES system.
        """

        assert t == 0.0, "only for t=0.0"
        x.setArray(self.h.getArray())

        return

    def diffuseSoil(self):
        r"""
        For river-transported sediments reaching the marine realm, this function computes the related marine deposition diffusion. It is based on a non-linear diffusion approach.

        .. math::
          \frac{\partial h}{\partial t}= \nabla \cdot \left( C_d \times (1.0 - e^{-h_s/H_0} \nabla h \right)

        It calls the following *private functions*:

        - _evalFunctionSoil
        - _evalJacobianSoil
        - _evalSolutionSoil

        .. note::

            PETSc SNES and time stepping TS approaches are used to solve the non-linear equation above over the considered time step.
        """

        t0 = process_time()

        # Get diffusion soil coefficient
        self.Cd = np.full(self.lpoints, self.Cda, dtype=np.float64)
        self.Cd[self.seaID] = self.Cdm
        # Dual-lithology (Phase 7): scale soil diffusivity by the surface
        # composition so fine-rich soil diffuses faster (neutral when single-
        # fraction / no contrast).
        if self.stratLith:
            self.Cd = self.Cd * self._surfaceLithoD()

        # Remove the soil thickness from the elevation
        self.hLocal.copy(result=self.hl)
        self.dm.localToGlobal(self.hl, self.h)
        self.gHbed.waxpy(-1.0, self.Gsoil, self.hGlobal)
        self.lHbed.waxpy(-1.0, self.Lsoil, self.hLocal)

        # Time stepping definition (cached across timesteps)
        if self._ts_soil is None:
            ts = petsc4py.PETSc.TS().create(comm=petsc4py.PETSc.COMM_WORLD)
            # arkimex: IMEX Runge-Kutta schemes | rosw: Rosenbrock W-schemes
            ts.setType("rosw")
            ts.setIFunction(self._evalFunctionSoil, self.tmp1)
            ts.setIJacobian(self._evalJacobianSoil, self.mat)
            ts.setExactFinalTime(petsc4py.PETSc.TS.ExactFinalTime.MATCHSTEP)
            # Allow an unlimited number of failures (step rejected and retried)
            ts.setMaxSNESFailures(-1)
            # SNES nonlinear solver
            snes = ts.getSNES()
            snes.setTolerances(max_it=10)
            # KSP linear solver
            ksp = snes.getKSP()
            ksp.setType("preonly")
            pc = ksp.getPC()
            pc.setType("gasm")
            ts.setFromOptions()
            self._ts_soil = ts
            self._ts_soil_x = self.tmp1.duplicate()
            self._ts_soil_f = self.tmp1.duplicate()

        ts = self._ts_soil
        x = self._ts_soil_x
        # Soil thicknesses are meters; mm-level absolute tolerance is plenty.
        ts.setTolerances(atol=1e-3, rtol=1e-3)
        ts.setTime(0.0)
        # Larger initial step (was self.dt / 1000.0).
        ts.setTimeStep(self.dt / 100.0)
        ts.setMaxTime(self.dt)
        ts.setMaxSteps(self.tsStep)

        tstart = ts.getTime()
        self._evalSolutionSoil(tstart, x)

        # Solve nonlinear equation
        ts.solve(x)

        if MPIrank == 0 and self.verbose:
            print(
                "Nonlinear soil diffusion solution (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )
            print(
                "steps %d (%d rejected, %d SNES fails), nonlinear its %d, linear its %d"
                % (
                    ts.getStepNumber(),
                    ts.getStepRejections(),
                    ts.getSNESFailures(),
                    ts.getSNESIterations(),
                    ts.getKSPIterations(),
                ),
                flush=True,
            )

        # Get diffused sediment thicknesses
        self.dh.waxpy(-1.0, self.hGlobal, x)
        self.dm.globalToLocal(self.dh, self.tmpL)
        chgSoil = self.tmpL.getArray().copy()
        self.tmpL.setArray(chgSoil)
        self.dm.localToGlobal(self.tmpL, self.tmp)

        petsc4py.PETSc.garbage_cleanup()

        # Update cumulative erosion and deposition as well as elevation
        self.cumED.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        # Update soil thickness
        self.updateSoilThickness()

        # Update erosion/deposition rates
        self.dm.globalToLocal(self.tmp, self.tmpL)
        add_rate = self.tmpL.getArray() / self.dt
        self.tmpL.setArray(add_rate)
        self.EbLocal.axpy(1.0, self.tmpL)

        # Update stratigraphic layer parameters
        if self.stratNb > 0:
            self.deposeStrat()

        if MPIrank == 0 and self.verbose:
            print(
                "Diffuse Soil Sediments (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        return
