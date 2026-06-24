import os
import gc
import sys
import petsc4py
from gospl.tools.petscgc import safe_garbage_cleanup
import numpy as np

from mpi4py import MPI
from time import process_time

from gospl.tools.constants import MARINE_SMOOTH_N_LAND, MARINE_SMOOTH_N_SEA

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import sethillslopecoeff
    from gospl._fortran import hillslp_nl
    from gospl._fortran import jacobiancoeff
    from gospl._fortran import fctcoeff

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()


class hillSLP(object):
    """
    This class encapsulates all the functions related to hillslope (soil creep) processes (both linear and non-linear cases are implemented).
    """

    def __init__(self, *args, **kwargs):
        """
        The initialisation of `hillSLP` class consists in the declaration of several PETSc vectors.
        """

        self.dlim = False
        self.Dlimit = 5.
        self.dexp = 0.05
        self.minDiff = 1.e-4

        self.mat = self.dm.createMatrix()
        self.mat.setOption(self.mat.Option.NEW_NONZERO_LOCATIONS, True)

        self.h = self.hGlobal.duplicate()
        self.hl = self.hLocal.duplicate()
        self.dh = self.hGlobal.duplicate()

        # Cached SNES + helper vectors for the non-linear hillslope solver;
        # created lazily on first call and reused across timesteps.
        self._snes_hill = None
        self._snes_hill_f = None
        self._snes_hill_x = None

        # Cached TS for marine non-linear diffusion (rosw scheme); created
        # lazily on first call and reused across timesteps.
        self._ts_marine = None
        self._ts_marine_x = None

        # Lagged-diffusivity (Picard) alternative to the marine/lake diffusion
        # TS (opt-in via self.marineSolver == 'picard'); cached linear KSP.
        self._ksp_picard = None

        # Cached marine-morphology smoothing operator (_hillSlope smooth=2) and
        # its KSP. The smoothing diffusivity is mesh-size-based and timestep-
        # independent, so the operator depends only on the (fixed) mesh and the
        # land/sea mask; it is built once and reused, and rebuilt only when the
        # coastline (seaID) moves. _smooth_pc_ready tracks whether the cached
        # preconditioner is valid for the current operator.
        self._smoothMat = None
        self._ksp_smooth = None
        self._smooth_seaID = None
        self._smooth_pc_ready = False

        # Cached LINEAR-hillslope operator (_hillSlope smooth=0) + its KSP. Same
        # invariance as the marine smoother: with constant Cda/Cdm (no dual
        # lithology) the operator (I - Cd*dt*L) depends only on the fixed mesh
        # and the land/sea mask, so it is built once, the PC reused, and rebuilt
        # only when the coastline moves. NOT used under dual lithology (Cd is
        # then step-varying via _surfaceLithoD).
        self._hillMat = None
        self._ksp_hill_lin = None
        self._hill_seaID = None
        self._hill_pc_ready = False

        return

    def _hillSlope(self, smooth=0):
        r"""
        This function computes hillslope using a **linear** diffusion law commonly referred to as **soil creep**:

        .. math::
          \frac{\partial z}{\partial t}= \kappa_{D} \nabla^2 z

        in which :math:`\kappa_{D}` is the diffusion coefficient and can be defined with different values for the marine and land environments (set with `hillslopeKa` and `hillslopeKm` in the YAML input file). It encapsulates, in a simple formulation, processes operating on superficial sedimentary layers. Main controls on variations of :math:`\kappa_{D}` include substrate, lithology, soil depth, climate and biological activity.

        .. note::
            The hillslope processes in goSPL are considered to be happening at the same rate for coarse and fine sediment sizes.

        :arg smooth: smoothing mode. ``0`` (default) computes linear soil-creep
            hillslope and updates the elevation in place. ``1`` returns a one-off
            smoothed copy of the ice-surface proxy used to derive glacial flow
            directions. ``2`` returns a smoothed bathymetry used **only** to
            derive coherent marine flow directions in ``_matOcean`` — it never
            alters the elevation. The ``smooth=2`` operator uses a mesh-size-
            based, timestep-independent diffusivity (see
            ``MARINE_SMOOTH_N_*``) and is cached + reused across steps, rebuilt
            only when the coastline moves.

        .. note::
            This method uses scratch Vecs ``self.tmp`` / ``self.tmpL`` (and, for
            ``smooth=1``, ``self.tmp1`` as the input field). ``smooth=0`` also
            uses ``self.hOld``.
        """

        if smooth == 0:
            if self.Cda == 0.0 and self.Cdm == 0.0:
                return

        t0 = process_time()

        # --- smooth=2: marine-morphology smoothing (flow-direction finding) ---
        # The diffusivity is mesh-size-based: per-node Kd = N * cell_area, with N
        # a dimensionless smoothing number (MARINE_SMOOTH_N_*). Because the FV
        # Laplacian stencil scales as Kd / Δx^2, this is both timestep- and
        # resolution-independent. The operator then depends only on the fixed
        # mesh and the land/sea mask, so it is cached and its preconditioner
        # reused; it is rebuilt only when the coastline (seaID) moves.
        if smooth == 2:
            # The rebuild test is rank-local (seaID is per-partition), but the
            # rebuild + PCSetUp it gates are COLLECTIVE — so the decision MUST be
            # reduced across ranks (logical OR), or one rank rebuilds while
            # another reuses and the collective Mat assembly / PCSetUp deadlocks.
            rebuild = self._smoothMat is None or not np.array_equal(
                self.seaID, self._smooth_seaID
            )
            rebuild = MPI.COMM_WORLD.allreduce(rebuild, op=MPI.LOR)
            if rebuild:
                Cd = np.full(self.lpoints, MARINE_SMOOTH_N_LAND, dtype=np.float64)
                Cd[self.seaID] = MARINE_SMOOTH_N_SEA
                if self._smoothMat is not None:
                    self._smoothMat.destroy()
                self._smoothMat = self._buildDiffMat(Cd * self.larea)
                self._smooth_seaID = self.seaID.copy()
                self._smooth_pc_ready = False
            self._solveSmooth(self.hGlobal, self.tmp)
            self.dm.globalToLocal(self.tmp, self.tmpL)
            return self.tmpL.getArray().copy()

        # --- smooth=1: ice-surface smoothing (always rebuilt; ice path) ---
        # The diffusivity is the physical Ka/Km (optionally lithology-scaled).
        if smooth == 1:
            Cd = np.full(self.lpoints, self.Cda, dtype=np.float64)
            Cd[self.seaID] = self.Cdm
            # Dual-lithology (Phase 7): scale the diffusivity by the surface
            # composition so fines diffuse faster (1.0 everywhere when single-
            # fraction / no contrast, so behaviour is unchanged).
            if self.stratLith:
                Cd = Cd * self._surfaceLithoD()
            diffMat = self._buildDiffMat(Cd * self.dt)
            if self.tmp1.max()[1] > 0:
                self._solve_KSP(True, diffMat, self.tmp1, self.tmp)
            else:
                self.tmp1.copy(result=self.tmp)
            diffMat.destroy()
            self.dm.globalToLocal(self.tmp, self.tmpL)
            return self.tmpL.getArray().copy()

        # --- smooth=0: linear soil creep (updates elevation) ---
        self.hGlobal.copy(result=self.hOld)
        if self.stratLith:
            # Dual lithology makes Cd step-varying (surface composition evolves),
            # so the operator is rebuilt every step and uses the shared flow KSP.
            Cd = np.full(self.lpoints, self.Cda, dtype=np.float64)
            Cd[self.seaID] = self.Cdm
            Cd = Cd * self._surfaceLithoD()
            diffMat = self._buildDiffMat(Cd * self.dt)
            self._solve_KSP(True, diffMat, self.hOld, self.hGlobal)
            diffMat.destroy()
        else:
            # Constant Cda/Cdm ⇒ the operator (I - Cd*dt*L) is mesh+mask-
            # invariant: cache it and reuse the factorised PC, rebuilding only
            # when the coastline (seaID) moves. Mirrors the marine smoother
            # (smooth=2); bit-faithful to the previous per-step build+solve.
            # Collective rebuild decision (see smooth=2 above): seaID is rank-
            # local but _buildDiffMat + PCSetUp are collective, so reduce the
            # "coastline moved" test across ranks or the solve deadlocks at np>1.
            rebuild = self._hillMat is None or not np.array_equal(
                self.seaID, self._hill_seaID
            )
            rebuild = MPI.COMM_WORLD.allreduce(rebuild, op=MPI.LOR)
            if rebuild:
                Cd = np.full(self.lpoints, self.Cda, dtype=np.float64)
                Cd[self.seaID] = self.Cdm
                if self._hillMat is not None:
                    self._hillMat.destroy()
                self._hillMat = self._buildDiffMat(Cd * self.dt)
                self._hill_seaID = self.seaID.copy()
                self._hill_pc_ready = False
            if self._ksp_hill_lin is None:
                self._ksp_hill_lin = self._makeDiffusionKSP("hilllin_")
            ksp = self._ksp_hill_lin
            ksp.setOperators(self._hillMat, self._hillMat)
            ksp.getPC().setReusePreconditioner(self._hill_pc_ready)
            self._hill_pc_ready = True
            ksp.solve(self.hOld, self.hGlobal)
            safe_garbage_cleanup()

        # Update cumulative erosion/deposition and elevation
        self.tmp.waxpy(-1.0, self.hOld, self.hGlobal)
        self.cumED.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        if self.memclear:
            del Cd
            gc.collect()

        if self.stratNb > 0:
            self.erodeStrat()
            self.deposeStrat()

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Linear Hillslope Processes (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )
        safe_garbage_cleanup()

        return

    def _buildDiffMat(self, Kd):
        """
        Assemble the implicit linear-diffusion operator :math:`(I - L)` for the
        finite-volume Laplacian, given per-node coefficients ``Kd`` (m², i.e. a
        diffusivity × time, or the mesh-size-based marine-smoothing coefficient
        ``N · cell_area``).

        Shared by the linear hillslope (``smooth=0``), the ice-surface smoothing
        (``smooth=1``) and the cached marine-smoothing operator (``smooth=2``).

        :arg Kd: per-node diffusion coefficient array (``lpoints``, m²).

        :return: the assembled PETSc diffusion matrix (caller owns it; the
            ``smooth=0/1`` callers destroy it after the solve, the ``smooth=2``
            caller caches it).
        """

        diffCoeffs = sethillslopecoeff(self.lpoints, Kd)
        if self.flatModel:
            diffCoeffs[self.outletIDs, 1:] = 0.0
            diffCoeffs[self.outletIDs, 0] = 1.0

        return self._assembleDiffMat(diffCoeffs)

    def _assembleDiffMat(self, coeffs):
        """
        Assemble a PETSc finite-volume stencil matrix from a per-node
        coefficient array ``coeffs`` of shape ``(lpoints, 1 + maxnb)``: column 0
        is the diagonal, columns ``1..maxnb`` are the off-diagonal entries for
        the ``maxnb`` neighbours in ``self.FVmesh_ngbID``. Used by
        ``_buildDiffMat`` (``sethillslopecoeff`` output) and the lagged-
        diffusivity marine solver (``jacobiancoeff`` output, scaled to
        ``I + dt·L``).

        :arg coeffs: ``(lpoints, 1+maxnb)`` diagonal + neighbour coefficients.

        :return: the assembled PETSc matrix (caller owns / destroys it).
        """

        diffMat = self._matrix_build_diag(coeffs[:, 0])
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)

        for k in range(0, self.maxnb):
            tmpMat = self._matrix_build()
            indices = self.FVmesh_ngbID[:, k].copy()
            data = coeffs[:, k + 1]
            ids = np.nonzero(data == 0.0)
            indices[ids] = ids
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                indices.astype(petsc4py.PETSc.IntType),
                data,
            )
            tmpMat.assemblyEnd()
            diffMat.axpy(1.0, tmpMat)
            tmpMat.destroy()

        return diffMat

    def _assembleDiffMatCSR(self, coeffs):
        """
        Single-pass CSR assembly of the FV stencil matrix from ``coeffs``
        ``(lpoints, 1+maxnb)`` (col 0 = diagonal, cols ``1..maxnb`` = neighbour
        entries for ``self.FVmesh_ngbID``). Equivalent to ``_assembleDiffMat``
        (diag + per-neighbour axpy) but builds the matrix in ONE
        ``setValuesLocalCSR`` pass — ``ADD_VALUES`` sums shared columns exactly,
        mirroring ``seaplex._matOcean`` (#450). Much cheaper than the 12-axpy
        loop, used in the Picard hot loop (`_diffuseImplicitPicard`), which
        rebuilds the operator every iteration.

        ``_buildDiffMat`` deliberately keeps the 12-axpy ``_assembleDiffMat`` so
        the bit-faithful cached smoother / linear-hillslope operators (#457 /
        #459) are not perturbed by a different floating-point summation order.

        :arg coeffs: ``(lpoints, 1+maxnb)`` diagonal + neighbour coefficients.

        :return: the assembled PETSc matrix (caller owns / destroys it).
        """

        nnz = 1 + self.maxnb
        nodes = np.arange(self.lpoints, dtype=petsc4py.PETSc.IntType)
        ngb = self.FVmesh_ngbID[:, : self.maxnb].astype(petsc4py.PETSc.IntType).copy()
        offvals = coeffs[:, 1 : 1 + self.maxnb].copy()
        # Redirect zero-valued neighbour entries to the diagonal so the column
        # index is always valid (adds 0 under ADD_VALUES) — matches the
        # `indices[ids] = ids` trick in _assembleDiffMat.
        zero = offvals == 0.0
        ngb[zero] = np.broadcast_to(nodes[:, None], ngb.shape)[zero]

        cols = np.empty((self.lpoints, nnz), dtype=petsc4py.PETSc.IntType)
        cols[:, 0] = nodes
        cols[:, 1:] = ngb
        vals = np.empty((self.lpoints, nnz), dtype=np.float64)
        vals[:, 0] = coeffs[:, 0]
        vals[:, 1:] = offvals
        indptr = np.arange(
            0, self.lpoints * nnz + 1, nnz, dtype=petsc4py.PETSc.IntType
        )

        M = self._matrix_build(nnz=(nnz, nnz))
        M.assemblyBegin()
        M.setValuesLocalCSR(
            indptr,
            cols.ravel(),
            vals.ravel(),
            addv=petsc4py.PETSc.InsertMode.ADD_VALUES,
        )
        M.assemblyEnd()

        return M

    def _makeDiffusionKSP(self, prefix):
        """
        Create a cached KSP for a linear diffusion / smoothing solve: PETSc
        Richardson + block-Jacobi (matching ``flowplex._solve_KSP``), nonzero
        initial guess, with the per-rank ILU pivots shifted so the block-Jacobi
        ``PCSetUp`` cannot fail on a degenerate / ocean-only partition. The
        options prefix scopes that shift (and any user override) to this solver
        only. Shared by the cached marine smoother (``_solveSmooth``) and the
        cached linear hillslope (``_hillSlope`` smooth=0).

        :arg prefix: PETSc options prefix for this KSP (e.g. ``"marsmooth_"``).

        :return: the configured PETSc KSP (caller caches and owns it).
        """

        ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
        ksp.setType("richardson")
        ksp.getPC().setType("bjacobi")
        ksp.setTolerances(rtol=self.rtol)
        ksp.setInitialGuessNonzero(True)
        ksp.setOptionsPrefix(prefix)
        petsc4py.PETSc.Options()[f"{prefix}sub_pc_factor_shift_type"] = "nonzero"
        ksp.setFromOptions()

        return ksp

    def _solveSmooth(self, rhs, sol):
        """
        Solve one implicit step of the cached marine-smoothing diffusion
        operator (see ``_hillSlope`` ``smooth=2``).

        The operator (``self._smoothMat``) is rebuilt only when the coastline
        moves, so the block-Jacobi/ILU preconditioner is factorised once and
        **reused** across every step with an unchanged land/sea mask — that
        factorisation (not the 1-iteration Richardson solve) was the dominant
        per-step cost of the previous build-from-scratch approach.

        :arg rhs: PETSc global Vec — the field to smooth (read).
        :arg sol: PETSc global Vec — the smoothed result (written; scratch
            ``self.tmp`` at the call site).
        """

        if self._ksp_smooth is None:
            self._ksp_smooth = self._makeDiffusionKSP("marsmooth_")

        ksp = self._ksp_smooth
        ksp.setOperators(self._smoothMat, self._smoothMat)
        # Reuse the factorised preconditioner unless the operator was just
        # rebuilt (coastline moved) — then PCSetUp runs once for the new matrix.
        # (setReusePreconditioner lives on the PC, not the KSP, in petsc4py.)
        ksp.getPC().setReusePreconditioner(self._smooth_pc_ready)
        self._smooth_pc_ready = True
        ksp.solve(rhs, sol)
        safe_garbage_cleanup()

        return

    def _diff_nl_monitor(self, snes, its, norm):
        """
        Non-linear diffusion solver convergence evaluation.
        """

        if MPIrank == 0 and its % 5 == 0:
            print(f"  ---  Non-linear hillslope solver iteration {its}, Residual norm: {norm}", flush=True)

    def _form_residual_nl_hillslope(self, snes, h, F):
        """
        The nonlinear system (SNES) at each time step is solved iteratively by assessing the residual of the hillslope equation.

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
        if self.K_sc == 0:
            val = hillslp_nl(self.lpoints, h_array, self.Cd_nl, self.K_nl, -1)
        else:
            val = hillslp_nl(self.lpoints, h_array, self.Cd_nl, self.K_sc, self.K_nb)

        if self.flatModel:
            val[self.outletIDs] = 0.

        # Residuals
        res = h_array - self.hOldArray - self.dt * val

        # Residual vector
        F.setArray(res[self.glIDs])

        return

    def _hillSlopeNL(self):
        r"""
        This function computes hillslope using a **non-linear** diffusion law:

        .. math::
          \frac{\partial h}{\partial t}= \nabla \cdot \left( C_d(h) \nabla h \right)

        Two flavors of non-linear diffusion are possible based on user defined parameters:

        1. a non-critical hillslope model following the work of `Wang et al. (2024) <https://www.sciencedirect.com/science/article/pii/S0169555X24001053>`_.
        2. a non-linear depth-dependent creep law as described in `Barnhart et al. (2019) <https://gmd.copernicus.org/articles/12/1267/2019/gmd-12-1267-2019.pdf>`_ (section 3.4.5).

        .. note::
            In this implementation of the SNES, we do not form the Jacobian and PETSc calculates it based on the residual function. Here, a Nonlinear Generalized Minimum Residual method is used ``ngmres``, a Preconditioned Conjugate Gradient ``cg`` method is defined for the KSP and the preconditioner allows for multi-grid methods based on the HYPRE BoomerAMG approach.
        """

        t0 = process_time()

        self.hGlobal.copy(result=self.hOld)
        self.Cd_nl = np.full(self.lpoints, self.Cda, dtype=np.float64)
        self.Cd_nl[self.seaID] = self.Cdm
        # Dual-lithology (Phase 7): scale the non-linear diffusivity by the
        # surface composition so fines diffuse faster (neutral when single-
        # fraction / no contrast).
        if self.stratLith:
            self.Cd_nl = self.Cd_nl * self._surfaceLithoD()
        self.hOldArray = self.hLocal.getArray().copy()

        if self._snes_hill is None:
            snes = petsc4py.PETSc.SNES().create(comm=petsc4py.PETSc.COMM_WORLD)
            snes.setTolerances(rtol=self.snes_rtol, atol=self.snes_atol,
                               max_it=self.snes_maxit)
            if self.verbose:
                snes.setMonitor(self._diff_nl_monitor)
            self._snes_hill_f = self.hGlobal.duplicate()
            snes.setFunction(self._form_residual_nl_hillslope, self._snes_hill_f)
            snes.setType('ngmres')
            ksp = snes.getKSP()
            ksp.setType("cg")
            pc = ksp.getPC()
            pc.setType("hypre")
            petsc_options = petsc4py.PETSc.Options()
            petsc_options['pc_hypre_type'] = 'boomeramg'
            ksp.getPC().setFromOptions()
            ksp.setTolerances(rtol=1.0e-6)
            self._snes_hill_x = self.hGlobal.duplicate()
            self._snes_hill = snes

        snes = self._snes_hill
        x = self._snes_hill_x
        self.hGlobal.copy(result=x)
        snes.solve(None, x)
        r = snes.getConvergedReason()
        if r < 0 and MPIrank == 0:
            print(
                "Non-linear hillslope SNES failed to converge after %d iterations (reason %d)"
                % (snes.getIterationNumber(), r),
                flush=True,
            )

        x.copy(result=self.hGlobal)

        # Update cumulative erosion/deposition and elevation
        self.tmp.waxpy(-1.0, self.hOld, self.hGlobal)
        self.cumED.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        if self.stratNb > 0:
            self.erodeStrat()
            self.deposeStrat()

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Non-Linear Hillslope Processes (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )
        safe_garbage_cleanup()

        return

    def getHillslope(self):
        """
        This function is the main entry point to compute linear and non-linear hillslope processes. 
        """

        h = self.hLocal.getArray().copy()
        self.seaID = np.where(h <= self.sealevel)[0]

        # Specify the hillslope diffusion law to use
        if self.cptSoil:
            self.diffuseSoil()
            return

        if self.K_nl == 1.0 and self.K_sc == 0.0:
            self._hillSlope(smooth=0)
        else:
            self._hillSlopeNL()

        # Update layer elevation
        if self.stratNb > 0:
            self.elevStrat()
        if self.memclear:
            del h
            gc.collect()

        # Update erosion/deposition rates
        self.dm.globalToLocal(self.tmp, self.tmpL)
        add_rate = self.tmpL.getArray() / self.dt
        self.tmpL.setArray(add_rate)
        self.EbLocal.axpy(1.0, self.tmpL)

        if self.memclear:
            del add_rate
            gc.collect()

        return

    def _evalFunctionMardDiff(self, ts, t, x, xdot, f):
        """
        The non-linear system for freshly-deposited marine sediments diffusion at each time step is solved iteratively using PETSc time stepping and SNES solution and is based on Rosenbrock W-scheme (``rosw``).

        Here we evaluate the residual function on a DMPlex for an implicit time-stepping method.

        Parameters:
        -----------
        ts : PETSc.TS: The time-stepper object.
        t : float: The current time.
        x : PETSc.Vec: The current solution vector (h^{n+1}) at the new time step.
        xdot : PETSc.Vec: The time derivative approximation (h^{n+1} - h^n) / dt.
        f : PETSc.Vec: The residual vector to be filled.
        """

        self.dm.globalToLocal(x, self.hl)
        with self.hl as hl, self.hLocal as zb, xdot as hdot:
            dh = hl - zb
            dh[dh < 0.1] = 0.0
            if self.dlim:
                Cd = self.minDiff_vec + np.multiply(self.Cd, dh / (dh + self.Dlimit))
            else:
                Cd = self.minDiff_vec + np.multiply(self.Cd, (1.0 - np.exp(-self.dexp * dh)))
            nlvec = fctcoeff(hl, Cd)
            f.setArray(hdot + nlvec[self.glIDs])

        return

    def _evalJacobianMardDiff(self, ts, t, x, xdot, a, A, B):
        """
        The non-linear system for freshly-deposited marine sediments diffusion at each time step is solved iteratively using PETSc time stepping and SNES solution and is based on Rosenbrock W-scheme (``rosw``).

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

        with self.hl as hl, self.hLocal as zb:
            dh = hl - zb
            dh[dh < 0.1] = 0.0
            if self.dlim:
                Cd = self.minDiff_vec + np.multiply(self.Cd, dh / (dh + self.Dlimit))
            else:
                Cd = self.minDiff_vec + np.multiply(self.Cd, (1.0 - np.exp(-self.dexp * dh)))

            # Coefficient derivatives
            if self.dlim:
                Cp = np.multiply(self.Cd, self.Dlimit / (dh + self.Dlimit)**2)
            else:
                Cp = np.multiply(self.Cd, self.dexp * np.exp(-self.dexp * dh))
            nlC = jacobiancoeff(hl, Cd, Cp)
            diag = a + nlC[:, 0]
            offdiag = nlC[:, 1:]
            # Combine the diagonal and off-diagonal columns into one array for setting values in the matrix
            ngb_cols = self.FVmesh_ngbID[self.glIDs, :]
            cols_2d = np.column_stack([self.glIDs[:, None], ngb_cols]).astype(
                petsc4py.PETSc.IntType
            )
            vals_2d = np.column_stack([diag[self.glIDs, None], offdiag[self.glIDs]])
            for i, row in enumerate(self.glIDs):
                B.setValuesLocal(row, cols_2d[i], vals_2d[i])
            B.assemble()

            if A != B:
                A.assemble()

        return True

    def _evalSolutionMardDiff(self, t, x):
        """
        Evaluate the initial solution of the SNES system.
        """

        assert t == 0.0, "only for t=0.0"
        x.setArray(self.h.getArray())

        return

    def _diffuseOcean(self, dh):
        r"""
        For river-transported sediments reaching the marine realm, this function computes the related marine deposition diffusion. It is based on a non-linear diffusion approach.

        .. math::
          \frac{\partial h}{\partial t}= \nabla \cdot \left( C_d(h) \nabla h \right)

        It calls the following *private functions*:

        - _evalFunctionMardDiff
        - _evalJacobianMardDiff
        - _evalSolutionMardDiff

        .. note::

            PETSc SNES and time stepping TS approaches are used to solve the non-linear equation above over the considered time step.

        :arg dh: Numpy Array of incoming marine depositional thicknesses
        """

        t0 = process_time()
        ndepo = self._diffuseImplicit(
            dh, self.seaID, self.nlK, label="marine"
        )
        self.tmpL.setArray(ndepo)
        self.dm.localToGlobal(self.tmpL, self.tmp)

        if MPIrank == 0 and self.verbose:
            print(
                "Marine diffusion total (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        if self.memclear:
            del ndepo
            gc.collect()
        safe_garbage_cleanup()

        return

    def _diffuseImplicit(self, dh, mask, Cd_val, label="diffusion"):
        """
        Solve the non-linear thickness diffusion equation

        .. math::
            \\partial_t h = \\nabla \\cdot (C_d(h) \\nabla h)

        on the cells indicated by ``mask`` for one global timestep ``self.dt``.
        Cells outside ``mask`` get a zero diffusion coefficient and therefore
        do not move. The cached PETSc TS (rosw) solver is reused across calls;
        ``self.Cd`` and ``self.minDiff_vec`` are swapped in-place so the same
        solver can be applied to marine sediments, lake sediments, or any
        other sub-domain.

        :arg dh: per-node initial deposit thickness (m)
        :arg mask: bool array (lpoints) or integer index array; cells where
                   diffusion is active.
        :arg Cd_val: scalar non-linear diffusion coefficient (m^2/yr)
        :arg label: log label (e.g. "marine" or "lake") for the timing print

        :return: ndepo - smoothed deposit thickness (m), zeroed where negative
        """

        # Opt-in: lagged-diffusivity (Picard) linear solver instead of the
        # adaptive non-linear TS (see _diffuseImplicitPicard). Default 'ts'.
        if getattr(self, "marineSolver", "ts") == "picard":
            return self._diffuseImplicitPicard(dh, mask, Cd_val, label=label)

        t0 = process_time()

        self.mat.zeroEntries()

        self.Cd = np.zeros(self.lpoints)
        self.Cd[mask] = Cd_val

        self.minDiff_vec = np.zeros(self.lpoints)
        self.minDiff_vec[mask] = self.minDiff

        self.hl.setArray(dh)
        self.hl.axpy(1.0, self.hLocal)
        self.dm.localToGlobal(self.hl, self.h)

        # Cached TS (rosw IMEX) shared across all _diffuseImplicit callers.
        if self._ts_marine is None:
            ts = petsc4py.PETSc.TS().create(comm=petsc4py.PETSc.COMM_WORLD)
            ts.setType("rosw")
            ts.setIFunction(self._evalFunctionMardDiff, self.tmp1)
            ts.setIJacobian(self._evalJacobianMardDiff, self.mat)
            ts.setExactFinalTime(petsc4py.PETSc.TS.ExactFinalTime.MATCHSTEP)
            ts.setMaxSNESFailures(-1)
            snes = ts.getSNES()
            snes.setTolerances(max_it=10)
            ksp = snes.getKSP()
            ksp.setType("preonly")
            pc = ksp.getPC()
            pc.setType("gasm")
            ts.setFromOptions()
            self._ts_marine = ts
            self._ts_marine_x = self.tmp1.duplicate()

        ts = self._ts_marine
        x = self._ts_marine_x
        # Note: do NOT call ts.reset() here. TSReset destroys the TS's
        # internal Jacobian reference (the link we set up via setIJacobian
        # at TS creation), causing PCSetUp_GASM to fail with PETSc error 56
        # on the next solve. The TS internal state is fine to reuse across
        # calls; setTime / setTimeStep below override the previous timing.
        # Per-step error bound = min(atol, rtol*|x|). For marine sediment
        # thicknesses we already drop sub-mm deposits as numerical noise
        # elsewhere in the pipeline, so a 1 cm absolute floor is physically
        # adequate. rtol=1e-4 (0.01 %) is loose enough for the controller to
        # take large steps but tight enough to suppress the "peak" overshoot
        # artefacts seen with 1e-3/1e-3.
        ts.setTolerances(atol=5.0e-3, rtol=1.0e-4)
        ts.setTime(0.0)
        # Reset the step COUNTER too (setTime only resets the clock). The TS is
        # cached and reused, and getStepNumber() is NOT reset by setTime, so
        # without this it accumulates across calls — making setMaxSteps below a
        # *cumulative* cap that, after ~tsStep total substeps (a few hundred
        # model steps at the default tsStep=2000), is already exceeded on entry,
        # so TSSolve returns immediately and the marine deposit is left
        # un-diffused (silent). It also made the verbose substep/iteration
        # counts grow without bound (they were cumulative, not per-call).
        ts.setStepNumber(0)
        # Initial dt close to the controller's typical equilibrium for
        # marine diffusion; minor over-large warmup is corrected by the
        # adaptive controller within a step or two.
        ts.setTimeStep(self.dt / 100.0)
        ts.setMaxTime(self.dt)
        ts.setMaxSteps(self.tsStep)

        tstart = ts.getTime()
        self._evalSolutionMardDiff(tstart, x)

        ts.solve(x)
        if MPIrank == 0 and self.verbose:
            print(
                "Nonlinear diffusion solution (%s, %0.02f seconds)"
                % (label, process_time() - t0),
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

        # Extract resulting deposit thickness
        x.copy(result=self.h)
        self.dh.waxpy(-1.0, self.hGlobal, self.h)
        self.dm.globalToLocal(self.dh, self.tmpL)
        ndepo = self.tmpL.getArray().copy()
        ndepo[ndepo < 0.0] = 0.0

        return ndepo

    def _diffuseImplicitPicard(self, dh, mask, Cd_val, label="diffusion"):
        r"""
        Lagged-diffusivity (Picard) alternative to ``_diffuseImplicit`` for the
        non-linear thickness diffusion :math:`\partial_t h = \nabla\cdot(C_d(h)
        \nabla h)` over one model step.

        Instead of the adaptive non-linear TS, this takes ``self.picardSub``
        backward-Euler sub-steps of size ``dt/picardSub``; within each sub-step
        it freezes the diffusivity :math:`C_d(h^k)` and solves the resulting
        **linear** system :math:`(I + \Delta t_{sub} L(C_d^k)) h^{k+1} = h^{start}`
        with ``self.picardIts`` Picard updates. The operator ``L`` is built from
        ``jacobiancoeff`` with a zero derivative term (``Kp=0``), so it is the
        exact linear :math:`\nabla\cdot(C_d\nabla)` operator consistent with the
        TS residual ``fctcoeff`` — including the "no flux across a zero-
        diffusivity face" marine-mask gating.

        Each solve is linear (a cached richardson+bjacobi KSP, no SNES) and,
        because :math:`C_d` is frozen, smooth — so it avoids the error-estimate
        rejections the adaptive TS hits at the :math:`C_d` kink (``dh<0.1``).
        ``picardSub`` is the accuracy/speed knob. Same signature / return as
        ``_diffuseImplicit``.

        :arg dh: per-node initial deposit thickness (m)
        :arg mask: cells where diffusion is active
        :arg Cd_val: scalar non-linear diffusion coefficient (m^2/yr)
        :arg label: log label

        :return: ndepo - smoothed deposit thickness (m), zeroed where negative
        """

        t0 = process_time()
        nsub = int(getattr(self, "picardSub", 10))
        npic = int(getattr(self, "picardIts", 2))
        dt_sub = self.dt / nsub

        Cd0 = np.zeros(self.lpoints)
        Cd0[mask] = Cd_val
        minDiff_vec = np.zeros(self.lpoints)
        minDiff_vec[mask] = self.minDiff
        zb = self.hLocal.getArray().copy()         # pre-deposition bed
        zeroKp = np.zeros(self.lpoints)

        # Deposited surface h = bed + dh (global self.h, local self.hl).
        self.hl.setArray(dh)
        self.hl.axpy(1.0, self.hLocal)
        self.dm.localToGlobal(self.hl, self.h)

        if self._ksp_picard is None:
            self._ksp_picard = self._makeDiffusionKSP("marpicard_")
        ksp = self._ksp_picard

        nsolve = 0
        for _ in range(nsub):
            self.h.copy(result=self.tmp1)          # h^start (BE right-hand side)
            for _ in range(npic):
                self.dm.globalToLocal(self.h, self.hl)
                hloc = self.hl.getArray()
                dhdep = hloc - zb
                dhdep[dhdep < 0.1] = 0.0
                Cd = minDiff_vec + np.multiply(
                    Cd0, 1.0 - np.exp(-self.dexp * dhdep)
                )
                # Linear operator L = div(Cd grad .) consistent with fctcoeff
                # (Kp=0 -> no derivative term); 13-col (diag, 12 neighbours).
                nlC = jacobiancoeff(hloc, Cd, zeroKp)
                coeffs = dt_sub * nlC
                coeffs[:, 0] += 1.0                # M = I + dt_sub * L
                M = self._assembleDiffMatCSR(coeffs)
                ksp.setOperators(M, M)
                ksp.solve(self.tmp1, self.h)
                M.destroy()
                nsolve += 1

        # Deposit thickness relative to the pre-deposition elevation.
        self.dh.waxpy(-1.0, self.hGlobal, self.h)
        self.dm.globalToLocal(self.dh, self.tmpL)
        ndepo = self.tmpL.getArray().copy()
        ndepo[ndepo < 0.0] = 0.0

        if MPIrank == 0 and self.verbose:
            print(
                "Nonlinear diffusion (Picard %s, %0.02f seconds): %d linear "
                "solves (%d sub-steps x %d Picard)"
                % (label, process_time() - t0, nsolve, nsub, npic),
                flush=True,
            )
        safe_garbage_cleanup()

        return ndepo
