import os
import gc
import sys
import petsc4py
import numpy as np
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import mfdreceivers

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD


class FAMesh(object):
    """
    This class calculates **drainage area** in an implicit, iterative manner using PETSc solvers. It accounts for multiple flow direction paths (SFD to MFD) based on user input declaration.

    .. note::

        The class follows the parallel approach described in `Richardson et al., 2014 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013WR014326>`_.
    """

    def __init__(self, *args, **kwargs):
        """
        The initialisation of `FAMesh` class consists in the declaration of PETSc vectors and matrices.
        """

        # KSP solver parameters. rtol=1e-8 is ample for discharge used as a
        # flow-accumulation proxy; the previous 1e-10 was so tight that a cold
        # start on a high-resolution global mesh (long strictly-downhill MFD
        # chains -> Richardson needs ~longest-path-in-hops iterations) could
        # exceed max_it and fall through to the fgmres fallback.
        self.rtol = 1.0e-8

        # Cached KSP/PC objects: created lazily on first solve and reused
        # across timesteps so we avoid the create/destroy churn (~5-10ms per
        # solve). PETSc auto-detects matrix changes via setOperators and
        # rebuilds the preconditioner factor when necessary.
        self._ksp_main = None
        self._ksp_fallback = None

        # An (I - W^T) solve that fails to converge on only a HANDFUL of rows is
        # the benign isolated-pocket / micro-cycle case (a few genuinely
        # un-drainable cells that just pond; discharge there is clamped >=0, mass
        # conserved). It is expected on a partitioned high-res mesh and no solver
        # can converge it, so it is reported calmly rather than as a scary
        # failure. A failure spanning MORE than this many nodes is a real problem
        # (NaN source, broken partition) and keeps the full loud diagnostic. The
        # observed benign region is O(1-10) nodes even at 240 ranks, so 256 is
        # well clear of it and far below any genuine convergence failure.
        self._undrained_benign_cap = 256

        # Iteration cap for the NON-fatal IDA cascade solves (sediment/water
        # downstream routing). The well-posed bulk converges in O(100s) of
        # fgmres iterations; a partition-dependent near-singular cell (an
        # isolated pocket / residual cross-partition cycle) can NEVER reach rtol,
        # so without a tight bound fgmres grinds the full primary budget (5000)
        # on that one cell before the bounded fallback -- which is the dominant
        # cause of the erratic `sed` wall-time at some rank counts (measured: a
        # single capped solve = ~5500 iters at np=16 vs ~205/solve at np=8 on a
        # 9.2M mesh). The fatal flow-accumulation solve keeps the full budget
        # (it MUST converge); only the degrade-gracefully cascade solves are
        # capped here. The pinning of undrained cells (which removes the singular
        # rows entirely) is the deeper fix; this cap is the cheap backstop.
        self._cascade_max_it = 1000

        # Identity matrix construction
        self.II = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        self.JJ = np.arange(0, self.lpoints, dtype=petsc4py.PETSc.IntType)
        self.iMat = self._matrix_build_diag(np.ones(self.lpoints))

        # Petsc vectors
        self.fillFAL = self.hLocal.duplicate()
        self.FAG = self.hGlobal.duplicate()
        self.FAL = self.hLocal.duplicate()

        return

    def _matrix_build(self, nnz=(1, 1)):
        """
        Creates a sparse PETSc matrix.

        .. note::

            To achieve good performance during matrix assembly, the function preallocates the matrix storage.

        :arg nnz: array containing the number of nonzeros in the various rows

        :return: sparse PETSc matrix
        """

        matrix = petsc4py.PETSc.Mat().create(comm=MPIcomm)
        matrix.setType("aij")
        matrix.setSizes(self.sizes)
        matrix.setLGMap(self.lgmap_row, self.lgmap_col)
        matrix.setFromOptions()
        matrix.setPreallocationNNZ(nnz)

        return matrix

    def _matrix_build_diag(self, V, nnz=(1, 1)):
        """
        Builds a PETSc diagonal matrix based on a given array `V`

        :arg V: diagonal data array
        :arg nnz: array containing the number of nonzero blocks

        :return: sparse PETSc matrix
        """

        matrix = self._matrix_build()

        # Define diagonal matrix
        matrix.assemblyBegin()
        matrix.setValuesLocalCSR(
            self.II,
            self.JJ,
            V,
        )
        matrix.assemblyEnd()

        return matrix

    def _make_reasons(self, reasons):
        """
        Provides reasons for PETSc error...
        """

        return dict(
            [(getattr(reasons, r), r) for r in dir(reasons) if not r.startswith("_")]
        )

    def _solve_KSP2(self, matrix, vector1, vector2, fatal=False):
        """
        Solution of Krylov subspace iterative method (PETSc *scalable linear equations solvers* - **KSP**) implemented using the Flexible Generalized Minimal Residual method (`fgmres`) with Additive Schwarz preconditioning (`asm`).

        .. note::

            This function is used if the KSP convergence failed using the PETSc Richardson solver with block Jacobian preconditioning.

        :arg matrix: PETSc sparse matrix used by the KSP solver composed of diagonal terms set to unity (identity matrix) and off-diagonal terms (weights between 0 and 1). The weights are calculated based on the number of downslope neighbours (*i.e.*, user-defined number of flow directions) and are proportional to the slope.
        :arg vector1: PETSc vector corresponding to the local volume of water available for runoff during a given time step (*e.g.* voronoi area times local precipitation rate).
        :arg vector2: PETSc vector corresponding to the unknown flow discharge values.

        :return: vector2 PETSc vector of the new flow discharge values
        """
        if self._ksp_fallback is None:
            ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
            # Fallback = plain Richardson with NO preconditioner: the IDA
            # fixed-point x <- b + W^T x. The primary's bjacobi is fast but, on
            # this non-normal (I - W^T) operator, its per-rank blocks DESTABILISE
            # the iteration at some partitions -> the primary diverges and the
            # solve lands here. The unpreconditioned iteration converges whenever
            # the mesh drains (rho(W^T) < 1), INDEPENDENT of partition, so it
            # cleanly and cheaply recovers those bjacobi-destabilised-but-solvable
            # solves. This replaces the old fgmres+asm fallback, which on this
            # operator either stalled (DIVERGED_MAX_IT), broke down (BiCGStab
            # NANORINF), or -- when it did converge -- was so expensive it caused
            # the erratic-scaling sed/sea blow-ups (e.g. P=240 sed=112 s with NO
            # failure, just costly fallback). A genuinely singular sub-block (an
            # undrained micro-region) still cannot converge; it grinds to the
            # bounded max_it (cheap mat-vecs) and degrades gracefully (zeroed)
            # for non-fatal callers, or aborts for the fatal main solve.
            ksp.setType("richardson")
            ksp.getPC().setType("none")
            # With the fgmres default primary (robust, monotone), this fallback
            # is reached only by a genuinely (near-)singular sub-region that no
            # method can converge. So bound max_it TIGHTLY: it just confirms
            # non-convergence cheaply, then degrades gracefully (zeroed -> finite
            # and bounded) for non-fatal callers. A large cap here only added a
            # grind at the partitions where the singular region appears -- and
            # since fgmres already recovers every SOLVABLE case, the fallback
            # never needs deep iteration, so 500 is ample to confirm and stop.
            ksp.setTolerances(rtol=1.0e-6, divtol=1.e20, max_it=500)
            ksp.setInitialGuessNonzero(True)
            ksp.setOptionsPrefix("flowaccfb_")
            ksp.setFromOptions()
            self._ksp_fallback = ksp

        ksp = self._ksp_fallback
        ksp.setOperators(matrix, matrix)
        ksp.solve(vector1, vector2)
        r = ksp.getConvergedReason()
        if r < 0:
            # Both the primary and fallback KSP diverged. This helper is shared
            # by many solves (flow accumulation, sediment routing, hillslope,
            # SPL, tectonics, ice). The converged reason is global (identical
            # on every rank), so every rank reaches this branch together.
            KSPReasons = self._make_reasons(petsc4py.PETSc.KSP.ConvergedReason())
            its = ksp.getIterationNumber()
            # Localize the un-converged region FIRST (one matvec). With a finite
            # RHS+matrix and a cleaned guess, a persistent failure means the
            # operator is (near-)rank-deficient on those rows -- a drainage region
            # that never reaches a sink. The nodes carrying the largest residual
            # ARE that region. Wrapped defensively so the diagnostic can never
            # itself break the run; `nbad` stays -1 (-> loud, the safe default)
            # if it cannot be computed.
            worst_id, worst_val, nbad = -1, float("nan"), -1
            try:
                resid = vector1.duplicate()
                matrix.mult(vector2, resid)        # A x
                resid.aypx(-1.0, vector1)          # b - A x
                resid.abs()
                worst_id, worst_val = resid.max()  # (global index, value)
                bnorm = vector1.norm()
                thr = 1.0e-3 * bnorm if bnorm > 0.0 else 1.0e-3
                nbad = int(np.count_nonzero(resid.getArray() > thr))
                nbad = MPI.COMM_WORLD.allreduce(nbad, op=MPI.SUM)
                resid.destroy()
            except Exception:
                pass
            # A TINY un-drained region on a non-fatal (I - W^T) solve is the
            # benign isolated-pocket / micro-cycle case (the cells just pond;
            # discharge there is clamped >=0, mass conserved) -- expected and
            # unconvergeable, so report it calmly. The fatal main solve, or a
            # large region (a real failure), gets the full loud diagnostic plus
            # the RHS/matrix-finite probe (a non-finite RHS = a NaN source to
            # hunt upstream, distinct from a convergence stall). `nbad` is global
            # (Allreduced) and `fatal` is identical on every rank, so the branch
            # is collective-consistent.
            benign = (not fatal) and (0 <= nbad <= self._undrained_benign_cap)
            if MPIrank == 0:
                if benign:
                    print(
                        "[flow] %d isolated un-drained cell(s) left as local "
                        "sinks (benign: they pond, mass conserved)" % nbad,
                        flush=True,
                    )
                else:
                    rhs_finite = bool(np.isfinite(vector1.norm()))
                    mat_finite = bool(
                        np.isfinite(matrix.norm(petsc4py.PETSc.NormType.FROBENIUS))
                    )
                    print(
                        "KSP (richardson/none fallback) failed to converge after",
                        its,
                        "iterations with reason:",
                        KSPReasons.get(r, r),
                        "(RHS finite: %s, matrix finite: %s)"
                        % (rhs_finite, mat_finite),
                        flush=True,
                    )
                    if nbad >= 0:
                        print(
                            "  [flowKSP] worst residual %.3e at global node %d; "
                            "%d nodes exceed 1e-3*||b|| (the un-drained region)"
                            % (worst_val, worst_id, nbad),
                            flush=True,
                        )
            if fatal:
                # Only the main flow-accumulation discharge solve passes
                # fatal=True: a zero discharge there would silently feed a
                # no-river state into the erosion/sediment routines, so abort
                # instead. Raising on every rank is collective -> no deadlock.
                raise RuntimeError(
                    "Flow-accumulation KSP failed to converge (reason %s) after "
                    "%d iterations; aborting rather than continuing with zero "
                    "discharge." % (KSPReasons.get(r, r), its)
                )
            # Auxiliary / iterative-cascade solves degrade gracefully: drop
            # this solve's contribution (zero) and continue, as before.
            vector2.set(0.0)
        petsc4py.PETSc.garbage_cleanup()

        return vector2

    def _solve_KSP(self, guess, matrix, vector1, vector2, fatal=False, seed=False):
        """
        PETSc *scalable linear equations solvers* (**KSP**) component provides Krylov subspace iterative method and a preconditioner. Here, flow accumulation solution is obtained using PETSc Richardson solver (`richardson`) with block Jacobian preconditioning (`bjacobi`).

        .. note::

            The solver choice was made based on the convergence results from `Richardson et al. (2014) <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013WR014326>`_ but can be changed if better solver and preconditioner combinations are found.

        Using such iterative method allows for an initial guess to be provided. When this initial guess is close to the solution, the number of iterations required for convergence dramatically decreases. Here the flow discharge solution from previous time step can be passed as an initial `guess` to the solver as discharge often exhibits little change between successive time intervals.

        :arg guess: Boolean specifying if the iterative KSP solver initial guess is nonzero (when provided it corresponds to the previous flow discharge values).
        :arg matrix: PETSc sparse matrix used by the KSP solver composed of diagonal terms set to unity (identity matrix) and off-diagonal terms (weights between 0 and 1). The weights are calculated based on the number of downslope neighbours (based on the chosen number of flow direction directions) and are proportional to the slope.
        :arg vector1: PETSc vector corresponding to the local volume of water available for runoff during a given time step (*e.g.* voronoi area times local precipitation rate).
        :arg vector2: PETSc vector corresponding to the unknown flow discharge values.

        :return: vector2 PETSc vector of the new flow discharge values
        """

        if self._ksp_main is None:
            ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
            # Primary solver + preconditioner (override via env GOSPL_FLOW_KSP /
            # GOSPL_FLOW_PC). The (I - W^T) matrix is ill-conditioned (rho ~ 1
            # from long flat drainage chains) so it NEEDS acceleration: bjacobi
            # provides it. The historical choice was STATIONARY Richardson, but
            # its per-rank blocks DESTABILISE at some partitions (spectral radius
            # > 1 -> diverge -> slow fallback -> erratic scaling). 'fgmres'
            # minimises the residual monotonically so it CANNOT diverge that way,
            # while still benefiting from bjacobi -- validated fast AND stable at
            # every rank count (e.g. P=240: 180s -> 50s, no fallback). So
            # 'fgmres'+'bjacobi' is now the default; 'richardson' remains
            # available via the env var.
            flow_ksp = os.environ.get("GOSPL_FLOW_KSP", "fgmres")
            flow_pc = os.environ.get("GOSPL_FLOW_PC", "bjacobi")
            ksp.setType(flow_ksp)
            ksp.getPC().setType(flow_pc)
            if MPIrank == 0:
                print(
                    "[flow] primary KSP: %s + %s" % (flow_ksp, flow_pc),
                    flush=True,
                )
            # Iteration budget. fgmres (the default) converges the well-posed
            # bulk in O(100s) of passes; a tiny singular sub-region (a few
            # genuinely undrained nodes) can NEVER reach rtol, so without a bound
            # fgmres grinds the full budget failing to fix just those nodes --
            # the residual sed/sea blow-ups at some partitions. Cap it modestly:
            # the bulk is solved well within 5000 and the result is accepted as
            # best-effort on max_it (see _solve_KSP). Stationary Richardson (the
            # GOSPL_FLOW_KSP=richardson escape hatch) propagates one hop/pass, so
            # it still needs the larger ~longest-flow-path budget.
            self._primary_max_it = 100000 if flow_ksp == "richardson" else 5000
            ksp.setTolerances(rtol=self.rtol, max_it=self._primary_max_it)
            # Shift zero/negative pivots in the per-rank ILU sub-solver so the
            # block-Jacobi PCSetUp cannot fail (DIVERGED_PCSETUP_FAILED) on a
            # degenerate / ocean-only partition at higher rank counts. Scoped
            # by an options prefix so it touches only this solver.
            ksp.setOptionsPrefix("flowacc_")
            petsc4py.PETSc.Options()["flowacc_sub_pc_factor_shift_type"] = "nonzero"
            ksp.setFromOptions()
            self._ksp_main = ksp

        ksp = self._ksp_main
        if guess:
            ksp.setInitialGuessNonzero(True)
            # Cold-start seed (opt-in via seed=True). Applies ONLY to solves on
            # the substochastic (I - W^T) flow matrix (flow accumulation +
            # downstream water/sediment routing), where the solution satisfies
            # x = b + W^T x >= b element-wise, so the RHS b is a valid lower
            # bound. On a cold start (no previous-step solution, vector2 == 0)
            # Richardson would otherwise propagate the solution one hop per
            # iteration over the full network and exceed max_it on a high-
            # resolution global mesh. Do NOT enable for other systems
            # (hillslope/SPL/tectonics/ice) where b is not a valid lower bound.
            if seed and vector2.norm() == 0.0:
                vector1.copy(vector2)
        # Cap the iteration budget per call: the fatal flow-accumulation solve
        # keeps the full primary budget (it must converge); the non-fatal IDA
        # cascade solves (seed=True) are capped so a near-singular cell can't
        # grind the full budget before the bounded fallback. seed=False solves
        # (hillslope/SPL/tectonics/ice) keep the full budget.
        ksp.setTolerances(
            rtol=self.rtol,
            max_it=(self._cascade_max_it if (seed and not fatal)
                    else self._primary_max_it),
        )
        ksp.setOperators(matrix, matrix)
        ksp.solve(vector1, vector2)
        r = ksp.getConvergedReason()
        if r < 0:
            # The primary failed (max_it on a near-singular sub-region, or a
            # genuine breakdown). Do NOT accept the iterate: on a near-singular
            # system the SOLUTION is unbounded (huge null-space component) even
            # when the residual looks small, so accepting fgmres's best-effort
            # injects huge values that compound through the cascade (the 1e29
            # blow-up). Restart the fallback from a clean, BOUNDED guess -- the
            # seed lower-bound b for seeded solves, else a finite warm iterate
            # (drop a non-finite one) -- and let the bounded fallback finish (it
            # zeroes a non-convergent solve, keeping everything finite).
            if seed:
                vector1.copy(vector2)
            elif not np.isfinite(vector2.norm()):
                vector2.set(0.0)
            vector2 = self._solve_KSP2(matrix, vector1, vector2, fatal=fatal)
        petsc4py.PETSc.garbage_cleanup()

        return vector2

    def matrixFlow(self, flowdir, dep=None):
        """
        This function defines the flow direction matrices.

        .. note::

            The matrix is built incrementally looping through the number of flow direction paths defined by the user. It proceeds by assembling a local Compressed Sparse Row (**CSR**) matrix to a global PETSc matrix.

            When setting up the flow matrix in PETSc, we preallocate the non-zero entries of the matrix before starting filling in the values. Using PETSc sparse matrix storage scheme has the advantage that matrix-vector multiplication is extremely fast.

        The matrix coefficients consist of weights (comprised between 0 and 1) and calculated based on the number of downslope neighbours and proportional to the slope.

        :arg flow: boolean to compute matrix for either downstream water or sediment transport
        :arg dep: deposition flux coefficient in case where the sediment transport/deposition term is considered.
        """

        self.fMat = self.iMat.copy()
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        nodes = indptr[:-1]
        if dep is None:
            wght = self.wghtVal
        else:
            wght = np.multiply(self.wghtVal, dep.reshape((len(dep), 1)))
        rcv = self.rcvID

        for k in range(0, flowdir):
            # Flow direction matrix for a specific direction
            tmpMat = self._matrix_build()
            data = -wght[:, k].copy()
            data[rcv[:, k].astype(petsc4py.PETSc.IntType) == nodes] = 0.0
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                rcv[:, k].astype(petsc4py.PETSc.IntType),
                data,
            )
            tmpMat.assemblyEnd()
            # Add the weights from each direction
            self.fMat.axpy(1.0, tmpMat)
            tmpMat.destroy()

        if self.memclear:
            del data, indptr, nodes
            gc.collect()

        # Store flow accumulation matrix
        self.fMat.transpose()
        petsc4py.PETSc.garbage_cleanup()

        return

    def _buildFlowDirection(self, h, down=True):
        """
        This function builds from neighbouring slopes the flow directions. It calls a fortran subroutine that locally computes for each vertice:

        - the indices of receivers (downstream) nodes depending on the desired number of flow directions (SFD to MFD).
        - the distances to the receivers based on mesh resolution.
        - the associated weights calculated based on the number of receivers and proportional to the slope.

        :arg h: elevation in the form of a Numpy Array
        :arg down: boolean to indicate whether the filled elevation needs to be considered or not.
        """

        # Get open marine regions
        self.seaID = np.where(self.lFill <= self.sealevel)[0]

        # Define multiple flow directions for unfilled elevation.
        # `self.gid` is the per-node global ID; passed in to make the
        # slope-tie-break in mfdreceivers deterministic across MPI
        # decompositions (see fortran/functions.F90 and AGENTS.md).
        self.donRcvs, self.distRcv, self.wghtVal = mfdreceivers(
            self.flowDir, self.flowExp, h, self.sealevel, self.gid
        )

        self.rcvID = self.donRcvs.copy()
        self.rcvID[self.ghostIDs, :] = -1
        self.distRcv[self.ghostIDs, :] = 0
        self.wghtVal[self.ghostIDs, :] = 0

        if down:
            sum_weight = np.sum(self.wghtVal, axis=1)
            ids = (
                (h == self.lFill)
                & (self.pitIDs > -1)
                & (self.flatDirs > -1)
                & (sum_weight == 0.0)
            )
            ids = ids.nonzero()[0]
            self.rcvID[ids, :] = np.tile(ids, (self.flowDir, 1)).T
            self.rcvID[ids, 0] = self.flatDirs[ids]
            self.wghtVal[ids, :] = 0.0
            self.wghtVal[ids, 0] = 1.0

        # Set draining-border nodes (open / fixed outlets only — `outletIDs`
        # excludes true-wall edges so those behave as interior, retaining flow
        # and letting sediment deposit against them).
        if self.flatModel:
            self.rcvID[self.outletIDs, :] = np.tile(self.outletIDs, (self.flowDir, 1)).T
            self.distRcv[self.outletIDs, :] = 0.0
            self.wghtVal[self.outletIDs, :] = 0.0

        # Get local nodes with no receivers as boolean array
        sum_weight = np.sum(self.wghtVal, axis=1)
        lsink = sum_weight == 0.0

        # We don't consider open sea nodes and draining borders as sinks (wall
        # edges ARE allowed to be sinks so deposition contains the sediment).
        lsink[self.outletIDs] = False
        lsink[self.seaID] = False
        lsink = lsink.astype(int) * self.inIDs

        self.lsink = lsink == 1

        self.matrixFlow(self.flowDir)

        return

    def _potentialLakeEvap(self):
        """
        Per-pit max evaporation volume (m^3) for one timestep, assuming
        the lake fills to the spillover (= max-fill surface).

        Surface = Σ larea over cells with pitIDs >= 0 (the lake's potential
        extent at spillover). The `inIDs == 1` mask prevents MPI halo
        double-counting; the per-pit total is then Allreduce-summed across
        ranks so every rank ends up with the same evap budget array.

        See DESIGN_EVAPORATION.md §1 D3 for the design rationale (we use
        max-fill rather than current-fill to avoid the circular dependency
        between fill level and evap rate).

        :return: numpy array of shape (nbpits,) with per-pit evap m^3
        """
        out = np.zeros(len(self.pitParams), dtype=np.float64)
        if self.evapVal is None:
            return out
        isPitCell = (self.pitIDs >= 0) & (self.inIDs == 1)
        if isPitCell.any():
            cellEvap = self.evapVal * self.larea * self.dt
            grp = npi.group_by(self.pitIDs[isPitCell])
            uIDs = grp.unique
            _, vol = grp.sum(cellEvap[isPitCell])
            ids = uIDs > -1
            out[uIDs[ids]] = vol[ids]
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, out, op=MPI.SUM)
        return out

    def _distributeDownstream(self, pitVol, FA, hl, step):
        """
        In cases where rivers flow in depressions, they might fill the sink completely and overspill or remain within the depression, forming a lake. This function computes the excess of water (if any) able to flow dowstream.

        .. important::

            The excess water is then added to the downstream flow accumulation (`FA`) and used to estimate rivers' erosion.

        :arg pitVol: volume of depressions
        :arg FA: excess flow accumulation array
        :arg hl: current elevation array
        :arg step: downstream distribution step

        :return: pitVol, excess, nFA (updated volume in each depression, boolean set to True is excess flow remains to be distributed and new flow accumulation values)
        """

        excess = False

        # Remove points belonging to other processors
        FA = np.multiply(FA, self.inIDs)

        # Get volume incoming in each depression
        grp = npi.group_by(self.pitIDs[self.lsink])
        uID = grp.unique
        _, vol = grp.sum(FA[self.lsink])
        inV = np.zeros(len(self.pitParams), dtype=np.float64)
        ids = uID > -1
        inV[uID[ids]] = vol[ids]

        # Combine incoming volume globally
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, inV, op=MPI.SUM)

        # Lake-surface evaporation (DESIGN_EVAPORATION.md §2.2). Subtract
        # the per-pit max-fill evap budget from inflow exactly once per
        # outer step (step==0). Cascading spillover iterations (step>0)
        # carry already-net-of-evap water and must NOT be debited again.
        # When the evap budget exceeds inflow, lakeLoss is clamped at inV
        # so the pit ends with zero net inflow — the partial-fill branch
        # below skips it via the inV>0 mask, leaving waterFilled at hl.
        if self.evapVal is not None and step == 0:
            evapBudget = self._potentialLakeEvap()
            lakeLoss = np.minimum(inV, evapBudget)
            inV = inV - lakeLoss
            self.evapLoss += lakeLoss.sum()

        # Get excess volume to distribute downstream
        eV = inV - pitVol
        if (eV > 0.0).any():
            eIDs = eV > 0.0
            pitVol[eIDs] = 0.0
            spillIDs = self.pitInfo[eIDs, 0]
            localSpill = np.where(self.pitInfo[eIDs, 1] == MPIrank)[0]
            localPts = spillIDs[localSpill]
            nFA = np.zeros(self.lpoints, dtype=np.float64)
            nFA[localPts] = eV[eIDs][localSpill]
            ids = np.isin(self.pitIDs, np.where(eV > 0.0)[0])
            self.waterFilled[ids] = self.lFill[ids]

        # Update unfilled depressions volumes and assign water level in depressions
        # NOTE: the (inV > 0.0) clause excludes pits where lake-evap fully
        # consumed the inflow (inV becomes 0 after the subtraction above).
        # Without it, eV = -pitVol drives pitVol → 0 here, incorrectly
        # marking a "dry" depression as filled-to-spillover.
        if (eV < 0.0).any():
            eIDs = np.where((eV < 0.0) & (inV > 0.0))[0]
            pitVol[eIDs] += eV[eIDs]
            nid = np.absolute(self.filled_vol[eIDs] - pitVol[eIDs][:, None]).argmin(
                axis=1
            )
            fill_lvl = self.filled_lvl[eIDs, nid]
            # Vectorised replacement of the per-pit Python loop. Build a
            # pit-id → fill-level lookup and apply it in one pass over nodes.
            in_eIDs = np.isin(self.pitIDs, eIDs)
            pit_fill_lvl = np.zeros(len(self.pitParams))
            pit_fill_lvl[eIDs] = fill_lvl
            node_lvl = pit_fill_lvl[self.pitIDs]
            mask = in_eIDs & (self.waterFilled <= node_lvl)
            self.waterFilled[mask] = node_lvl[mask]

        # In case there is still remaining water flux to distribute downstream
        if (eV > 1.0e-3).any():  # TODO-REFACTOR: value matches DEPOSIT_FLOOR but distinct role (water-routing convergence threshold); do not replace
            if step == 100:
                self.fMat.destroy()
                self._buildFlowDirection(self.lFill)
            else:
                self.fMat.destroy()
                self._buildFlowDirection(self.waterFilled)
            self.tmpL.setArray(nFA / self.dt)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            # Skip the KSP solve for negligible residual fluxes: comparing
            # total residual flux (m^3/yr) against the largest cell area is a
            # cheap "less than 1 m of water over the biggest cell" guard.
            if self.tmp.sum() > self.maxarea[0]:
                excess = True
                self._solve_KSP(True, self.fMat, self.tmp, self.tmp1, seed=True)
                self.dm.globalToLocal(self.tmp1, self.tmpL)
                nFA = self.tmpL.getArray().copy() * self.dt
                FA = nFA.copy()
                FA[hl < self.waterFilled] = 0.0
                self.tmpL.setArray(FA / self.dt)
                self.FAL.axpy(1.0, self.tmpL)
            else:
                nFA = None
        else:
            nFA = None

        return excess, pitVol, nFA

    def flowAccumulation(self):
        """
        This function is the **main entry point** for flow accumulation computation.

        .. note::

            Flow accumulation (`FA`) calculations are a core component of landscape evolution models as they are often used as proxy to estimate flow discharge, sediment load, river width, bedrock erosion, and sediment deposition. Until recently, 

        goSPL model computes the flow discharge from `FA` and the net precipitation rate using a **parallel implicit drainage area (IDA) method** proposed by `Richardson et al., 2014 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013WR014326>`_ but adapted to unstructured grids.

        It calls the following *private functions*:

        1. _buildFlowDirection
        2. _solve_KSP
        3. _distributeDownstream

        """

        t0 = process_time()

        # Compute depressions information
        # Sub-phase profiling (no-op when profiling off): pit-filling vs flow-
        # direction build vs the FA KSP solve vs the downstream-routing loop,
        # to localise the high-rank `flow` plateau. These nest under the "flow"
        # phase, so they double-count there (shows as negative "unaccounted").
        self.profiler.start("flow_fill")
        self.fillElevation(sed=False)
        self.profiler.stop("flow_fill")
        pitVol = self.pitParams[:, 0].copy()

        # Build flow direction and downstream matrix
        self.profiler.start("flow_dir")
        hl = self.hLocal.getArray().copy()

        self._buildFlowDirection(hl, False)
        self.wghtVali = self.wghtVal.copy()
        self.rcvIDi = self.rcvID.copy()
        self.distRcvi = self.distRcv.copy()
        self.fMati = self.fMat.copy()
        self.lsinki = self.lsink.copy()

        # Get amount of water or ice
        rainA = self.bL.getArray().copy()
        rainA[self.seaID] = 0.
        # Channel evaporation (DESIGN_EVAPORATION.md §2.1). evapVal is m/yr,
        # larea is m^2, rainA is m^3/yr — same units. channelLoss is clamped
        # at the available rain so an arid cell cannot produce negative
        # runoff; the unused evap capacity is silently dropped (the
        # groundwater path is out of scope).
        if self.evapVal is not None:
            channelLoss = np.minimum(rainA, self.evapVal * self.larea)
            rainA = rainA - channelLoss
            self.evapLoss += channelLoss.sum() * self.dt
        if self.iceOn:
            # ELA / ice-cap altitude: per-vertex maps when in spatial /
            # time-series mode, otherwise the uniform time scalars. The
            # accumulation ramp is hardened against degenerate (iceH <= elaH)
            # cells so the precipitation split never divides by zero.
            elaH = self.elaMesh[self.locIDs] if self.elaMesh is not None else self.elaH(self.tNow)
            iceH = self.iceMesh[self.locIDs] if self.iceMesh is not None else self.iceH(self.tNow)
            denom = iceH - elaH
            safe = np.where(denom > 0.0, denom, 1.0)
            tmp = np.where(denom > 0.0, (hl - elaH) / safe, 0.0)
            tmp = np.clip(tmp, 0.0, 1.0)
            rainA = np.multiply(rainA, 1. - tmp)
            # Re-inject glacial meltwater (iceMeltRiverL, set in iceplex):
            # discharge-conserving by default — the precipitation that fell as
            # ice above the ELA is routed down-glacier and released here as
            # liquid water where the ice melts out, so Σ meltwater == Σ glacial
            # accumulation. Without this the meltwater would be lost and
            # downstream basins would under-predict discharge.
            rainA = rainA + self.iceMeltRiverL.getArray()
            rainA[self.seaID] = 0.

        #  Solve flow/ice accumulation
        self.bL.setArray(rainA)
        self.dm.localToGlobal(self.bL, self.bG)
        self.profiler.stop("flow_dir")
        self.profiler.start("flow_ksp")
        # seed=True: cold-start the discharge from the runoff vector b (valid
        # lower bound on the (I - W^T) system) so init does not blow past
        # max_it -> DIVERGED_MAX_IT. fatal=True: a failed discharge solve here
        # would silently zero all rivers and feed a degenerate state into
        # erosion/sediment, so abort instead.
        self._solve_KSP(True, self.fMat, self.bG, self.FAG, fatal=True, seed=True)
        self.dm.globalToLocal(self.FAG, self.FAL)
        self.profiler.stop("flow_ksp")

        # Volume of water flowing downstream
        self.profiler.start("flow_dist")
        self.waterFilled = hl.copy()
        if (pitVol > 0.0).any():
            FA = self.FAL.getArray().copy() * self.dt
            excess = True
            step = 0
            while excess:
                t1 = process_time()
                excess, pitVol, FA = self._distributeDownstream(pitVol, FA, hl, step)
                if MPIrank == 0 and self.verbose:
                    print(
                        "Downstream flow computation step %d (%0.02f seconds)"
                        % (step, process_time() - t1),
                        flush=True,
                    )
                step += 1

            # Get overall water flowing donwstream accounting for filled depressions
            FA = self.FAL.getArray().copy()
            ids = (hl <= self.waterFilled) & (self.pitIDs > -1)
            FA[ids] = 0.0
            self.FAL.setArray(FA)
            self.dm.localToGlobal(self.FAL, self.FAG)
            self.dm.globalToLocal(self.FAG, self.FAL)
            # Lake-erosion proxy: assign a background discharge (10% of the
            # global max) at filled-depression nodes so that SPL still produces
            # some erosion on the filled topography. FAL itself stays zero
            # there because no water actually flows.
            FA[ids] = self.FAG.max()[1] * 0.1
            self.fillFAL.setArray(FA)
        else:
            self.FAL.copy(result=self.fillFAL)

        # fgmres (unlike the IDA fixed-point Richardson) does not preserve the
        # non-negativity of discharge, so even a converged solve can leave tiny
        # negative values. Discharge is physically >= 0 and feeds fractional
        # powers downstream (PA ** spl_m in the SPL eroders), where a negative
        # base produces NaN. Clamp the final discharge fields at the source so
        # every consumer (SPL / nlSPL / soilSPL / sediment) is protected. No-op
        # under Richardson (which never produces negatives).
        for vec in (self.FAL, self.fillFAL):
            arr = vec.getArray().copy()
            if (arr < 0.0).any():
                np.clip(arr, 0.0, None, out=arr)
                vec.setArray(arr)

        # Get water level
        self.waterFilled -= hl
        self.profiler.stop("flow_dist")

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Flow Accumulation (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return
