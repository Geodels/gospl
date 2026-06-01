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

        # KSP solver parameters
        self.rtol = 1.0e-10

        # Cached KSP/PC objects: created lazily on first solve and reused
        # across timesteps so we avoid the create/destroy churn (~5-10ms per
        # solve). PETSc auto-detects matrix changes via setOperators and
        # rebuilds the preconditioner factor when necessary.
        self._ksp_main = None
        self._ksp_fallback = None

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

    def _solve_KSP2(self, matrix, vector1, vector2):
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
            ksp.setType("fgmres")
            ksp.getPC().setType("asm")
            ksp.setTolerances(rtol=1.0e-6, divtol=1.e20)
            ksp.setInitialGuessNonzero(True)
            self._ksp_fallback = ksp

        ksp = self._ksp_fallback
        ksp.setOperators(matrix, matrix)
        ksp.solve(vector1, vector2)
        r = ksp.getConvergedReason()
        if r < 0:
            KSPReasons = self._make_reasons(petsc4py.PETSc.KSP.ConvergedReason())
            if MPIrank == 0:
                print(
                    "LinearSolver failed to converge after iterations",
                    ksp.getIterationNumber(),
                    flush=True,
                )
                print("with reason: ", KSPReasons[r], flush=True)
            # If the fgmres+asm fallback also diverges, return a zero discharge
            # for this step. Operators should monitor the warning above.
            vector2.set(0.0)
        petsc4py.PETSc.garbage_cleanup()

        return vector2

    def _solve_KSP(self, guess, matrix, vector1, vector2):
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
            ksp.setType("richardson")
            ksp.getPC().setType("bjacobi")
            ksp.setTolerances(rtol=self.rtol)
            self._ksp_main = ksp

        ksp = self._ksp_main
        if guess:
            ksp.setInitialGuessNonzero(True)
        ksp.setOperators(matrix, matrix)
        ksp.solve(vector1, vector2)
        r = ksp.getConvergedReason()
        if r < 0:
            vector2 = self._solve_KSP2(matrix, vector1, vector2)
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

        # Define multiple flow directions for unfilled elevation
        self.donRcvs, self.distRcv, self.wghtVal = mfdreceivers(
            self.flowDir, self.flowExp, h, self.sealevel
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

        # Set borders nodes
        if self.flatModel:
            self.rcvID[self.idBorders, :] = np.tile(self.idBorders, (self.flowDir, 1)).T
            self.distRcv[self.idBorders, :] = 0.0
            self.wghtVal[self.idBorders, :] = 0.0

        # Get local nodes with no receivers as boolean array
        sum_weight = np.sum(self.wghtVal, axis=1)
        lsink = sum_weight == 0.0

        # We don't consider open sea nodes and borders as sinks
        lsink[self.idBorders] = False
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

    def _distributeDownstream(self, pitVol, FA, hl, step, ice=False):
        """
        In cases where rivers flow in depressions, they might fill the sink completely and overspill or remain within the depression, forming a lake. This function computes the excess of water (if any) able to flow dowstream.

        .. important::

            The excess water is then added to the downstream flow accumulation (`FA`) and used to estimate rivers' erosion.

        :arg pitVol: volume of depressions
        :arg FA: excess flow accumulation array
        :arg hl: current elevation array
        :arg step: downstream distribution step
        :arg ice: boolean indicating where the ice flow is considered or not.

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
                self._solve_KSP(True, self.fMat, self.tmp, self.tmp1)
                self.dm.globalToLocal(self.tmp1, self.tmpL)
                nFA = self.tmpL.getArray().copy() * self.dt
                FA = nFA.copy()
                FA[hl < self.waterFilled] = 0.0
                self.tmpL.setArray(FA / self.dt)
                if ice:
                    self.iceFAL.axpy(1.0, self.tmpL)
                else:
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
        self.fillElevation(sed=False)
        pitVol = self.pitParams[:, 0].copy()

        # Build flow direction and downstream matrix
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
            elaH = self.elaH(self.tNow)
            iceH = self.iceH(self.tNow)
            tmp = (hl - elaH) / (iceH - elaH)
            tmp[tmp > 1.] = 1.0
            tmp[tmp < 0.] = 0.0
            rainA = np.multiply(rainA, 1. - tmp)
            # Re-inject glacial meltwater captured during iceAccumulation:
            # sub-ELA cells with ice present release the local ablation
            # rate as liquid water into the river source. Without this,
            # melt computed in the ice solve is discarded by the negative
            # clamp on iceFAL and downstream basins under-predict
            # discharge.
            rainA = rainA + self.iceMeltL.getArray()
            rainA[self.seaID] = 0.

        #  Solve flow/ice accumulation
        self.bL.setArray(rainA)
        self.dm.localToGlobal(self.bL, self.bG)
        self._solve_KSP(True, self.fMat, self.bG, self.FAG)
        self.dm.globalToLocal(self.FAG, self.FAL)

        # Volume of water flowing downstream
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

        # Get water level
        self.waterFilled -= hl

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Flow Accumulation (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return
