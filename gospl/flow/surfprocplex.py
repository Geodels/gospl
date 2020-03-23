import os
import gc
import sys
import petsc4py
import numpy as np

from time import clock
from mpi4py import MPI
from scipy import sparse
from petsc4py import PETSc

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import fillPIT
    from gospl._fortran import MFDreceivers
    from gospl._fortran import setHillslopeCoeff
    from gospl._fortran import setDiffusionCoeff

petsc4py.init(sys.argv)
MPIrank = PETSc.COMM_WORLD.Get_rank()
MPIcomm = PETSc.COMM_WORLD


class SPMesh(object):
    """
    Performing surface evolution induced by considered processes:

     - hillslope processes
     - rivers using stream power law
    """

    def __init__(self, *args, **kwargs):

        # KSP solver parameters
        self.rtol = 1.0e-8
        self.hbot = -500.0

        # Identity matrix construction
        self.II = np.arange(0, self.npoints + 1, dtype=PETSc.IntType)
        self.JJ = np.arange(0, self.npoints, dtype=PETSc.IntType)
        self.iMat = self._matrix_build_diag(np.ones(self.npoints))

        # Petsc vectors
        self.FillG = self.hGlobal.duplicate()
        self.FillL = self.hLocal.duplicate()
        self.FAG = self.hGlobal.duplicate()
        self.FAL = self.hLocal.duplicate()
        self.hOld = self.hGlobal.duplicate()
        self.hOldLocal = self.hLocal.duplicate()
        self.stepED = self.hGlobal.duplicate()
        self.tmp = self.hGlobal.duplicate()
        self.tmpL = self.hLocal.duplicate()
        self.vGlob = self.hGlobal.duplicate()
        self.Eb = self.hGlobal.duplicate()
        self.EbLocal = self.hLocal.duplicate()
        self.vSed = self.hGlobal.duplicate()
        self.vSedLocal = self.hLocal.duplicate()

        # Diffusion matrix construction
        diffCoeffs, self.maxnb = setHillslopeCoeff(self.npoints, self.Cd * self.dt)
        self.Diff = self._matrix_build_diag(diffCoeffs[:, 0])

        for k in range(0, self.maxnb):
            tmpMat = self._matrix_build()
            indptr = np.arange(0, self.npoints + 1, dtype=PETSc.IntType)
            indices = self.FVmesh_ngbID[:, k].copy()
            data = np.zeros(self.npoints)
            ids = np.nonzero(indices < 0)
            indices[ids] = ids
            data = diffCoeffs[:, k + 1]
            ids = np.nonzero(data == 0.0)
            indices[ids] = ids
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                indices.astype(PETSc.IntType),
                data,
                PETSc.InsertMode.INSERT_VALUES,
            )
            tmpMat.assemblyEnd()
            self.Diff += tmpMat
            tmpMat.destroy()
        del ids, indices, indptr, diffCoeffs
        gc.collect()

        return

    def _matrix_build(self, nnz=(1, 1)):
        """
        Define PETSC Matrix.

        :arg nnz: array containing the number of nonzero blocks
        """

        matrix = PETSc.Mat().create(comm=MPIcomm)
        matrix.setType("aij")
        matrix.setSizes(self.sizes)
        matrix.setLGMap(self.lgmap_row, self.lgmap_col)
        matrix.setFromOptions()
        matrix.setPreallocationNNZ(nnz)

        return matrix

    def _matrix_build_diag(self, V, nnz=(1, 1)):
        """
        Define PETSC Diagonal Matrix.

        :arg V: diagonal data array
        :arg nnz: array containing the number of nonzero blocks
        """

        matrix = self._matrix_build()

        # Define diagonal matrix
        matrix.assemblyBegin()
        matrix.setValuesLocalCSR(self.II, self.JJ, V, PETSc.InsertMode.INSERT_VALUES)
        matrix.assemblyEnd()

        return matrix

    def _make_reasons(self, reasons):
        """
        Provide reasons for PETSC error if possible...
        """

        return dict(
            [(getattr(reasons, r), r) for r in dir(reasons) if not r.startswith("_")]
        )

    def _solve_KSP(self, guess, matrix, vector1, vector2):
        """
        Set PETSC KSP solver.

        :arg guess: Boolean specifying if the iterative KSP solver initial guess is nonzero
        :arg matrix: PETSC matrix used by the KSP solver
        :arg vector1: PETSC vector corresponding to the initial values
        :arg vector2: PETSC vector corresponding to the new values
        :return vector2: PETSC vector of the new values
        """

        ksp = PETSc.KSP().create(PETSc.COMM_WORLD)
        if guess:
            ksp.setInitialGuessNonzero(guess)
        ksp.setOperators(matrix, matrix)
        ksp.setType("richardson")
        pc = ksp.getPC()
        pc.setType("bjacobi")
        ksp.setTolerances(rtol=self.rtol)
        ksp.solve(vector1, vector2)
        r = ksp.getConvergedReason()
        if r < 0:
            KSPReasons = self._make_reasons(PETSc.KSP.ConvergedReason())
            print(
                "LinearSolver failed to converge after %d iterations",
                ksp.getIterationNumber(),
                flush=True,
            )
            print("with reason: %s", KSPReasons[r], flush=True)
            raise RuntimeError("LinearSolver failed to converge!")
        ksp.destroy()

        return vector2

    def _buildFlowDirection(self, h1):
        """
        Build multiple flow direction based on neighbouring slopes.

        :arg h1: elevation array
        """

        t0 = clock()

        # Define multiple flow directions
        self.rcvID, self.slpRcv, self.distRcv, self.wghtVal = MFDreceivers(
            self.flowDir, self.inIDs, h1
        )

        # Account for marine regions
        self.seaID = np.where(h1 <= self.sealevel)[0]

        # Set deep ocean nodes
        self.rcvID0 = self.rcvID.copy()
        self.wghtVal0 = self.wghtVal.copy()

        deepID = np.where(h1 <= self.hbot)[0]
        self.rcvID0[deepID, :] = np.tile(deepID, (self.flowDir, 1)).T
        # self.distRcv0[deepID,:] = 0.
        self.wghtVal0[deepID, :] = 0.0

        # Set marine nodes
        self.rcvID[self.seaID, :] = np.tile(self.seaID, (self.flowDir, 1)).T
        self.distRcv[self.seaID, :] = 0.0
        self.wghtVal[self.seaID, :] = 0.0

        if MPIrank == 0 and self.verbose:
            print(
                "Flow Direction declaration (%0.02f seconds)" % (clock() - t0),
                flush=True,
            )

        return

    def sedChange(self):
        """
        Perform erosion deposition changes.

        This function contains methods for the following operations:

         - erosion/deposition induced by stream power law
         - depression identification and pit filling
         - stream induced deposition diffusion
         - hillslope diffusion
        """

        # Compute Erosion using Stream Power Law
        self._cptErosion()

        # Compute Deposition and Sediment Flux
        self._cptSedFlux()

        # Compute Sediment Deposition
        self._sedimentDeposition()

        # Compute Fresh Sediment Diffusion
        self._sedimentDiffusion()

        # Compute Hillslope Diffusion Law
        self._hillSlope()

    def FlowAccumulation(self):
        """
        Compute multiple flow accumulation.
        """

        self.dm.globalToLocal(self.hGlobal, self.hLocal)
        # Get global elevation for pit filling...
        t0 = clock()
        hl = self.hLocal.getArray().copy()
        gZ = np.zeros(self.gpoints)
        gZ = hl[self.lgIDs]
        gZ[self.outIDs] = -1.0e8
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gZ, op=MPI.MAX)

        # Perform pit filling
        nZ = fillPIT(self.sealevel + self.hbot, gZ)
        self._buildFlowDirection(nZ[self.glIDs])

        # Update fill elevation
        id = nZ < self.sealevel
        nZ[id] = gZ[id]
        self.pitID = np.where(nZ[self.glIDs] > hl)[0]
        self.FillL.setArray(nZ[self.glIDs])
        self.dm.localToGlobal(self.FillL, self.FillG)
        del hl, nZ, gZ, id
        if MPIrank == 0 and self.verbose:
            print("Compute pit filling (%0.02f seconds)" % (clock() - t0), flush=True)

        t0 = clock()
        # Build transport direction matrices
        WAMat = self.iMat.copy()
        WAMat0 = self.iMat.copy()
        indptr = np.arange(0, self.npoints + 1, dtype=PETSc.IntType)
        nodes = indptr[:-1]

        for k in range(0, self.flowDir):

            # Drainage area matrix
            tmpMat = self._matrix_build()
            data = -self.wghtVal[:, k].copy()
            data[self.rcvID[:, k].astype(PETSc.IntType) == nodes] = 0.0
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                self.rcvID[:, k].astype(PETSc.IntType),
                data,
                PETSc.InsertMode.INSERT_VALUES,
            )
            tmpMat.assemblyEnd()
            WAMat += tmpMat
            tmpMat.destroy()

            # Marine sediment matrix
            tmpMat0 = self._matrix_build()
            data0 = -self.wghtVal0[:, k].copy()
            data0[self.rcvID0[:, k].astype(PETSc.IntType) == nodes] = 0.0
            tmpMat0.assemblyBegin()
            tmpMat0.setValuesLocalCSR(
                indptr,
                self.rcvID0[:, k].astype(PETSc.IntType),
                data0,
                PETSc.InsertMode.INSERT_VALUES,
            )
            tmpMat0.assemblyEnd()
            WAMat0 += tmpMat0
            tmpMat0.destroy()

        del data, data0, indptr, nodes
        gc.collect()

        # Solve flow accumulation
        self.wMat = WAMat.transpose().copy()
        self.wMat0 = WAMat0.transpose().copy()
        if self.tNow == self.rStart:
            self._solve_KSP(False, self.wMat, self.bG, self.FAG)
        else:
            self._solve_KSP(True, self.wMat, self.bG, self.FAG)
        WAMat.destroy()
        WAMat0.destroy()

        self.dm.globalToLocal(self.FAG, self.FAL, 1)

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Flow Accumulation (%0.02f seconds)" % (clock() - t0),
                flush=True,
            )

        return

    def _getErosionRate(self):
        """
        Compute sediment and bedrock erosion rates.
        """

        Kcoeff = self.FAL.getArray()
        Kbr = np.sqrt(Kcoeff) * self.K * self.dt
        Kbr[self.seaID] = 0.0
        Kbr[self.pitID] = 0.0

        # Initialise identity matrices...
        EbedMat = self.iMat.copy()
        wght = self.wghtVal.copy()

        # Define erosion coefficients
        for k in range(0, self.flowDir):

            indptr = np.arange(0, self.npoints + 1, dtype=PETSc.IntType)
            nodes = indptr[:-1]
            # Define erosion limiter to prevent formation of flat
            dh = self.hOldArray - self.hOldArray[self.rcvID[:, k]]
            limiter = np.divide(dh, dh + 1.0e-3, out=np.zeros_like(dh), where=dh != 0)

            # Bedrock erosion processes SPL computation (maximum bedrock incision)
            data = np.divide(
                Kbr * limiter,
                self.distRcv[:, k],
                out=np.zeros_like(Kcoeff),
                where=self.distRcv[:, k] != 0,
            )
            tmpMat = self._matrix_build()
            wght[self.seaID, k] = 0.0
            data = np.multiply(data, -wght[:, k])
            data[self.rcvID[:, k].astype(PETSc.IntType) == nodes] = 0.0
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                self.rcvID[:, k].astype(PETSc.IntType),
                data,
                PETSc.InsertMode.INSERT_VALUES,
            )
            tmpMat.assemblyEnd()
            EbedMat += tmpMat
            EbedMat -= self._matrix_build_diag(data)
            tmpMat.destroy()

        del dh, limiter, wght, data
        gc.collect()

        # Solve bedrock erosion thickness
        self._solve_KSP(True, EbedMat, self.hOld, self.vGlob)
        EbedMat.destroy()
        self.stepED.waxpy(-1.0, self.hOld, self.vGlob)

        # Define erosion rate (positive for incision)
        E = -self.stepED.getArray().copy()
        E = np.divide(E, self.dt)
        E[E < 0.0] = 0.0
        self.Eb.setArray(E)
        self.dm.globalToLocal(self.Eb, self.EbLocal, 1)
        E = self.EbLocal.getArray().copy()
        E[self.seaID] = 0.0
        E[self.pitID] = 0.0
        self.EbLocal.setArray(E)
        self.dm.localToGlobal(self.EbLocal, self.Eb, 1)

        del E, Kcoeff, Kbr
        gc.collect()

        return

    def _cptErosion(self):
        """
        Compute erosion using stream power law.
        """

        t0 = clock()

        # Constant local & global vectors/arrays
        self.Eb.set(0.0)
        self.hGlobal.copy(result=self.hOld)
        self.dm.globalToLocal(self.hOld, self.hOldLocal, 1)
        self.hOldArray = self.hOldLocal.getArray().copy()

        self._getErosionRate()

        # Update bedrock thicknesses due to erosion
        Eb = self.Eb.getArray().copy()
        self.stepED.setArray(-Eb * self.dt)
        self.cumED.axpy(1.0, self.stepED)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)

        self.hGlobal.axpy(1.0, self.stepED)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

        if MPIrank == 0 and self.verbose:
            print(
                "Get Erosion Thicknesses (%0.02f seconds)" % (clock() - t0), flush=True
            )

        return

    def _cptSedFlux(self):
        """
        Compute sediment flux.
        """

        # Build sediment load matrix
        t0 = clock()
        SLMat = self.wMat.copy()
        SLMat -= self.iMat
        SLMat.scale(1.0 - self.wgth)
        SLMat += self.iMat
        Eb = self.Eb.getArray().copy()
        Eb = np.multiply(Eb, 1.0 - self.frac_fine)

        self.stepED.setArray(Eb)
        self.stepED.pointwiseMult(self.stepED, self.areaGlobal)
        if self.tNow == self.rStart:
            self._solve_KSP(False, SLMat, self.stepED, self.vSed)
        else:
            self._solve_KSP(True, SLMat, self.stepED, self.vSed)
        SLMat.destroy()

        # Update local vector
        self.dm.globalToLocal(self.vSed, self.vSedLocal, 1)
        if MPIrank == 0 and self.verbose:
            print("Update Sediment Load (%0.02f seconds)" % (clock() - t0), flush=True)
        del Eb
        gc.collect()

        return

    def _hillSlope(self):
        """
        Perform hillslope diffusion.
        """

        t0 = clock()
        if self.Cd > 0.0:
            # Get erosion values for considered time step
            self.hGlobal.copy(result=self.hOld)
            self._solve_KSP(True, self.Diff, self.hOld, self.hGlobal)

            # Update cumulative erosion/deposition and elevation
            self.stepED.waxpy(-1.0, self.hOld, self.hGlobal)
            self.cumED.axpy(1.0, self.stepED)
            self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
            self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Hillslope Processes (%0.02f seconds)" % (clock() - t0),
                flush=True,
            )

        return

    def _sedimentDeposition(self):
        """
        Perform sediment deposition from incoming river flux.
        """

        t0 = clock()
        # Get the marine volumetric sediment rate (m3/yr) to diffuse
        # during the time step...
        tmp = self.vSedLocal.getArray().copy()
        depVol = np.zeros(self.npoints)
        depVol[self.seaID] = tmp[self.seaID] * self.dt  # volume m3

        # Build sediment load matrix
        h0 = self.hLocal.getArray().copy()
        hf = self.FillL.getArray().copy()
        Db = np.zeros(self.npoints)
        Db[self.seaID] = -0.9 * (self.sealevel - h0[self.seaID])
        Db[Db > 0] = 0.0
        Db = np.multiply(Db, self.area)
        Db += depVol
        self.tmpL.setArray(Db)
        self.dm.localToGlobal(self.tmpL, self.stepED)
        del tmp, depVol, Db
        gc.collect()

        self._solve_KSP(False, self.wMat0, self.stepED, self.tmp)
        self.tmp.pointwiseDivide(self.tmp, self.areaGlobal)
        self.dm.globalToLocal(self.tmp, self.tmpL, 1)
        tmp = self.tmpL.getArray().copy()

        id = tmp >= 0.0
        depo = np.zeros(self.npoints)
        depo[id] = 0.9 * (self.sealevel - h0[id])
        depo[depo < 0] = 0.0
        tmp[hf > self.sealevel] = -1.0e6
        tmp[id] = -1.0e6

        id = 0.9 * (self.sealevel - h0) + tmp > self.rtol
        depo[id] += 0.9 * (self.sealevel - h0[id]) + tmp[id]
        self.tmpL.setArray(depo)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        del depo, tmp, h0, hf
        gc.collect()

        # Get deposition thickness
        self.stepED.set(0.0)
        if self.wgth > 0.0:
            self.stepED.axpy(-(1.0 - self.frac_fine), self.Eb)
            self.stepED.pointwiseMult(self.stepED, self.areaGlobal)
            self.stepED.axpy(1.0, self.vSed)
            self.stepED.scale(self.wgth * self.dt / (1.0 - self.wgth))
            self.stepED.pointwiseDivide(self.stepED, self.areaGlobal)

        self.stepED.axpy(1.0, self.tmp)
        self.cumED.axpy(1.0, self.stepED)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)

        self.hGlobal.axpy(1.0, self.stepED)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

        if MPIrank == 0 and self.verbose:
            print(
                "Perform Sediment Deposition (%0.02f seconds)" % (clock() - t0),
                flush=True,
            )

        return

    def _sedimentDiffusion(self):
        """
        Perform freshly deposited sediment diffusion.
        """

        t0 = clock()
        limit = 1.0e-1
        h0 = self.hOldLocal.getArray().copy()

        h = self.hLocal.getArray().copy()
        dh = h - h0
        dh[dh < 0.0] = 0.0
        sedCoeffs = setDiffusionCoeff(self.sedimentK * self.dt, limit, h, dh)
        sedDiff = self._matrix_build_diag(sedCoeffs[:, 0])
        del h, dh
        gc.collect()

        for k in range(0, self.maxnb):
            tmpMat = self._matrix_build()
            indptr = np.arange(0, self.npoints + 1, dtype=PETSc.IntType)
            indices = self.FVmesh_ngbID[:, k].copy()
            data = np.zeros(self.npoints)
            ids = np.nonzero(indices < 0)
            indices[ids] = ids
            data = sedCoeffs[:, k + 1]
            ids = np.nonzero(data == 0.0)
            indices[ids] = ids
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                indices.astype(PETSc.IntType),
                data,
                PETSc.InsertMode.INSERT_VALUES,
            )
            tmpMat.assemblyEnd()
            sedDiff += tmpMat
            tmpMat.destroy()
        del ids, indices, indptr
        gc.collect()

        # Get erosion values for considered time step
        self.hGlobal.copy(result=self.hOld)
        self._solve_KSP(True, sedDiff, self.hOld, self.hGlobal)
        sedDiff.destroy()

        # Update cumulative erosion/deposition and elevation
        self.stepED.waxpy(-1.0, self.hOld, self.hGlobal)
        self.cumED.axpy(1.0, self.stepED)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

        if MPIrank == 0 and self.verbose:
            print("Diffuse Top Sediment (%0.02f seconds)" % (clock() - t0), flush=True)

        return
