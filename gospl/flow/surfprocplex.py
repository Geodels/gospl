import os
import gc
import sys
import vtk
import warnings
import petsc4py
import numpy as np
from sklearn.neighbors import BallTree


from mpi4py import MPI
from time import process_time
from vtk.util import numpy_support

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import fillPIT
    from gospl._fortran import setMaxNb
    from gospl._fortran import MFDreceivers
    from gospl._fortran import sedReceivers
    from gospl._fortran import setHillslopeCoeff
    from gospl._fortran import setDiffusionCoeff

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD


class SPMesh(object):
    """
    Performing surface evolution induced by considered processes:

     - hillslope processes
     - rivers using stream power law
    """

    def __init__(self, *args, **kwargs):

        # KSP solver parameters
        self.rtol = 1.0e-8

        # Identity matrix construction
        self.II = np.arange(0, self.npoints + 1, dtype=petsc4py.PETSc.IntType)
        self.JJ = np.arange(0, self.npoints, dtype=petsc4py.PETSc.IntType)
        self.iMat = self._matrix_build_diag(np.ones(self.npoints))
        self.iMat0 = self._matrix_build_diag(np.zeros(self.npoints))

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
        self.Qs = self.hGlobal.duplicate()
        self.Eb = self.hGlobal.duplicate()
        self.EbLocal = self.hLocal.duplicate()
        self.vSed = self.hGlobal.duplicate()
        self.vSedLocal = self.hLocal.duplicate()

        self.maxnb = setMaxNb(self.npoints)

        return

    def _matrix_build(self, nnz=(1, 1)):
        """
        Define PETSC Matrix.

        :arg nnz: array containing the number of nonzero blocks
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
        Define PETSC Diagonal Matrix.

        :arg V: diagonal data array
        :arg nnz: array containing the number of nonzero blocks
        """

        matrix = self._matrix_build()

        # Define diagonal matrix
        matrix.assemblyBegin()
        matrix.setValuesLocalCSR(
            self.II, self.JJ, V, petsc4py.PETSc.InsertMode.INSERT_VALUES
        )
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

        ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
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
            KSPReasons = self._make_reasons(petsc4py.PETSc.KSP.ConvergedReason())
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

        t0 = process_time()

        # Define multiple flow directions
        self.rcvID, self.distRcv, self.wghtVal, self.maxSlp = MFDreceivers(
            self.flowDir, self.inIDs, h1, self.sealevel
        )

        # Get marine regions
        self.seaID = np.where(h1 <= self.sealevel)[0]

        if MPIrank == 0 and self.verbose:
            print(
                "Flow Direction declaration (%0.02f seconds)" % (process_time() - t0),
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
        # for k in range(1):
        #     self._sedimentDiffusion(self.dt * 1.0)

        # Compute Hillslope Diffusion Law
        self._hillSlope()

        return

    def _getCoastlines(self, data):
        """
        Update the elevation and set it up in the VTK mesh. Then perform contour filtering
        to extract coastline positions globally.
        """

        pointData = self.vtkMesh.GetPointData()
        array = numpy_support.numpy_to_vtk(data, deep=1)
        array.SetName("z")
        pointData.AddArray(array)

        cf = vtk.vtkContourFilter()
        cf.SetInputData(self.vtkMesh)
        cf.SetValue(0, self.sealevel)
        cf.SetInputArrayToProcess(0, 0, 0, 0, "z")
        cf.GenerateTrianglesOff()
        cf.Update()
        coasts = numpy_support.vtk_to_numpy(cf.GetOutput().GetPoints().GetData())

        # Convert from cartesian to lat / lon (radians)
        coastLatLon = np.zeros((len(coasts), 2), dtype=np.float64)
        coastLatLon[:, 0] = np.arcsin(coasts[:, 2] / self.radius)
        coastLatLon[:, 1] = np.arctan2(coasts[:, 1], coasts[:, 0])

        del coasts

        return coastLatLon

    def _dist2coast(self, oceanPts, coastPts, k_neighbors=1):
        """
        Find nearest neighbors for all source points from a set of candidate points.
        """

        # Create tree from the candidate points
        tree = BallTree(coastPts, leaf_size=15, metric="haversine")

        # Find closest points and distances
        distances, indices = tree.query(oceanPts, k=k_neighbors)
        del tree, indices

        return distances[:, 0] * self.radius

    def flowAccumulation(self):
        """
        Compute multiple flow accumulation.
        """

        self.dm.globalToLocal(self.hGlobal, self.hLocal)
        # Get global elevation for pit filling...
        t0 = process_time()
        hl = self.hLocal.getArray().copy()
        gZ = hl[self.lgIDs]
        gZ[self.outIDs] = -1.0e8
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gZ, op=MPI.MAX)

        # Perform pit filling on process rank 0
        if MPIrank == 0:
            nZ = fillPIT(self.sealevel - 5000.0, gZ)
        else:
            nZ = np.zeros(self.gpoints, dtype=np.float64)
        nZ = MPI.COMM_WORLD.bcast(nZ, root=0)

        if MPIrank == 0 and self.verbose:
            print(
                "Compute pit filling (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        # Build flow direction
        t0 = process_time()
        self._buildFlowDirection(nZ[self.glIDs])

        # Get marine weight coefficients
        self.wgtCoeffs = sedReceivers(nZ[self.glIDs], self.sealevel)

        # Define coastal points
        with warnings.catch_warnings():
            # ignore all caught warnings
            warnings.filterwarnings("ignore")
            coastLatLon = self._getCoastlines(gZ)
        self.coastDist = np.zeros(self.npoints)
        self.coastDist[self.seaID] = self._dist2coast(
            self.lLatLon[self.seaID, :], coastLatLon
        )

        # Update fill elevation
        id = nZ < self.sealevel
        nZ[id] = gZ[id]
        self.pitID = np.where(nZ[self.glIDs] > hl)[0]
        self.FillL.setArray(nZ[self.glIDs])
        self.dm.localToGlobal(self.FillL, self.FillG)
        del hl, nZ, gZ, id, coastLatLon
        gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "Build flow directions and coastlines (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        t0 = process_time()
        # Build transport direction matrices
        WAMat = self.iMat.copy()
        indptr = np.arange(0, self.npoints + 1, dtype=petsc4py.PETSc.IntType)
        nodes = indptr[:-1]

        for k in range(0, self.flowDir):

            # Drainage area matrix
            tmpMat = self._matrix_build()
            data = -self.wghtVal[:, k].copy()
            data[self.rcvID[:, k].astype(petsc4py.PETSc.IntType) == nodes] = 0.0
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                self.rcvID[:, k].astype(petsc4py.PETSc.IntType),
                data,
                petsc4py.PETSc.InsertMode.INSERT_VALUES,
            )
            tmpMat.assemblyEnd()
            WAMat += tmpMat
            tmpMat.destroy()

        del data, indptr, nodes
        gc.collect()

        # Solve flow accumulation
        self.wMat = WAMat.transpose().copy()
        if self.tNow == self.rStart:
            self._solve_KSP(False, self.wMat, self.bG, self.FAG)
        else:
            self._solve_KSP(True, self.wMat, self.bG, self.FAG)
        WAMat.destroy()

        self.dm.globalToLocal(self.FAG, self.FAL, 1)

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Flow Accumulation (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def _getErosionRate(self):
        """
        Compute sediment and bedrock erosion rates in metres per year.
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

            indptr = np.arange(0, self.npoints + 1, dtype=petsc4py.PETSc.IntType)
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
            data[self.rcvID[:, k].astype(petsc4py.PETSc.IntType) == nodes] = 0.0
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                self.rcvID[:, k].astype(petsc4py.PETSc.IntType),
                data,
                petsc4py.PETSc.InsertMode.INSERT_VALUES,
            )
            tmpMat.assemblyEnd()
            EbedMat += tmpMat
            EbedMat -= self._matrix_build_diag(data)
            tmpMat.destroy()

        del dh, limiter, wght, data
        gc.collect()

        # Solve bedrock erosion thickness
        self._solve_KSP(True, EbedMat, self.hOld, self.tmp)
        EbedMat.destroy()
        self.stepED.waxpy(-1.0, self.hOld, self.tmp)

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

        t0 = process_time()

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
                "Get Erosion Thicknesses (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def _cptSedFlux(self):
        """
        Compute sediment flux in cubic metres per year.
        """

        # Build sediment load matrix
        t0 = process_time()
        SLMat = self.wMat.copy()
        SLMat -= self.iMat
        SLMat.scale(1.0 - self.wgth)
        SLMat += self.iMat

        # Get erosion rate (m/yr)
        Eb = self.Eb.getArray().copy()
        Eb = np.multiply(Eb, 1.0 - self.frac_fine)
        self.stepED.setArray(Eb)
        self.stepED.pointwiseMult(self.stepED, self.areaGlobal)

        # Get the volume of sediment transported in m3 per year
        if self.tNow == self.rStart:
            self._solve_KSP(False, SLMat, self.stepED, self.vSed)
        else:
            self._solve_KSP(True, SLMat, self.stepED, self.vSed)
        SLMat.destroy()

        # Update local vector
        self.dm.globalToLocal(self.vSed, self.vSedLocal, 1)
        if MPIrank == 0 and self.verbose:
            print(
                "Update Sediment Load (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )
        del Eb
        gc.collect()

        return

    def _hillSlope(self):
        """
        Perform hillslope diffusion.
        """

        if self.Cda == 0.0 and self.Cdm == 0.0:
            return

        t0 = process_time()

        # Diffusion matrix construction
        Cd = np.full(self.npoints, self.Cda, dtype=np.float64)
        Cd[self.seaID] = self.Cdm

        diffCoeffs = setHillslopeCoeff(self.npoints, Cd * self.dt)
        self.Diff = self._matrix_build_diag(diffCoeffs[:, 0])
        indptr = np.arange(0, self.npoints + 1, dtype=petsc4py.PETSc.IntType)

        for k in range(0, self.maxnb):
            tmpMat = self._matrix_build()
            indices = self.FVmesh_ngbID[:, k].copy()
            data = diffCoeffs[:, k + 1]
            ids = np.nonzero(data == 0.0)
            indices[ids] = ids
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                indices.astype(petsc4py.PETSc.IntType),
                data,
                petsc4py.PETSc.InsertMode.INSERT_VALUES,
            )
            tmpMat.assemblyEnd()
            self.Diff += tmpMat
            tmpMat.destroy()
        del ids, indices, indptr, diffCoeffs, Cd
        gc.collect()

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
                "Compute Hillslope Processes (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def _sedimentDeposition(self):
        """
        Perform sediment deposition from incoming river flux.
        """

        t0 = process_time()

        # Get the marine volumetric sediment rate (m3 / yr) to diffuse
        # during the time step...
        tmp = self.vSedLocal.getArray().copy()
        Qs = np.zeros(self.npoints, dtype=np.float64)

        # Convert in volume (m3) for considered timestep
        Qs[self.seaID] = tmp[self.seaID] * self.dt

        if np.sum(Qs) == 0:
            del Qs, tmp
            gc.collect()

            return

        # Build sediment fluxes vector
        self.tmpL.setArray(Qs)
        self.dm.localToGlobal(self.tmpL, self.Qs)

        # Sediment fluxes weight matrix construction
        wgtQs = self._matrix_build_diag(self.wgtCoeffs[:, 0])
        indptr = np.arange(0, self.npoints + 1, dtype=petsc4py.PETSc.IntType)

        for k in range(0, self.maxnb):
            tmpMat = self._matrix_build()
            indices = self.FVmesh_ngbID[:, k].copy()
            data = self.wgtCoeffs[:, k + 1]
            ids = np.nonzero(data == 0.0)
            indices[ids] = ids
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                indices.astype(petsc4py.PETSc.IntType),
                data,
                petsc4py.PETSc.InsertMode.INSERT_VALUES,
            )
            tmpMat.assemblyEnd()
            wgtQs += tmpMat
            tmpMat.destroy()

        self.wgtQs = wgtQs.transpose().copy()
        wgtQs.destroy()
        del indices, data, ids, indptr
        gc.collect()

        # Build maximum sediment volume (in m3)
        # Here we fix the slope for the shelf (-5.e-4)
        h0 = self.hLocal.getArray().copy()
        toplimit = self.sealevel - self.coastDist * 5.0e-4
        maxD = 0.9 * (toplimit - h0) * self.area
        maxD[maxD < 0.0] = 0.0

        # Initialise local deposit volume
        Db = np.zeros(self.npoints, dtype=np.float64)

        while self.Qs.sum() > 0:

            # Check local cell deposits
            excess = np.where(np.logical_and(Qs > 0, Qs >= maxD))[0]
            added = np.where(np.logical_and(Qs > 0, Qs < maxD))[0]
            # Overfilling cells
            Db[excess] += maxD[excess]
            Qs[excess] -= maxD[excess]
            maxD[excess] = 0.0
            # Underfilling cells
            Db[added] += Qs[added]
            maxD[added] -= Qs[added]
            Qs[added] = 0.0

            # Distribute overfilling fluxes downstream
            self.tmpL.setArray(Qs)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.wgtQs.mult(self.tmp, self.Qs)

            # Update downstream nodes receiving fluxes
            self.dm.globalToLocal(self.Qs, self.tmpL)
            Qs = self.tmpL.getArray().copy()

        # Get deposition thickness
        self.stepED.set(0.0)
        if self.wgth > 0.0:
            self.stepED.axpy(-(1.0 - self.frac_fine), self.Eb)
            self.stepED.pointwiseMult(self.stepED, self.areaGlobal)
            self.stepED.axpy(1.0, self.vSed)
            self.stepED.scale(self.wgth * self.dt / (1.0 - self.wgth))
            self.stepED.pointwiseDivide(self.stepED, self.areaGlobal)

        self.tmpL.setArray(Db / self.area)
        self.dm.localToGlobal(self.tmpL, self.tmp)

        self.stepED.axpy(1.0, self.tmp)
        self.cumED.axpy(1.0, self.stepED)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)

        self.hGlobal.axpy(1.0, self.stepED)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

        del Qs, tmp, h0, maxD, Db, excess, added
        gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "Perform Sediment Deposition (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def _sedimentDiffusion(self, dt):
        """
        Perform freshly deposited sediment diffusion.
        """

        t0 = process_time()
        limit = 1.0
        h0 = self.hOldLocal.getArray().copy()

        h = self.hLocal.getArray().copy()
        dh = h - h0
        dh[dh < 0.0] = 0.0

        # Get local sediment flux (m/yr)
        Qs = self.tmpL.getArray().copy()

        sedCoeffs = setDiffusionCoeff(self.sedimentK * dt, limit, h, Qs)
        sedDiff = self._matrix_build_diag(sedCoeffs[:, 0])
        del h, dh
        gc.collect()

        for k in range(0, self.maxnb):
            tmpMat = self._matrix_build()
            indptr = np.arange(0, self.npoints + 1, dtype=petsc4py.PETSc.IntType)
            indices = self.FVmesh_ngbID[:, k].copy()
            ids = np.nonzero(indices < 0)
            indices[ids] = ids
            data = sedCoeffs[:, k + 1]
            ids = np.nonzero(data == 0.0)
            indices[ids] = ids
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                indices.astype(petsc4py.PETSc.IntType),
                data,
                petsc4py.PETSc.InsertMode.INSERT_VALUES,
            )
            tmpMat.assemblyEnd()
            sedDiff += tmpMat
            tmpMat.destroy()
        del ids, indices, indptr, sedCoeffs, data
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
            print(
                "Diffuse Top Sediment (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return
