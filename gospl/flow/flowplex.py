import os
import gc
import sys
import vtk
import warnings
import petsc4py
import numpy as np
from scipy import spatial
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time
from vtk.util import numpy_support

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import fillPIT
    from gospl._fortran import MFDreceivers

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD


class FAMesh(object):
    """
    Compute flow over surface and river erosion.
    """

    def __init__(self, *args, **kwargs):

        # KSP solver parameters
        self.rtol = 1.0e-8

        # Identity matrix construction
        self.II = np.arange(0, self.npoints + 1, dtype=petsc4py.PETSc.IntType)
        self.JJ = np.arange(0, self.npoints, dtype=petsc4py.PETSc.IntType)
        self.iMat = self._matrix_build_diag(np.ones(self.npoints))

        # Petsc vectors
        self.FAG = self.hGlobal.duplicate()
        self.FAL = self.hLocal.duplicate()
        self.fillFAG = self.hGlobal.duplicate()
        self.fillFAL = self.hLocal.duplicate()
        self.hOld = self.hGlobal.duplicate()
        self.hOldLocal = self.hLocal.duplicate()
        self.Eb = self.hGlobal.duplicate()
        self.stepED = self.hGlobal.duplicate()
        self.EbLocal = self.hLocal.duplicate()

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

    def _buildFlowDirection(self, h):
        """
        Build multiple flow direction based on neighbouring slopes.

        :arg h1: elevation array
        """

        t0 = process_time()

        sl = self.sealevel
        if self.isfill:
            sl = self.sealevel - 1000.0

        # Define multiple flow directions for unfilled elevation
        self.rcvID, self.distRcv, self.wghtVal = MFDreceivers(
            self.flowDir, self.inIDs, h, sl
        )

        if not self.isfill:
            # Get nodes that have no receivers
            sum_weight = np.sum(self.wghtVal, axis=1)
            self.pitPts = sum_weight == 0.0

        # Get marine regions
        self.seaID = np.where(h <= self.sealevel)[0]

        # Set marine nodes
        self.rcvID[self.seaID, :] = np.tile(self.seaID, (self.flowDir, 1)).T
        self.distRcv[self.seaID, :] = 0.0
        self.wghtVal[self.seaID, :] = 0.0

        if MPIrank == 0 and self.verbose:
            print(
                "Flow Direction declaration (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def _distanceCoasts(self, data, k_neighbors=1):
        """
        Update the elevation and set it up in the VTK mesh. Then perform contour filtering
        to extract coastline positions globally and compute distance to shore for ocean nodes.
        """

        t0 = process_time()

        self.coastDist = np.zeros(self.npoints)
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
        coastXYZ = numpy_support.vtk_to_numpy(cf.GetOutput().GetPoints().GetData())
        tree = spatial.cKDTree(coastXYZ, leafsize=10)
        self.coastDist[self.seaID], indices = tree.query(
            self.lcoords[self.seaID, :], k=k_neighbors
        )
        del array, pointData, cf, tree, indices, coastXYZ
        gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "Construct distance to coast (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def _pitInformation(self, gZ, hFill):
        """
        Define depression characteristics.
        """

        # Compute pit volumes
        groupPits = npi.group_by(self.pits[:, 0])
        pitNb, self.pitVol = groupPits.sum((hFill - gZ) * self.garea)
        _, outids, _ = np.intersect1d(self.pits[:, 0], pitNb, return_indices=True)
        self.outFlows = self.pits[outids, 1]

        # For pits that are below a given volume threshold, we impose the filled
        # elevation value (here we set the threshold to 10m * min area)
        id = np.where(self.pitVol / self.minArea[1] < 10.0)[0]
        self.pitVol[id] = 0.0
        mask = np.in1d(self.pits[:, 0], id)
        gZ[mask] = hFill[mask]

        del groupPits, pitNb, id, mask, outids
        gc.collect()

        return gZ

    def _depressionlessSurface(self):
        """
        Compute depression less surface.
        """

        # Get global elevations for pit filling...
        t0 = process_time()
        hl = self.hLocal.getArray().copy()
        gZ = hl[self.lgIDs]
        gZ[self.outIDs] = -1.0e8
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gZ, op=MPI.MAX)

        if self.isfill:
            # Perform pit filling on process rank 0
            if MPIrank == 0:
                # Fillpit returns:
                # - hFill: filled elevation values
                # - pits: 2D array containing each pit ID and
                #         corresponding overspilling point ID
                hFill, pits = fillPIT(self.sealevel - 1000.0, gZ)
            else:
                hFill = np.zeros(self.gpoints, dtype=np.float64)
                pits = np.zeros((self.gpoints, 2), dtype=np.int64)
            hFill = MPI.COMM_WORLD.bcast(hFill, root=0)
            self.pits = MPI.COMM_WORLD.bcast(pits, root=0)
            self.pitID = np.where(self.pits[self.glIDs, 0] >= 0)[0]
            self.hFill = hFill[self.glIDs]

            # Get depressions information
            gZ = self._pitInformation(gZ, hFill)

            if MPIrank == 0 and self.verbose:
                print(
                    "Get pit filling information (%0.02f seconds)"
                    % (process_time() - t0),
                    flush=True,
                )

            del hl, pits
            gc.collect()
        else:
            hFill = None

            del hl
            gc.collect()

        return gZ, hFill

    def _matrixFA(self):
        """
        Build transport direction matrices for both filled and unfilled elevations.
        """

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

        if self.isfill:
            self.fillMat = self.wMat.copy()

        WAMat.destroy()

        return

    def flowAccumulation(self, filled=False):
        """
        Compute multiple flow accumulation.

         - depression identification and pit filling
         - flow accumulation based on filled and unfilled surfaces
        """

        self.isfill = filled
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        # Fill surface to remove pits
        gZ, hFill = self._depressionlessSurface()

        # Build flow direction
        if self.isfill:
            self._buildFlowDirection(hFill[self.glIDs])
            # Define coastal distance for marine points
            if self.vtkMesh is not None:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore")
                    self._distanceCoasts(gZ)
            del hFill, gZ
            gc.collect()

        else:
            self._buildFlowDirection(gZ[self.glIDs])
            del gZ
            gc.collect()

        # Build transport direction matrices
        t0 = process_time()
        self._matrixFA()

        if self.isfill:
            # Solve flow accumulation for filled elevation
            self._solve_KSP(True, self.wMat, self.bG, self.fillFAG)
            self.dm.globalToLocal(self.fillFAG, self.fillFAL, 1)
        else:
            # Solve flow accumulation for unfilled elevation
            if self.tNow == self.rStart:
                self._solve_KSP(False, self.wMat, self.bG, self.FAG)
            else:
                self._solve_KSP(True, self.wMat, self.bG, self.FAG)
            self.dm.globalToLocal(self.FAG, self.FAL, 1)

        if MPIrank == 0 and self.verbose:
            if self.isfill:
                print(
                    "Compute Filled Flow Accumulation (%0.02f seconds)"
                    % (process_time() - t0),
                    flush=True,
                )
            else:
                print(
                    "Compute Flow Accumulation (%0.02f seconds)"
                    % (process_time() - t0),
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
        self._solve_KSP(True, EbedMat, self.hOld, self.stepED)
        EbedMat.destroy()
        self.tmp.waxpy(-1.0, self.hOld, self.stepED)

        # Define erosion rate (positive for incision)
        E = -self.tmp.getArray().copy()
        E = np.divide(E, self.dt)
        E[E < 0.0] = 0.0
        self.Eb.setArray(E)
        self.dm.globalToLocal(self.Eb, self.EbLocal, 1)
        E = self.EbLocal.getArray().copy()
        E[self.seaID] = 0.0
        ids = np.where(
            np.logical_and(
                self.hOldArray > self.sealevel + 1.0e-2,
                self.hOldArray - E * self.dt < self.sealevel + 1.0e-2,
            )
        )[0]
        E[ids] = (self.hOldArray[ids] - self.sealevel - 1.0e-2) / self.dt
        self.EbLocal.setArray(E)
        self.dm.localToGlobal(self.EbLocal, self.Eb, 1)

        del E, Kcoeff, Kbr, ids
        gc.collect()

        return

    def riverIncision(self):
        """
        Compute stream erosion using stream power law.
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
        self.tmp.setArray(-Eb * self.dt)
        self.cumED.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)

        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

        # Update stratigraphic layer
        if self.stratNb > 0:
            self.erodeStrat()

        if MPIrank == 0 and self.verbose:
            print(
                "Get Erosion Thicknesses (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return
