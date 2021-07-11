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
    from gospl._fortran import fillpit
    from gospl._fortran import mfdreceivers

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD


class FAMesh(object):
    """
    This class calculates **drainage area** in an implicit, iterative manner using PETSc solvers. It accounts  for multiple flow direction paths (SFD to MFD) based on user input declaration.

    .. note::

        The class follows the parallel approach described in `Richardson et al., 2014 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013WR014326>`_ where the iterative nature of the computational algorithms used to solve the linear system creates the possibility of accelerating the solution by providing an initial guess.

    For drainage computation, the class requires to compute depressionless surfaces and the *priority-flood + ϵ* variant of the algorithm proposed in `Barnes et al. (2014) <https://doi.org/10.1016/j.cageo.2013.04.024>`_ is used. It provides a solution to remove flat surfaces, and it produces surfaces for which each cell has a defined gradient from which flow directions can be determined.

    Finally, the class computes river incision expressed using a **stream power formulation** function of river discharge and slope.

    """

    def __init__(self, *args, **kwargs):
        """
        The initialisation of `FAMesh` class consists in the declaration of PETSc vectors and matrices.
        """

        # KSP solver parameters
        self.rtol = 1.0e-8

        # Identity matrix construction
        self.II = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        self.JJ = np.arange(0, self.lpoints, dtype=petsc4py.PETSc.IntType)
        self.iMat = self._matrix_build_diag(np.ones(self.lpoints))

        # Empty matrix construction
        self.II = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        self.JJ = np.arange(0, self.lpoints, dtype=petsc4py.PETSc.IntType)
        self.zMat = self._matrix_build_diag(np.zeros(self.lpoints))

        # Petsc vectors
        # self.fillFAG = self.hGlobal.duplicate()
        self.fillFAL = self.hLocal.duplicate()
        self.FAG = self.hGlobal.duplicate()
        self.FAL = self.hLocal.duplicate()
        self.hOld = self.hGlobal.duplicate()
        self.hOldLocal = self.hLocal.duplicate()
        self.Eb = self.hGlobal.duplicate()
        self.stepED = self.hGlobal.duplicate()
        self.EbLocal = self.hLocal.duplicate()

        # Clinoforms slopes for each sediment type
        self.slps = [6.0, 5.0, 3.0, 1.0]

        return

    def _matrix_build(self, nnz=(1, 1)):
        """
        Creates a sparse PETSc matrix.

        .. note::

            To achieve good performance during matrix assembly, the function preallocates the matrix storage by setting the array nnz.

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

    def _solve_KSP(self, guess, matrix, vector1, vector2):
        """

        PETSc *scalable linear equations solvers* (**KSP**) component provides Krylov subspace iterative method and a preconditioner. Here, flow accumulation solution is obtained using PETSc Richardson solver (`richardson`) with block Jacobian preconditioning (`bjacobi`).

        .. note::

            The solver choice was made based on the convergence results from `Richardson et al. (2014) <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013WR014326>`_ but can be changed if better solver and preconditioner combinations are found.

        Using such iterative method allows for an initial guess to be provided. When this initial guess is close to the solution, the number of iterations required for convergence dramatically decreases. Here the flow discharge solution from previous time step can be passed as an initial `guess` to the solver as discharge often exhibits little change between successive time intervals.

        :arg guess: Boolean specifying if the iterative KSP solver initial guess is nonzero (when provided it corresponds to the previous flow discharge values)
        :arg matrix: PETSc sparse matrix used by the KSP solver composed of diagonal terms set to unity (identity matrix) and off-diagonal terms (weights between 0 and 1). The weights are calculated based on the number of downslope neighbours (based on the chosen number of flow direction directions) and are proportional to the slope.
        :arg vector1: PETSc vector corresponding to the local volume of water available for runoff during a given time step (*e.g.* voronoi area times local precipitation rate)
        :arg vector2: PETSc vector corresponding to the unknown flow discharge values

        :return: vector2 PETSc vector of the new flow discharge values
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
        This function builds from neighbouring slopes the flow directions. It calls a fortran subroutine that locally computes for each vertice:

        - the indices of receivers (downstream) nodes depending on the desired number of flow directions (SFD to MFD).
        - the distances to the receivers based on mesh resolution.
        - the associated weights calculated based on the number of receivers and proportional to the slope.

        :arg h: elevation numpy array
        """

        t0 = process_time()

        # Define multiple flow directions for unfilled elevation
        self.rcvID, self.distRcv, self.wghtVal = mfdreceivers(
            self.flowDir, self.flowExp, self.inIDs, h, self.sealevel
        )
        # Get marine regions
        self.seaID = np.where(h <= self.sealevel)[0]

        # Set marine nodes
        if self.flatModel:
            self.rcvID[self.idBorders, :] = np.tile(self.idBorders, (self.flowDir, 1)).T
            self.distRcv[self.idBorders, :] = 0.0
            self.wghtVal[self.idBorders, :] = 0.0

        # Get global nodes that have no receivers
        sum_weight = np.sum(self.wghtVal, axis=1)
        tmpw = np.zeros(self.mpoints, dtype=np.float64)
        tmpw[self.locIDs] = sum_weight
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, tmpw, op=MPI.MAX)
        self.pitPts = tmpw == 0.0

        if MPIrank == 0 and self.verbose:
            print(
                "Flow Direction declaration (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def _distanceCoasts(self, data, k_neighbors=1):
        """
        This function computes for every marine vertices the distance to the closest coastline.

        .. important::

            The calculation takes advantage of the `vtkContourFilter` function from VTK library
            which is performed on the **global** VTK mesh. Once the coastlines have been extracted,
            the distances are obtained by querying a kd-tree (initialised with the coastal nodes) for
            marine vertices contained within each partition.

        :arg data: global elevation numpy array
        :arg k_neighbors: number of nodes to use when querying the kd-tree
        """

        t0 = process_time()

        self.coastDist = np.zeros(self.lpoints)
        pointData = self.vtkMesh.GetPointData()
        array = numpy_support.numpy_to_vtk(data, deep=1)
        array.SetName("elev")
        pointData.AddArray(array)
        self.vtkMesh.SetFieldData(pointData)

        cf = vtk.vtkContourFilter()
        cf.SetInputData(self.vtkMesh)
        cf.SetValue(0, self.sealevel)
        cf.SetInputArrayToProcess(0, 0, 0, 0, "elev")
        cf.GenerateTrianglesOff()
        cf.Update()
        coastXYZ = numpy_support.vtk_to_numpy(cf.GetOutput().GetPoints().GetData())
        tree = spatial.cKDTree(coastXYZ, leafsize=10)
        self.coastDist[self.seaID], indices = tree.query(
            self.lcoords[self.seaID, :], k=k_neighbors
        )

        if self.memclear:
            del array, pointData, cf, tree, indices, coastXYZ
            gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "Construct Distance to Coast (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def _flowSurface(self):
        """
        This function computes the depression less surface.

        .. note::

            In most landscape evolution models, internally draining regions (*e.g.* depressions and pits) are usually filled before the calculation of flow discharge and erosion–deposition rates. This ensures that all flows conveniently reach the coast or the boundary of the simulated domain. In models intended to simulate purely erosional features, such depressions are usually treated as transient features and often ignored.

            However, `gospl` is designed to not only address erosion problems but also to simulate source-to-sink transfer and sedimentary basins formation and evolution in potentially complex tectonic settings. In such cases, depressions may be formed at different periods during runtime and may be filled or remain internally drained (*e.g.* endorheic basins) depending on the volume of sediment transported by upstream catchments.

        The function is **not parallelised** and is performed on the master processor. It calls a fortran subroutine `fillpit` that uses a *priority-flood + ϵ* variant of the algorithm proposed in `Barnes et al. (2014) <https://doi.org/10.1016/j.cageo.2013.04.024>`_ for unstructured meshes.

        The initialisation step consists of pushing all the marine nodes which are neighbours to a deep ocean node (*i.e.* nodes at water depth below 8000 m) onto a priority queue. The priority queue rearranges these nodes so that the ones with the lowest elevations in the queue are always processed first. The algorithm then proceeds incrementally by adding new vertices in the queue and defines filled elevations in the ascending order of their spill elevations + ϵ.

        At the end of the function, the global fill-limited elevations are returned as well as the depressions charcateristics (*i.e.* a unique indice for each depressions, their volumes and the positions of the spill-over global vertices).

        :return: gZ hFill (global and filled elevations numpy arrays)
        """

        # Get global elevations for pit filling...
        t0 = process_time()
        hl = self.hLocal.getArray().copy()
        gZ = np.zeros(self.mpoints, dtype=np.float64)
        gZ[self.locIDs] = hl
        gZ[self.outIDs] = -1.0e8
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gZ, op=MPI.MAX)

        # Perform pit filling on process rank 0
        if MPIrank == 0:
            # fillpit returns:
            # - hFill: filled elevation values
            # - pits: 2D array containing each pit ID and
            #         corresponding overspilling point ID
            hFill, _ = fillpit(self.sealevel - 1.0, gZ, self.waterfill)

        else:
            hFill = np.zeros(self.mpoints, dtype=np.float64)
        hFill = MPI.COMM_WORLD.bcast(hFill, root=0)
        self.hFill = hFill[self.locIDs]

        if MPIrank == 0 and self.verbose:
            print(
                "Get pit filling information (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        if self.memclear:
            del hl
            gc.collect()

        return gZ, hFill

    def matrixFlow(self, flow=True):
        """
        This function defines the flow direction matrices for filled-limited elevations.

        .. note::

            The matrix is built incrementally looping through the number of flow direction paths defined by the user. It proceeds by assembling a local Compressed Sparse Row (**CSR**) matrix to a global PETSc matrix.

            When setting up the flow matrix in PETSc, we preallocate the non-zero entries of the matrix before starting filling in the values. Using PETSc sparse matrix storage scheme has the advantage that matrix-vector multiplication is extremely fast.

        The  matrix coefficients consist of weights (comprised between 0 and 1) calculated based on the number of downslope neighbours and proportional to the slope.

        :arg flow: boolean to compute matrix for either downstream water or sediment transport

        """

        flowMat = self.iMat.copy()
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        nodes = indptr[:-1]

        if flow:
            wght = self.wghtVal
            rcv = self.rcvID
        else:
            wght = self.lWght
            rcv = self.lRcv

        for k in range(0, self.flowDir):

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
            flowMat += tmpMat
            tmpMat.destroy()

        if self.memclear:
            del data, indptr, nodes
            gc.collect()

        # Store flow accumulation matrix
        self.fMat = flowMat.transpose().copy()

        flowMat.destroy()

        return

    def _distributeFlowExcess(self):
        """
        In cases where rivers flow in depressions, they might fill the sink completely and overspill or a portion of it forming a lake. This function computes the lake elevation and the excess of water (if any) able to flow dowstream.

        .. important::

            The excess water is then added to the downstream flow accumulation (`FA`) and used to estimate rivers' erosion.

        """

        t0 = process_time()
        self.hierarchicalSink(True, True, False, None)
        lFA = self.FAL.getArray().copy()
        nFA = np.zeros(self.mpoints, dtype=np.float64)
        nFA[self.locIDs] = lFA
        nFA[self.outIDs] = 0.0
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, nFA, op=MPI.MAX)
        FA = np.zeros(self.mpoints, dtype=np.float64)
        FA[self.pitPts] = nFA[self.pitPts]
        FA[self.sinkID] = 0.0
        if self.flatModel:
            FA[self.glBorders] = 0.0
        self.depo = np.zeros(self.mpoints, dtype=np.float64)
        inFA = np.where(FA > 0)[0]

        if len(inFA) == 0:
            if self.memclear:
                del lFA, FA, nFA, inFA
                gc.collect()
            return

        step = 0
        self.flowDist = True
        while (FA > 0).any():

            # Check if FA are in closed sea
            insea = np.zeros(1, dtype=np.int64)
            if MPIrank == 0:
                if (self.sIns[inFA] > -1).any():
                    insea[0] = 1
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, insea, op=MPI.MAX)

            if insea[0] > 0:
                # Distribute sediments in closed sea
                FA = self.distributeSea(FA)
                # Remove flux that flow into open sea
                FA[self.sinkID] = 0.0
                # In case of 2D model remove flow on mesh borders
                if self.flatModel:
                    FA[self.glBorders] = 0.0

            # At this point excess water are all in continental regions
            if (FA > 0).any():
                FA = self.distributeContinent(FA, None)
                if self.flatModel:
                    FA[self.glBorders] = 0.0
            step += 1
            if step > 20:
                break

        # Store lake water volume for each node
        self.fillFAL.setArray(self.depo[self.locIDs] * self.marea[self.locIDs])
        self.flowDist = False
        if MPIrank == 0 and self.verbose:
            print(
                "Distribute Excess Flow Downstream (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        if self.memclear:
            del FA, nFA, inFA, lFA
            gc.collect()

        return

    def flowAccumulation(self):
        """
        This function is the **main entry point** for flow accumulation computation.

        .. note::

            Flow accumulation (`FA`) calculations are a core component of landscape evolution models as they are often used as proxy to estimate flow discharge, sediment load, river width, bedrock erosion, and sediment deposition. Until recently, conventional `FA` algorithms were serial and limited to small spatial problems.

        `gospl` model computes the flow discharge from `FA` and the net precipitation rate using a **parallel implicit drainage area (IDA) method** proposed by `Richardson et al., 2014 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013WR014326>`_ but adapted to unstructured grids.

        It calls the following *private functions*:

        1. _flowSurface
        2. _buildFlowDirection
        3. _distanceCoasts
        4. _solve_KSP
        5. _distributeFlowExcess

        """

        # self.dm.localToGlobal(self.hLocal, self.hGlobal)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        # Filled-limited surface used to move river flow downstream
        gZ, hFill = self._flowSurface()

        # Build flow direction
        self._buildFlowDirection(hFill[self.locIDs])

        # Define coastal distance for marine points
        if self.vtkMesh is not None:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                self._distanceCoasts(gZ)

        if self.memclear:
            del hFill, gZ
            gc.collect()

        # Build transport direction matrices
        t0 = process_time()
        self.matrixFlow(True)

        # Solve flow accumulation for filled-limited elevation
        if self.tNow == self.rStart:
            self._solve_KSP(False, self.fMat, self.bG, self.FAG)
        else:
            self._solve_KSP(True, self.fMat, self.bG, self.FAG)
        self.dm.globalToLocal(self.FAG, self.FAL)

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Flow Accumulation (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        # Distribute excess flow
        self._distributeFlowExcess()

        return

    def _getErosionRate(self):
        r"""
        This function computes erosion rates in metres per year. This is done on the filled elevation. We use the filled-limited elevation to ensure that erosion is not going to be underestimate by small depressions which are likely to be filled (either by sediments or water) during a single time step.

        The simplest law to simulate fluvial incision is based on the detachment-limited stream power law, in which erosion rate  depends on drainage area :math:`A`, net precipitation :math:`P` and local slope :math:`S` and takes the form:

        .. math::

          E = − \kappa P^d (PA)^m S^n

        :math:`\kappa` is a dimensional coefficient describing the erodibility of the channel bed as a function of rock strength, bed roughness and climate, :math:`d`, :math:`m` and :math:`n` are dimensionless positive constants.

        Default formulation assumes :math:`d = 0`, :math:`m = 0.5` and :math:`n = 1`. The precipitation exponent :math:`d` allows for representation of climate-dependent chemical weathering of river bed across non-uniform rainfall.

        .. important::

            In `gospl`, the coefficients `m` and `n` are fixed and the only variables that the user can tune are the coefficient `d` and the erodibility :math:`\kappa`. E is in m/yr and therefore the erodibility dimension is :math:`(m yr)^{−0.5}`.

        The erosion rate is solved by an implicit time integration method, the matrix system is based on the receiver distributions and is assembled from local Compressed Sparse Row (**CSR**) matrices into a global PETSc matrix. The PETSc *scalable linear equations solvers* (**KSP**) is used with both an iterative method and a preconditioner and erosion rate solution is obtained using PETSc Richardson solver (`richardson`) with block Jacobian preconditioning (`bjacobi`).

        """

        # Upstream-averaged mean annual precipitation rate and the drainage area
        PA = self.FAL.getArray()

        # Incorporate the effect of local mean annual precipitation rate on erodibility
        Kbr = self.K * (self.rainVal ** self.coeffd)
        Kbr *= np.sqrt(PA) * self.dt
        Kbr[self.seaID] = 0.0

        # Initialise identity matrices...
        EbedMat = self.iMat.copy()
        wght = self.wghtVal.copy()

        # Define erosion coefficients
        for k in range(0, self.flowDir):

            indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
            nodes = indptr[:-1]
            # Define erosion limiter to prevent formation of flat
            dh = self.hOldArray - self.hOldArray[self.rcvID[:, k]]
            # Set places which have been filled (e.g. depression/lakes) as no slope regions
            # dh[self.hFill > self.hOldArray] = 0.0
            # dh[dh < 0.0] = 0.0
            limiter = np.divide(dh, dh + 1.0e-3, out=np.zeros_like(dh), where=dh != 0)

            # Bedrock erosion processes SPL computation (maximum bedrock incision)
            data = np.divide(
                Kbr * limiter,
                self.distRcv[:, k],
                out=np.zeros_like(PA),
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
            )
            tmpMat.assemblyEnd()
            EbedMat += tmpMat
            EbedMat -= self._matrix_build_diag(data)
            tmpMat.destroy()

        if self.memclear:
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
        self.dm.globalToLocal(self.Eb, self.EbLocal)
        E = self.EbLocal.getArray().copy()

        E[self.seaID] = 0.0
        ids = np.where(
            (self.hOldArray > self.sealevel + 1.0e-2)
            & (self.hOldArray - E * self.dt < self.sealevel + 1.0e-2)
        )[0]
        E[ids] = (self.hOldArray[ids] - self.sealevel - 1.0e-2) / self.dt
        if self.flatModel:
            E[self.idBorders] = 0.0
        # E[self.hFill > self.hOldArray] = 0.0
        self.EbLocal.setArray(E)
        self.dm.localToGlobal(self.EbLocal, self.Eb)

        if self.memclear:
            del E, PA, Kbr, ids
            gc.collect()

        return

    def riverIncision(self):
        """
        River incision is based on a standard form of the **stream power law** assuming
        **detachment-limited behaviour**.

        It calls the private function `_getErosionRate <https://gospl.readthedocs.io/en/latest/api.html#flow.flowplex.FAMesh._getErosionRate>`_ described above. Once erosion rates have been calculated, the function computes local eroded thicknesses for the considered time step and update local elevation and cumulative erosion, deposition values.

        If multiple lithologies are considered, the stratigraphic record is updated based on eroded thicknesses.

        .. important::

            The approach assumes that the volume of rock eroded using the stream power law accounts for both the solid and void phase.

        """

        t0 = process_time()

        # Local & global vectors/arrays
        self.Eb.set(0.0)
        self.hGlobal.copy(result=self.hOld)
        self.dm.globalToLocal(self.hOld, self.hOldLocal)
        self.hOldArray = self.hOldLocal.getArray().copy()

        self._getErosionRate()

        # Get eroded thicknesses
        Eb = self.Eb.getArray().copy()
        self.tmp.setArray(-Eb * self.dt)
        self.cumED.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)

        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        # Update stratigraphic layers
        if self.stratNb > 0:
            self.erodeStrat()

        if MPIrank == 0 and self.verbose:
            print(
                "Get Erosion Thicknesses (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        if self.memclear:
            del Eb
            gc.collect()

        return
