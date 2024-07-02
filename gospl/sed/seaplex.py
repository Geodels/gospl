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
from vtk.util import numpy_support  # type: ignore

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import donorslist
    from gospl._fortran import donorsmax
    from gospl._fortran import mfdrcvrs

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()


class SEAMesh(object):
    """
    This class encapsulates all the functions related to sediment transport and deposition in the marine environment for **river delivered sediments**.

    .. note::
        All of these functions are run in parallel using the underlying PETSc library.

    For an overview of solution to nonlinear ODE and PDE problems, one might found the online book from `Langtangen (2016) <http://hplgit.github.io/num-methods-for-PDEs/doc/pub/nonlin/pdf/nonlin-4print-A4-2up.pdf>`_ relevant.

    """

    def __init__(self, *args, **kwargs):
        """
        The initialisation of `SEAMesh` class consists in the declaration of several PETSc vectors.
        """

        self.coastDist = None

        self.zMat = self._matrix_build_diag(np.zeros(self.lpoints))
        self.dh = self.hGlobal.duplicate()

        return

    def _globalCoastsTree(self, coastXYZ, k_neighbors=1):
        """
        This function takes all local coastline points and computes locally the distance of all marine points to the coastline.

        :arg coastXYZ: local coastline coordinates
        :arg k_neighbors: number of nodes to use when querying the kd-tree
        """

        label_offset = np.zeros(MPIsize + 1, dtype=int)
        label_offset[MPIrank + 1] = len(coastXYZ)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, label_offset, op=MPI.MAX)
        offset = np.cumsum(label_offset)[:-1]
        totpts = np.sum(label_offset)

        coastpts = np.zeros(((totpts) * 3,), dtype=np.float64)
        MPI.COMM_WORLD.Allgatherv(
            [coastXYZ, MPI.DOUBLE],
            [coastpts, label_offset[1:] * 3, offset * 3, MPI.DOUBLE],
        )

        # Get coastlines points globally
        tree = spatial.cKDTree(coastpts.reshape((totpts, 3)), leafsize=10)
        self.coastDist[self.seaID], _ = tree.query(
            self.lcoords[self.seaID, :], k=k_neighbors
        )

        if self.memclear:
            del tree, coastpts, offset, totpts, label_offset
            gc.collect()

        return

    def _distanceCoasts(self, data, k_neighbors=1):
        """
        This function computes for every marine vertices the distance to the closest coastline. It calls the private functions:

        - _globalCoastsTree

        .. important::

            The calculation takes advantage of the `vtkContourFilter` function from VTK library which is performed on the **global** VTK mesh. Once the coastlines have been extracted, the distances are obtained by querying a kd-tree (initialised with the coastal nodes) for marine vertices contained within each partition.

        :arg data: local elevation numpy array
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
        if cf.GetOutput().GetPoints() is not None:
            coastXYZ = numpy_support.vtk_to_numpy(cf.GetOutput().GetPoints().GetData())
        else:
            coastXYZ = np.zeros((0, 3), dtype=np.float64)

        # Get coastlines points globally
        self._globalCoastsTree(coastXYZ, k_neighbors)

        if self.memclear:
            del array, pointData, cf, coastXYZ
            gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "Construct Distance to Coast (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def _upstreamDeposition(self, flowdir):
        """
        In cases where too much deposition occurs in depressions and/or flat regions, this function enforces a maximum deposition rate based on upstream elevations computed from the implicit erosion/deposition equation.

        :arg flowdir: number of flow direction
        """
        niter = 0
        tolerance = 1.e-3
        equal = False
        h = self.hLocal.getArray().copy()

        # Get donors list based on rcvs
        donors = donorslist(flowdir, self.inIDs, self.donRcvs)
        topIDs = np.where(np.max(donors, axis=1) == -1)[0]
        while not equal and niter < 1000:

            dh = self.tmpL.getArray()  # -self.EbLocal.getArray().copy() * self.dt
            elev = (dh + h).copy()
            elev[topIDs] = h[topIDs]
            maxh = donorsmax(elev, donors)
            maxh[maxh == -1.e8] = elev[maxh == -1.e8]
            self.tmpL.setArray(maxh)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.dm.globalToLocal(self.tmp, self.tmpL)
            maxh = self.tmpL.getArray().copy() - 1.e-6
            maxh[maxh > elev] = elev[maxh > elev]

            self.tmpL.setArray(maxh - elev)
            self.dm.localToGlobal(self.tmpL, self.dh)
            self.dh.abs()
            maxdh = self.dh.max()[1]

            equal = maxdh < tolerance
            self.tmpL.setArray(maxh - h)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.dm.globalToLocal(self.tmp, self.tmpL)
            niter += 1

        if self.memclear:
            del donors, h, elev, maxh, Eb
            gc.collect()

        return

    def _marineFluxes(self, sedflux):
        """
        Based on the incoming marine volumes of sediment and maximum clinoforms slope we distribute sediments downslope.

        :arg sedflux: volumetric marine sediment rate
        """

        # Define multiple flow directions under water
        hl = self.hLocal.getArray().copy()
        self.donRcvs, self.distRcv, self.wghtVal = mfdrcvrs(
            self.flowDir, 0.01, hl, self.sealevel
        )
        self.rcvID = self.donRcvs.copy()
        self.rcvID[self.ghostIDs, :] = -1
        self.distRcv[self.ghostIDs, :] = 0
        self.wghtVal[self.ghostIDs, :] = 0

        # Define the flow direction matrix
        self.matrixFlow(8)
        FAL = self.fillFAL.getArray().copy()
        FAL[np.invert(self.sinkIDs)] = 0.0
        self.tmpL.setArray(FAL)
        self.dm.localToGlobal(self.tmpL, self.tmp)

        # Compute fluxes
        self._solve_KSP(True, self.fMat, self.tmp, self.tmp1)
        self.dm.globalToLocal(self.tmp1, self.tmpL)

        # Dimensionless depositional coefficient
        PA = self.tmpL.getArray().copy()
        fDep = np.divide(self.fDepm * self.larea, PA, out=np.zeros_like(PA), where=PA > 1.e-6)
        fDep[fDep > 0.99] = 0.99
        self.fDep = fDep
        ids = np.where(hl > self.clinoH)[0]
        self.fDep[ids] = 0.

        # Set the RHS vector to the incoming marine sediment flux
        sedflux[self.idBorders] = 0.
        self.tmpL.setArray(sedflux)
        self.dm.localToGlobal(self.tmpL, self.tmp1)

        tolerance = 1.e-2
        niter = 0
        sumExcess = tolerance + 1.0
        upH = hl.copy()
        self.fDep[upH >= self.clinoH] = 0.

        while sumExcess > tolerance and niter < 100:
            self.matrixFlow(8, 1. - self.fDep)
            self._solve_KSP(True, self.fDepMat, self.tmp1, self.tmp)
            self.fDepMat.destroy()

            # Get the corresponding sedimentation flux
            self.dm.globalToLocal(self.tmp, self.tmpL)
            qs = self.tmpL.getArray().copy()
            self.fMat.transpose().mult(self.tmp, self.tmp1)
            self.dm.globalToLocal(self.tmp1, self.tmpL)
            sedDep = (qs - self.tmpL.getArray()) * self.fDep
            sedDep[sedDep < 0] = 0.
            self.tmpL.setArray(sedDep * self.dt / self.larea)
            depH = self.tmpL.getArray().copy()
            depH[depH < 0] = 0.

            # Get over-deposition
            currH = depH + upH
            ids = np.where(currH > self.clinoH)[0]
            currH[ids] = self.clinoH[ids]
            self.fDep[ids] = 0.
            excess = (depH + upH - currH) * self.larea / self.dt
            excess[excess < 0.] = 0.
            excess[hl > self.sealevel] = 0
            self.tmpL.setArray(excess)
            upH = currH.copy()
            self.dm.localToGlobal(self.tmpL, self.tmp1)
            niter += 1
            sumExcess = self.tmp1.sum()
            if MPIrank == 0 and self.verbose:
                print(
                    "  --- Marine excess (sum) %0.01f m"
                    % (sumExcess),
                    flush=True
                )

        # upH[lsink] = hl[lsink]
        self.tmpL.setArray(upH - hl)

        # Limiting deposition rates based on upstream conditions
        self._upstreamDeposition(8)

        # Smoothing marine deposition
        if self.smthD > 0:
            self.dm.localToGlobal(self.tmpL, self.tmp1)
            smthH = self._hillSlope(smooth=1) + hl
            smthH[smthH >= self.clinoH] = self.clinoH[smthH >= self.clinoH]
            self.tmpL.setArray(smthH - hl)
            self.dm.localToGlobal(self.tmpL, self.tmp)
        else:
            self.dm.localToGlobal(self.tmpL, self.tmp)

        return

    def seaChange(self):
        """
        This function is the main entry point to perform marine river-induced deposition. It calls the private functions:

        - _distanceCoasts
        - _marineFluxes

        """

        t0 = process_time()

        # Set all nodes below sea-level as sinks
        self.sinkIDs = self.lFill <= self.sealevel

        # Define coastal distance for marine points
        if self.clinSlp > 0.0:
            self.dm.globalToLocal(self.hGlobal, self.hLocal)
            hl = self.hLocal.getArray().copy()
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                self._distanceCoasts(hl)
            # From the distance to coastline define the upper limit of the shelf to ensure a maximum slope angle
            self.clinoH = self.sealevel - 1.0e-3 - self.coastDist * self.clinSlp
        else:
            self.clinoH = np.full(self.lpoints, self.sealevel - 1.0e-3, dtype=np.float64)
        self.clinoH[hl >= self.sealevel] = hl[hl >= self.sealevel]
        self.maxDepQs = (self.clinoH - hl) * self.larea / self.dt

        # Get the volumetric marine sediment rate (m3/yr) to distribute during the time step and convert it in volume (m3)
        self.vSedLocal.copy(result=self.QsL)
        sedFlux = self.QsL.getArray().copy()
        sedFlux[np.invert(self.sinkIDs)] = 0.0
        sedFlux[sedFlux < 0] = 0.0
        self._marineFluxes(sedFlux)

        # Update cumulative erosion and deposition as well as elevation
        self.cumED.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        # Update stratigraphic layer parameters
        if self.stratNb > 0:
            self.deposeStrat()

        if MPIrank == 0 and self.verbose:
            print(
                "Distribute River Sediments in the Ocean (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        return
