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
    from gospl._fortran import mfdrcvrs

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()


class SEAMesh(object):
    """
    This class encapsulates all the functions related to sediment transport and deposition in the marine environment for **river delivered sediments**. `gospl` has the ability to track three types of clastic sediment size and one type of carbonate (still under development).

    .. note::
        All of these functions are ran in parallel using the underlying PETSc library.

    For an overview of solution to nonlinear ODE and PDE problems, one might found the online book from `Langtangen (2016) <http://hplgit.github.io/num-methods-for-PDEs/doc/pub/nonlin/pdf/nonlin-4print-A4-2up.pdf>`_ relevant.

    """

    def __init__(self, *args, **kwargs):
        """
        The initialisation of `SEAMesh` class consists in the declaration of several PETSc vectors.
        """

        self.coastDist = None

        self.zMat = self._matrix_build_diag(np.zeros(self.lpoints))

        self.dh = self.hGlobal.duplicate()

        if self.excessIn:
            self.ePit = np.zeros(self.lpoints)

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
        This function computes for every marine vertices the distance to the closest coastline.

        .. important::

            The calculation takes advantage of the `vtkContourFilter` function from VTK library
            which is performed on the **global** VTK mesh. Once the coastlines have been extracted,
            the distances are obtained by querying a kd-tree (initialised with the coastal nodes) for
            marine vertices contained within each partition.

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

    def _waterFluxes(self, sedflux):
        """
        Based on the incoming marine volumes of sediment and maximum clinoforms slope we distribute
        sediments downslope.

        """
        
        # Define multiple flow directions under water
        self.donRcvs, self.distRcv, self.wghtVal = mfdrcvrs(
            self.flowDir, 0.01, self.oceanFill, self.sealevel
        )
        
        self.rcvID = self.donRcvs.copy()
        self.rcvID[self.ghostIDs,:] = -1
        self.distRcv[self.ghostIDs,:] = 0
        self.wghtVal[self.ghostIDs,:] = 0

        # Set borders nodes
        if self.flatModel:
            self.rcvID[self.idBorders, :] = np.tile(self.idBorders, (8, 1)).T
            self.distRcv[self.idBorders, :] = 0.0
            self.wghtVal[self.idBorders, :] = 0.0

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
        fDep = np.divide(self.fDepm*self.larea, PA, out=np.zeros_like(PA), where=PA > 1.e-6)
        if self.dmthd == 1:
            fDep[fDep>0.99] = 0.99
            self.matrixFlow(8, 1.-fDep)
        else:
            dMat = self._matrix_build_diag(fDep)
            dMat += self.fMat

        # Implicit sediment fluxes combining upstream flux and deposition
        self.tmpL.setArray(sedflux)
        self.dm.localToGlobal(self.tmpL, self.tmp1)
        if self.dmthd == 1:
            self._solve_KSP(False, self.fDepMat, self.tmp1, self.tmp)
            self.fDepMat.destroy()
        else:
            self._solve_KSP(True, dMat, self.tmp1, self.tmp)
            dMat.destroy()

        # Destroy temporary arrays

        if self.memclear:
            del PA, FAL 
            gc.collect()

        # Extract local sediment deposition thickness
        self.dm.globalToLocal(self.tmp, self.tmpL)
        if self.dmthd == 1:
            scale = np.divide(fDep, 1.0-fDep, out=np.zeros_like(fDep), where=fDep != 0)
            sedDep = self.tmpL.getArray()*scale
        else:
            sedDep = self.tmpL.getArray()*fDep
        if self.flatModel:
            sedDep[self.idBorders] = 0.0
        self.EbLocal.setArray(-sedDep/self.larea)
        self.dm.localToGlobal(self.EbLocal, self.Eb)

        # Limiting deposition rates based on upstream conditions
        self._upstreamDeposition(8)

        # Get deposition thicknesses
        Eb = self.EbLocal.getArray().copy()
        Eb[Eb>0] = 0.0

        # Define coastal distance for marine points
        self.dm.globalToLocal(self.hGlobal, self.hLocal)
        hl = self.hLocal.getArray().copy()
        if self.clinSlp > 0.0:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                self._distanceCoasts(hl)

        # From the distance to coastline define the upper limit of the shelf to ensure a maximum slope angle
        if self.clinSlp > 0.0:
            clinoH = self.sealevel - 1.0e-3 - self.coastDist * self.clinSlp
        else:
            clinoH = np.full(self.lpoints, self.sealevel - 1.0e-3, dtype=np.float64)
        clinoH[hl >= self.sealevel] = hl[hl >= self.sealevel]
        
        # Update the marine maximal depositional thicknesses
        updateH = hl - Eb * self.dt
        updateH[updateH >= clinoH] = clinoH[updateH >= clinoH]
        self.tmpL.setArray(updateH - hl)

        # Smoothing marine deposition
        if self.smthD > 0:
            self.dm.localToGlobal(self.tmpL, self.tmp1)
            smthH = self._hillSlope(smooth=1) + hl
            smthH[smthH >= clinoH] = clinoH[smthH >= clinoH]
            self.tmpL.setArray(smthH - hl)
            self.dm.localToGlobal(self.tmpL, self.tmp)
        else:
            self.dm.localToGlobal(self.tmpL, self.tmp)

        if self.memclear:
            del updateH, Eb, hl, clinoH, fDep
            gc.collect()

        return 

    def seaChange(self):
        """
        This function is the main entry point to perform marine river-induced deposition. It calls the private functions:

        - _distanceCoasts
        - _waterFluxes

        """

        t0 = process_time()
        
        # Set all nodes below sea-level as sinks
        self.sinkIDs = self.lFill <= self.sealevel
        
        # Get the volumetric marine sediment rate (m3/yr) to distribute during the time step and convert it in volume (m3)
        self.vSedLocal.copy(result=self.QsL)
        sedFlux = self.QsL.getArray().copy() 
        sedFlux[np.invert(self.sinkIDs)] = 0.0
        flxStp = sedFlux/self.diffNb
        for k in range(self.diffNb):
            # Compute marine directions and fluxes
            self._waterFluxes(flxStp)

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
