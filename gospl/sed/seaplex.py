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
    from gospl._fortran import sethillslopecoeff
    from gospl._fortran import setjacobiancoeff
    from gospl._fortran import mfdreceivers
    from gospl._fortran import distocean

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

        self.nleps = 1.0e-4
        self.max_iter = 500
        self.coastDist = None

        self.zMat = self._matrix_build_diag(np.zeros(self.lpoints))

        # Clinoforms slopes for each sediment type
        self.slps = [6.0, 5.0, 3.0, 1.0]

        self.dhl = self.hLocal.duplicate()
        self.dh = self.hGlobal.duplicate()
        self.h0 = self.hGlobal.duplicate()
        self.h = self.hGlobal.duplicate()
        self.hl = self.hLocal.duplicate()
        self.zb = self.hGlobal.duplicate()

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

    def _matOcean(self):
        """
        This function builds from neighbouring slopes the downstream directions in the marine environment. It calls a fortran subroutine that locally computes for each vertice:

        - the indices of receivers (downstream) nodes depending on the desired number of flow directions (SFD to MFD).
        - the distances to the receivers based on mesh resolution.
        - the associated weights calculated based on the number of receivers and proportional to the slope.

        From these downstream directions, a local ocean downstream matrix is computed.

        """

        # Define multiple flow directions for filled + eps elevations
        rcv, _, wght = mfdreceivers(8, 1.0e-2, self.inIDs, self.oceanFill, -1.0e5)

        # Set borders nodes
        if self.flatModel:
            rcv[self.idBorders, :] = np.tile(self.idBorders, (8, 1)).T
            wght[self.idBorders, :] = 0.0

        # Define downstream matrix based on filled + eps elevations
        dMat = self.zMat.copy()
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        nodes = indptr[:-1]

        for k in range(0, 8):

            # Flow direction matrix for a specific direction
            tmpMat = self._matrix_build()
            data = wght[:, k].copy()
            data[rcv[:, k].astype(petsc4py.PETSc.IntType) == nodes] = 0.0
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                rcv[:, k].astype(petsc4py.PETSc.IntType),
                data,
            )
            tmpMat.assemblyEnd()

            # Add the weights from each direction
            dMat += tmpMat
            tmpMat.destroy()

        if self.memclear:
            del data, indptr, nodes
            gc.collect()

        # Store flow direction matrix
        self.dMat = dMat.transpose().copy()

        dMat.destroy()

        return

    def _nlCoefficents(self):
        r"""
        Nonlinear diffusion coefficients based on each iteration freshly deposited thicknesses. In addition diffusion coefficient derivatives with respect to the independent variable `h` in the function below is computed and used in the Jacobian calculation.

        .. math::
          C_d = \kappa_{Ds} (h-z_b)^2

        in which :math:`\kappa_{Ds}` is the diffusion coefficient for river sediment transport (set with `sedK` in the YAML input file) and :math:`z_b` the basement elevation at the start of the iteration.

        """

        self.dh.waxpy(-1.0, self.hGlobal, self.h)
        self.dm.globalToLocal(self.dh, self.dhl)
        dh = self.dhl.getArray().copy()

        # Non linear coefficients
        self.C = np.multiply(self.Cd, dh ** 2)

        # Derivative from h
        self.Cp = 2.0 * dh * self.Cd.copy()

        return

    def _RHSvector(self):
        """
        The nonlinear system at each time step is solved iteratively using a Newton/GMRES method. In the solution process, the right-hand side vector has to be evaluated.

        """

        self._nlCoefficents()
        rhsCoeffs = sethillslopecoeff(self.lpoints, self.C)
        if self.flatModel:
            rhsCoeffs[self.idBorders, 1:] = 0.0
            rhsCoeffs[self.idBorders, 0] = 1.0

        rhs = self._matrix_build_diag(rhsCoeffs[:, 0])
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)

        for k in range(0, self.maxnb):
            tmpMat = self._matrix_build()
            indices = self.FVmesh_ngbID[:, k].copy()
            data = rhsCoeffs[:, k + 1]
            ids = np.nonzero(data == 0.0)
            indices[ids] = ids
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                indices.astype(petsc4py.PETSc.IntType),
                data,
            )
            tmpMat.assemblyEnd()
            rhs += tmpMat
            tmpMat.destroy()

        # Compute RHS vector
        rhs.mult(self.h, self.tmp1)
        self.tmp.waxpy(-1.0, self.tmp1, self.h0)
        rhs.destroy()

        return

    def _Jacobian(self):
        """
        The nonlinear system at each time step is solved iteratively using a Newton/GMRES method. In the solution process, the Jacobian matrix–vector products have to be evaluated.

        """

        h = self.hl.getArray().copy()
        jCoeffs = setjacobiancoeff(self.lpoints, h, self.C, self.Cp)
        if self.flatModel:
            jCoeffs[self.idBorders, 1:] = 0.0
            jCoeffs[self.idBorders, 0] = 1.0

        J = self._matrix_build_diag(jCoeffs[:, 0])
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)

        for k in range(0, self.maxnb):
            tmpMat = self._matrix_build()
            indices = self.FVmesh_ngbID[:, k].copy()
            data = jCoeffs[:, k + 1]
            ids = np.nonzero(data == 0.0)
            indices[ids] = ids
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                indices.astype(petsc4py.PETSc.IntType),
                data,
            )
            tmpMat.assemblyEnd()
            J += tmpMat
            tmpMat.destroy()

        ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
        ksp.setInitialGuessNonzero(True)
        ksp.setOperators(J, J)
        ksp.setType("gmres")
        pc = ksp.getPC()
        pc.setType("asm")
        ksp.setTolerances(rtol=self.rtol)
        ksp.solve(self.tmp, self.dh)

        ksp.destroy()
        J.destroy()

        return

    def _diffuseOcean(self, dh, stype):
        r"""
        For sediment reaching the marine realm, this function computes the related marine deposition diffusion. The approach is based on a nonlinear diffusion.

        .. math::
          \frac{\partial h}{\partial t}= \nabla \cdot \left( C_d(h) \nabla h \right)

        .. note::

            PETSc **KSP** component provides Krylov subspace iterative method and several preconditioner. Here, we use Newton's method to solve the nonlinear equation using generalized minimal residual solver (`gmres`) with the Additive Schwarz method preconditioning (`asm`).

        It calls the following *private functions*:

        - _RHSvector
        - _Jacobian

        .. todo::

            SNES and time stepping approach might provide a better convergence rate than the one implemented here.

        The nondiffusion equation is ran for the coarser sediment first and for the finest ones afterwards. This mimicks the standard behaviour observed in stratigraphic architectures where the fine fraction are generally transported over longer distances.

        :arg dh: numpy array of incoming marine depositional thicknesses
        :arg stype: sediment type (integer)

        :return: ndepo (updated deposition numpy arrays)
        """

        t0 = process_time()

        # Get diffusion coefficients based on sediment type
        sedK = self.sedimentK
        if sedK > 0.0:
            if stype == 1:
                sedK = self.sedimentKf
            if stype == 3:
                sedK = self.sedimentKw

        # Matrix coefficients
        self.Cd = np.full(self.lpoints, sedK, dtype=np.float64) * self.dt

        self.dhl.setArray(dh)
        self.dm.localToGlobal(self.dhl, self.dh)
        self.dm.globalToLocal(self.dh, self.dhl)
        self.h0.waxpy(1.0, self.dh, self.hGlobal)
        self.h0.copy(result=self.h)
        self.dm.globalToLocal(self.h, self.hl)

        # Non linear diffusion
        residual = 1.0
        step = 0
        while residual > self.nleps and step < self.max_iter:
            self._RHSvector()
            self._Jacobian()
            self.h.axpy(1.0, self.dh)
            h = self.h.getArray().copy()
            zb = self.hGlobal.getArray().copy()
            h[zb > h] = zb[zb > h]
            self.h.setArray(h)
            self.dm.globalToLocal(self.h, self.hl)
            residual = np.abs(self.dh.max()[1])
            step += 1

        if MPIrank == 0 and self.verbose:
            print(
                "Non-linear convergence steps %d with residual %0.05f"
                % (step, residual),
                flush=True,
            )
            print(
                "Diffuse Marine Sediments (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        self.dh.waxpy(-1.0, self.hGlobal, self.h)
        self.dm.globalToLocal(self.dh, self.dhl)
        ndepo = self.dhl.getArray().copy()
        if self.flatModel:
            ndepo[self.idBorders] = 0.0
        ndepo[ndepo < 0.0] = 0.0

        return ndepo

    def _distOcean(self):
        """
        Based on the incoming marine volumes of sediment and maximum clinoforms slope we distribute
        sediments downslope.

        """

        vol = self.marVol.copy()
        self.tmpL.setArray(self.sinkVol)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        vdep = np.zeros(self.lpoints, dtype=np.float64)

        step = 0
        while self.tmp.sum() > 1.0e-6:

            # Move to downstream nodes
            self.dMat.mult(self.tmp, self.tmp1)
            self.dm.globalToLocal(self.tmp1, self.tmpL)

            # In case there is to much sediment coming in
            self.sinkVol = self.tmpL.getArray().copy()
            excess = self.sinkVol >= vol
            self.sinkVol[excess] -= vol[excess]
            vdep[excess] = self.marVol[excess]
            vol[excess] = 0.0

            # In case there is some room to deposit sediments
            noexcess = np.invert(excess)
            vol[noexcess] -= self.sinkVol[noexcess]
            vdep[noexcess] += self.sinkVol[noexcess]
            vdep[self.idBorders] = 0.0
            self.sinkVol[noexcess] = 0.0
            self.sinkVol[self.idBorders] = 0.0
            self.tmpL.setArray(self.sinkVol)
            self.dm.localToGlobal(self.tmpL, self.tmp)

            step += 1

        return vdep / self.larea

    def _getSeaVol(self, hl, stype):
        """
        Pick the relevant PETSc array for the specified sediment type and distribute incoming
        river delivered sediment volume in the marine environment.

        The function relies on two private functions from the class:

        - _distOcean
        - _diffuseOcean

        """

        # From the distance to coastline define the upper limit of the shelf to ensure a maximum slope angle
        clinoH = self.sealevel - 1.0e-4 - self.coastDist * self.slps[stype] * 1.0e-5

        # Get the maximum marine deposition heights and volume
        # self.marVol = depo * self.larea
        self.marVol = (clinoH - hl) * self.larea
        self.marVol[self.marVol < 0.0] = 0.0

        # Get the volumetric marine sediment rate (m3/yr) to distribute during the time step and convert it in volume (m3)
        if stype == 0:
            self.vSedLocal.copy(result=self.QsL)
        elif stype == 1:
            self.vSedfLocal.copy(result=self.QsL)
        elif stype == 2:
            self.vSedcLocal.copy(result=self.QsL)
        elif stype == 3:
            self.vSedwLocal.copy(result=self.QsL)

        self.sinkVol = self.QsL.getArray().copy() * self.dt
        self.sinkVol[np.invert(self.sinkIDs)] = 0.0

        # Distribute sediment downstream in the marine environment
        marDep = self._distOcean()
        marDep = self._diffuseOcean(marDep, stype)

        # Update cumulative erosion and deposition as well as elevation
        self.tmpL.setArray(marDep)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.cumED.axpy(1.0, self.tmp)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        # Update stratigraphic layer parameters
        if self.stratNb > 0:
            self.deposeStrat(stype)
            self.elevStrat()

        return

    def seaChange(self):
        """
        This function is the main entry point to perform marine river-induced deposition. It calls the private functions:

        -  _matOcean
        - _distanceCoasts
        - _getSeaVol

        """

        t0 = process_time()

        # Downstream direction matrix for ocean distribution
        self._matOcean()

        # Define coastal distance for marine points
        self.dm.globalToLocal(self.hGlobal, self.hLocal)
        hl = self.hLocal.getArray().copy()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self._distanceCoasts(hl)

        # Set all nodes below sea-level as sinks
        self.sinkIDs = self.epsFill < self.sealevel

        # Get sediment volume to distribute
        if self.stratNb > 0:
            self._getSeaVol(hl, stype=0)
            if self.stratF is not None:
                self._getSeaVol(hl, stype=1)
            if self.carbOn:
                self._getSeaVol(hl, stype=2)
            if self.stratW is not None:
                self._getSeaVol(hl, stype=3)
        else:
            self._getSeaVol(hl, stype=0)

        if MPIrank == 0 and self.verbose:
            print(
                "Distribute River Sediments in the Ocean (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        return
