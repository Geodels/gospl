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
    from gospl._fortran import jacobiancoeff
    from gospl._fortran import fctcoeff
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

        self.coastDist = None

        self.zMat = self._matrix_build_diag(np.zeros(self.lpoints))
        self.mat = self.dm.createMatrix()
        self.mat.setOption(self.mat.Option.NEW_NONZERO_LOCATIONS, True)

        # Clinoforms slopes for each sediment type
        # self.slps = [12.0, 8.0, 6.0, 1.0]

        self.dh = self.hGlobal.duplicate()
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
        # rcv, _, wght = mfdreceivers(8, 1.0, self.inIDs, self.smthH, -1.0e5)
        rcv, _, wght = mfdreceivers(8, 1.0e-2, self.inIDs, self.oceanFill, -1.0e5)
        sum_wght = np.sum(wght, axis=1)
        ids = (self.pitIDs > -1) & (self.flatOcean > -1) & (sum_wght == 0.0)
        ids = ids.nonzero()[0]
        rcv[ids, :] = np.tile(ids, (8, 1)).T
        rcv[ids, 0] = self.flatOcean[ids]
        wght[ids, :] = 0.0
        wght[ids, 0] = 1.0

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

    def _evalFunction(self, ts, t, x, xdot, f):
        """
        The nonlinear system at each time step is solved iteratively using PETSc time stepping and SNES method and is based on a Newton/GMRES method. . Here we define the function for the nonlinear solve.

        """

        self.dm.globalToLocal(x, self.hl)
        with self.hl as hl, self.hLocal as zb, xdot as hdot:
            dh = hl - zb
            dh[dh < 0] = 0.0
            Cd = np.multiply(self.Cd, dh ** self.cexp)
            nlvec = fctcoeff(hl, Cd)
            # Set borders nodes
            if self.flatModel:
                nlvec[self.idBorders] = 0.0
            f.setArray(hdot + nlvec[self.glIDs])

    def _evalJacobian(self, ts, t, x, xdot, a, J, P):
        """
        The nonlinear system at each time step is solved iteratively using PETSc time stepping and SNES method and is based on a Newton/GMRES method. Here we define the Jacobian for the nonlinear solve.

        """

        self.dm.globalToLocal(x, self.hl)

        with self.hl as hl, self.hLocal as zb:
            dh = hl - zb
            dh[dh < 0] = 0.0
            Cd = np.multiply(self.Cd, dh ** self.cexp)

            # Coefficient derivatives
            Cp = self.cexp * np.multiply(self.Cd, dh ** (self.cexp - 1.0))

            nlC = jacobiancoeff(hl, Cd, Cp)
            if self.flatModel:
                nlC[self.idBorders, :] = 0.0

            P.zeroEntries()
            for row in range(self.lpoints):
                P.setValuesLocal(row, row, a + nlC[row, 0])
                cols = self.FVmesh_ngbID[row, :]
                P.setValuesLocal(row, cols, nlC[row, 1:])
            P.assemble()

            if J != P:
                J.assemble()

        return petsc4py.PETSc.Mat.Structure.SAME_NONZERO_PATTERN

    def _diffuseOcean(self, dh, stype):
        r"""
        For sediment reaching the marine realm, this function computes the related marine deposition diffusion. The approach is based on a nonlinear diffusion.

        .. math::
          \frac{\partial h}{\partial t}= \nabla \cdot \left( C_d(h) \nabla h \right)

        It calls the following *private functions*:

        - _evalFunction
        - _evalJacobian

        .. note::

            PETSc SNES and time stepping TS approaches are used to solve the nonlinear equation above over the considered time step.

        The nonlinear diffusion equation is ran for the coarser sediment first and for the finest ones afterwards. This mimicks the standard behaviour observed in stratigraphic architectures where the fine fraction are generally transported over longer distances.

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
        self.Cd = np.full(self.lpoints, sedK, dtype=np.float64)
        self.hl.setArray(dh)
        self.hl.axpy(1.0, self.hLocal)
        self.dm.localToGlobal(self.hl, self.h)

        # Time stepping definition
        ts = petsc4py.PETSc.TS().create(comm=petsc4py.PETSc.COMM_WORLD)
        # ARKIMEX: implicit nonlinearl time stepping
        ts.setType("arkimex")
        ts.setIFunction(self._evalFunction, self.tmp1)
        ts.setIJacobian(self._evalJacobian, self.mat)

        ts.setTime(0.0)
        ts.setTimeStep(self.dt / 1000.0)
        ts.setMaxTime(self.dt)
        ts.setMaxSteps(50)
        ts.setExactFinalTime(petsc4py.PETSc.TS.ExactFinalTime.MATCHSTEP)
        # Allow an unlimited number of failures (step will be rejected and retried)
        ts.setMaxSNESFailures(-1)
        ts.setTolerances(rtol=1.0e-2)

        # SNES nonlinear solver definition
        snes = ts.getSNES()
        # Newton linear search
        snes.setType("newtonls")
        # Stop nonlinear solve after 10 iterations (TS will retry with shorter step)
        snes.setTolerances(rtol=1.0e-2, max_it=100)

        # KSP linear solver definition
        ksp = snes.getKSP()
        ksp.setType("fgmres")
        # Preconditioner for linear solution
        pc = ksp.getPC()
        pc.setType("asm")
        pc.setASMOverlap(3)
        ksp.setTolerances(rtol=1.0e-2, max_it=100)
        ksp.setFromOptions()
        snes.setFromOptions()
        ts.setFromOptions()

        # Solve nonlinear equation
        ts.solve(self.h)

        if MPIrank == 0 and self.verbose:
            print(
                "Nonlinear diffusion solution (%0.02f seconds)" % (process_time() - t0),
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
                )
            )
        ts.destroy()

        # Get diffused sediment thicknesses
        self.dh.waxpy(-1.0, self.hGlobal, self.h)
        self.dm.globalToLocal(self.dh, self.tmpL)
        ndepo = self.tmpL.getArray().copy()

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
        if self.clinSlp > 0.0:
            clinoH = self.sealevel - 1.0e-3 - self.coastDist * self.clinSlp
        else:
            clinoH = np.full(self.lpoints, self.sealevel - 1.0e-3, dtype=np.float64)
        # Update the marine maximal depositional thicknesses
        clinoH = np.minimum(clinoH, self.smthH + self.offset)
        clinoH[hl >= self.sealevel] = hl[hl >= self.sealevel]

        # Get the maximum marine deposition heights and volume
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
        self.tmpL.setArray(hl)
        self.dm.localToGlobal(self.tmpL, self.tmp1)
        smthH1 = self._hillSlope(smooth=2)
        self.tmpL.setArray(marDep + hl)
        self.dm.localToGlobal(self.tmpL, self.tmp1)
        smthH2 = self._hillSlope(smooth=2)
        marDep = smthH2 - smthH1
        marDep[marDep < 0.0] = 0.0
        marDep[hl > self.sealevel] = 0.0
        ids = marDep > clinoH - hl
        marDep[ids] = clinoH[ids] - hl[ids]
        marDep[marDep < 0.0] = 0.0
        if self.marineNl:
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

        # Get the smooth regional elevation
        self.smthH = self._hillSlope(smooth=1)

        # Downstream direction matrix for ocean distribution
        self._matOcean()

        # Define coastal distance for marine points
        self.dm.globalToLocal(self.hGlobal, self.hLocal)
        hl = self.hLocal.getArray().copy()

        if self.clinSlp > 0.0:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                self._distanceCoasts(hl)

        # Set all nodes below sea-level as sinks
        self.sinkIDs = self.lFill < self.sealevel

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
