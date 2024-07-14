import os
import gc
import sys
import vtk
import warnings
import petsc4py
import numpy as np
from scipy import spatial

from mpi4py import MPI
from time import process_time
from vtk.util import numpy_support  # type: ignore

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import donorslist
    from gospl._fortran import donorsmax
    from gospl._fortran import mfdrcvrs
    from gospl._fortran import jacobiancoeff
    from gospl._fortran import fctcoeff

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

        self.mat = self.dm.createMatrix()
        self.mat.setOption(self.mat.Option.NEW_NONZERO_LOCATIONS, True)

        self.zMat = self._matrix_build_diag(np.zeros(self.lpoints))
        self.dh = self.hGlobal.duplicate()
        self.h = self.hGlobal.duplicate()
        self.hl = self.hLocal.duplicate()

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

    def _matOcean(self):
        """
        This function builds from neighbouring slopes the downstream directions in the marine environment. It calls a fortran subroutine that locally computes for each vertice:

        - the indices of receivers (downstream) nodes depending on the desired number of flow directions (SFD to MFD).
        - the distances to the receivers based on mesh resolution.
        - the associated weights calculated based on the number of receivers and proportional to the slope.

        From these downstream directions, a local ocean downstream matrix is computed.
        """

        # Define multiple flow directions for filled + eps elevations
        rcv, _, wght = mfdrcvrs(8, 1.0e-2, self.oceanFill, -1.0e6)
        sum_wght = np.sum(wght, axis=1)
        ids = (self.pitIDs > -1) & (self.flatOcean > -1) & (sum_wght == 0.0)
        ids = ids.nonzero()[0]
        rcv[ids, :] = np.tile(ids, (8, 1)).T
        rcv[ids, 0] = self.flatOcean[ids]
        wght[ids, :] = 0.0
        wght[ids, 0] = 1.0

        # Get local nodes with no receivers as boolean array
        sum_weight = np.sum(wght, axis=1)
        msink = sum_weight == 0.0

        # We don't consider borders as sinks
        msink[self.idBorders] = False
        msink = msink.astype(int) * self.inIDs
        self.msink = msink == 1

        # Set borders nodes
        if self.flatModel:
            rcv[self.idBorders, :] = np.tile(self.idBorders, (8, 1)).T
            wght[self.idBorders, :] = 0.0

        # Define downstream matrix based on filled + dir elevations
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
        The nonlinear system at each time step is solved iteratively using PETSc time stepping and SNES solution and is based on a Nonlinear Generalized Minimum Residual method (``NGMRES``) .
        
        Here we define the function for the nonlinear solve.
        """

        self.dm.globalToLocal(x, self.hl)
        with self.hl as hl, self.hLocal as zb, xdot as hdot:
            dh = hl - zb
            dh[dh < 0] = 0.0
            Cd = np.multiply(self.Cd, dh / (dh + 5.))
            nlvec = fctcoeff(hl, Cd)
            f.setArray(hdot + nlvec[self.glIDs])

        return

    def _evalJacobian(self, ts, t, x, xdot, a, J, P):
        """
        The nonlinear system at each time step is solved iteratively using PETSc time stepping and SNES solution and is based on a Nonlinear Generalized Minimum Residual method (``NGMRES``) .
         
        Here we define the Jacobian for the nonlinear solve.
        """

        self.dm.globalToLocal(x, self.hl)

        with self.hl as hl, self.hLocal as zb:
            dh = hl - zb
            dh[dh < 0] = 0.0
            Cd = np.multiply(self.Cd, dh / (dh + 5.))

            # Coefficient derivatives
            Cp = np.multiply(self.Cd, 5. / (dh + 5.)**2)
            nlC = jacobiancoeff(hl, Cd, Cp)

            P.zeroEntries()
            for row in range(self.lpoints):
                P.setValuesLocal(row, row, a + nlC[row, 0])
                cols = self.FVmesh_ngbID[row, :]
                P.setValuesLocal(row, cols, nlC[row, 1:])
            P.assemble()

            if J != P:
                J.assemble()

        return petsc4py.PETSc.Mat.Structure.SAME_NONZERO_PATTERN

    def _diffuseOcean(self, dh):
        r"""
        For sediment reaching the marine realm, this function computes the related marine deposition diffusion. The approach is based on a nonlinear diffusion.

        .. math::
          \frac{\partial h}{\partial t}= \nabla \cdot \left( C_d(h) \nabla h \right)

        It calls the following *private functions*:

        - _evalFunction
        - _evalJacobian

        .. note::

            PETSc SNES and time stepping TS approaches are used to solve the nonlinear equation above over the considered time step.


        :arg dh: numpy array of incoming marine depositional thicknesses

        :return: ndepo (updated deposition numpy arrays)
        """

        t0 = process_time()

        # Get diffusion coefficients based on sediment type
        sedK = np.zeros(self.lpoints)
        sedK[self.seaID] = self.nlK

        # Matrix coefficients
        self.Cd = np.full(self.lpoints, sedK, dtype=np.float64)
        self.hl.setArray(dh)
        self.hl.axpy(1.0, self.hLocal)
        self.dm.localToGlobal(self.hl, self.h)

        # Time stepping definition
        ts = petsc4py.PETSc.TS().create(comm=petsc4py.PETSc.COMM_WORLD)
        # ARKIMEX: implicit nonlinear time stepping
        ts.setType("arkimex")
        ts.setIFunction(self._evalFunction, self.tmp1)
        ts.setIJacobian(self._evalJacobian, self.mat)

        ts.setTime(0.0)
        ts.setTimeStep(self.dt / 100.0)
        ts.setMaxTime(self.dt)
        ts.setMaxSteps(50)
        ts.setExactFinalTime(petsc4py.PETSc.TS.ExactFinalTime.MATCHSTEP)
        # Allow an unlimited number of failures (step will be rejected and retried)
        ts.setMaxSNESFailures(-1)
        ts.setTolerances(rtol=1.0e-2)

        # SNES nonlinear solver definition
        snes = ts.getSNES()
        # Newton linear search
        snes.setType("ngmres")
        # Stop nonlinear solve after 10 iterations (TS will retry with shorter step)
        snes.setTolerances(rtol=1.0e-2, max_it=50)

        # KSP linear solver definition
        ksp = snes.getKSP()
        ksp.setType("richardson")
        # Preconditioner for linear solution
        pc = ksp.getPC()
        pc.setType("bjacobi")
        ksp.setTolerances(rtol=1.0e-2, max_it=50)
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
        ndepo[ndepo < 0.0] = 0.0
        self.tmpL.setArray(ndepo)
        self.dm.localToGlobal(self.tmpL, self.tmp)

        return

    def _distOcean(self, sedflux):
        """
        Based on the incoming marine volumes of sediment and maximum clinoforms slope we distribute
        locally sediments downslope.

        :arg sedflux: incoming marine sediment volumes

        :return: seddep (the deposied thickness of the distributed sediments)
        """

        marVol = self.maxDepQs * self.dt
        sinkVol = sedflux * self.dt
        vol = marVol.copy()
        self.tmpL.setArray(sinkVol)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        vdep = np.zeros(self.lpoints, dtype=float)
        step = 0
        while self.tmp.sum() > 1.0:

            # Move to downstream nodes
            self.dMat.mult(self.tmp, self.tmp1)
            self.dm.globalToLocal(self.tmp1, self.tmpL)

            # In case there is too much sediment coming in
            sinkVol = self.tmpL.getArray().copy()
            excess = sinkVol >= vol
            sinkVol[excess] -= vol[excess]
            vdep[excess] = marVol[excess]
            vol[excess] = 0.0

            # In case there is some room to deposit sediments
            noexcess = np.invert(excess)
            vol[noexcess] -= sinkVol[noexcess]
            vdep[noexcess] += sinkVol[noexcess]
            vdep[self.idBorders] = 0.0
            sinkVol[noexcess] = 0.0
            sinkVol[self.idBorders] = 0.0

            # Find where excess and sink
            self.tmpL.setArray(sinkVol)
            self.dm.localToGlobal(self.tmpL, self.tmp)

            sumExcess = self.tmp.sum()
            if MPIrank == 0 and self.verbose:
                if step % 50 == 0:
                    print(
                        "  --- Marine excess (sum) %0.01f m | iter %d"
                        % (sumExcess, step),
                        flush=True
                    )

            step += 1

        self.dMat.destroy()

        return vdep / self.larea

    def seaChange(self):
        """
        This function is the main entry point to perform marine river-induced deposition. It calls the private functions:

        - _distanceCoasts
        - _matOcean
        - _diffuseOcean

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

        # Get the volumetric marine sediment rate (m3/yr) to distribute during the time step.
        self.vSedLocal.copy(result=self.QsL)
        sedFlux = self.QsL.getArray().copy()
        sedFlux[np.invert(self.sinkIDs)] = 0.0
        sedFlux[sedFlux < 0] = 0.0

        # Downstream direction matrix for ocean distribution
        self._matOcean()
        marDep = self._distOcean(sedFlux)
        self._diffuseOcean(marDep)

        # Update cumulative erosion and deposition as well as elevation
        self.cumED.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        # Update erosion/deposition rates
        self.dm.globalToLocal(self.tmp, self.tmpL)
        add_rate = self.tmpL.getArray() / self.dt
        self.tmpL.setArray(add_rate)
        self.EbLocal.axpy(1.0, self.tmpL)

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
