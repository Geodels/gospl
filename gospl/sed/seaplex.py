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
    from gospl._fortran import mfdreceivers
    from gospl._fortran import donorslist
    from gospl._fortran import donorsmax
    from gospl._fortran import mfdrcvrs
    from gospl._fortran import jacobiancoeff
    from gospl._fortran import fctcoeff
    from gospl._fortran import epsfill

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

        self.dlim = False
        self.Dlimit = 5.
        self.dexp = 0.05
        self.minDiff = 1.e-4
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
        hl = self.hLocal.getArray().copy()
        if not self.flatModel:
            # Only consider filleps in the first kms offshore
            hsmth = self._hillSlope(smooth=2)
            hsmth[self.coastDist > self.offshore] = -1.e6
        else:
            hsmth = hl.copy()

        fillz = np.zeros(self.mpoints, dtype=np.float64) - 1.0e8
        fillz[self.locIDs] = hsmth
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, fillz, op=MPI.MAX)
        if MPIrank == 0:
            minh = np.min(fillz) + 0.1
            if not self.flatModel:
                minh = min(minh, self.oFill)
            fillz = epsfill(minh, fillz)
        # Send elevation + eps globally
        fillEPS = MPI.COMM_WORLD.bcast(fillz, root=0)
        fillz = fillEPS[self.locIDs]
        if not self.flatModel:
            fillz[self.coastDist > self.offshore] = hl[self.coastDist > self.offshore]
        rcv, _, wght = mfdrcvrs(12, self.flowExp, fillz, -1.0e6)

        # Set borders nodes
        if self.flatModel:
            rcv[self.idBorders, :] = np.tile(self.idBorders, (12, 1)).T
            wght[self.idBorders, :] = 0.0

        # Define downstream matrix based on filled + dir elevations
        self.dMat1 = self.zMat.copy()
        if not self.flatModel and self.Gmar > 0.:
            self.dMat2 = self.iMat.copy()
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        nodes = indptr[:-1]

        for k in range(0, 12):
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
            self.dMat1.axpy(1.0, tmpMat)
            if not self.flatModel and self.Gmar > 0.:
                self.dMat2.axpy(-1.0, tmpMat)
            tmpMat.destroy()

        if self.memclear:
            del data, indptr, nodes
            del hl, fillz, fillEPS, rcv, wght
            gc.collect()

        # Store flow direction matrix
        self.dMat1.transpose()
        if not self.flatModel and self.Gmar > 0.:
            self.dMat2.transpose()

        return

    def _evalFunction(self, ts, t, x, xdot, f):
        """
        The nonlinear system at each time step is solved iteratively using PETSc time stepping and SNES solution and is based on a Nonlinear Generalized Minimum Residual method (``NGMRES``) .

        Here we define the function for the nonlinear solve.
        Evaluate the residual function on a DMPlex for an implicit time-stepping method.

        Parameters:
        -----------
        ts : PETSc.TS: The time-stepper object.
        t : float: The current time.
        x : PETSc.Vec: The current solution vector (h^{n+1}) at the new time step.
        xdot : PETSc.Vec: The time derivative approximation (h^{n+1} - h^n) / dt.
        f : PETSc.Vec: The residual vector to be filled.
        """

        self.dm.globalToLocal(x, self.hl)
        with self.hl as hl, self.hLocal as zb, xdot as hdot:
            dh = hl - zb
            dh[dh < 0.1] = 0.0
            if self.dlim:
                Cd = self.minDiff + np.multiply(self.Cd, dh / (dh + self.Dlimit))
            else:
                Cd = self.minDiff + np.multiply(self.Cd, (1.0 - np.exp(-self.dexp * dh)))
            nlvec = fctcoeff(hl, Cd)
            f.setArray(hdot + nlvec[self.glIDs])

        return

    def _evalJacobian(self, ts, t, x, xdot, a, A, B):
        """
        The nonlinear system at each time step is solved iteratively using PETSc time stepping and SNES solution and is based on a Nonlinear Generalized Minimum Residual method (``NGMRES``) .

        Here we define the Jacobian for the nonlinear solve.

        Evaluate the Jacobian matrix J and the preconditioner matrix P on a DMPlex.

        Parameters:
        -----------
        ts : PETSc.TS: The time-stepper object.
        t : float: The current time.
        x : PETSc.Vec: The current solution vector (h^{n+1}) at the new time step.
        xdot : PETSc.Vec: The time derivative approximation (h^{n+1} - h^n) / dt.
        a : float:  The shift factor for implicit methods.
        A : PETSc.Mat: The Jacobian matrix to be filled.
        B : PETSc.Mat: The preconditioner matrix to be filled.

        """

        self.dm.globalToLocal(x, self.hl)

        with self.hl as hl, self.hLocal as zb:
            dh = hl - zb
            dh[dh < 0.1] = 0.0
            if self.dlim:
                Cd = self.minDiff + np.multiply(self.Cd, dh / (dh + self.Dlimit))
            else:
                Cd = self.minDiff + np.multiply(self.Cd, (1.0 - np.exp(-self.dexp * dh)))

            # Coefficient derivatives
            if self.dlim:
                Cp = np.multiply(self.Cd, self.Dlimit / (dh + self.Dlimit)**2)
            else:
                Cp = np.multiply(self.Cd, self.dexp * np.exp(-self.dexp * dh))
            nlC = jacobiancoeff(hl, Cd, Cp)

            for row in range(self.lpoints):
                B.setValuesLocal(row, row, a + nlC[row, 0])
                cols = self.FVmesh_ngbID[row, :]
                B.setValuesLocal(row, cols, nlC[row, 1:])
            B.assemble()

            if A != B:
                A.assemble()

        return True

    def _evalSolution(self, t, x):

        assert t == 0.0, "only for t=0.0"
        x.setArray(self.h.getArray())

        return

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

        x = self.tmp1.duplicate()
        f = self.tmp1.duplicate()

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
        # arkimex: IMEX Runge-Kutta schemes | rosw: Rosenbrock W-schemes
        ts.setType("rosw")

        ts.setIFunction(self._evalFunction, self.tmp1)
        ts.setIJacobian(self._evalJacobian, self.mat)

        ts.setTime(0.0)
        ts.setTimeStep(self.dt / 1000.0)
        ts.setMaxTime(self.dt)
        ts.setMaxSteps(self.tsStep)
        ts.setExactFinalTime(petsc4py.PETSc.TS.ExactFinalTime.MATCHSTEP)

        # Allow an unlimited number of failures
        ts.setMaxSNESFailures(-1)  # (step will be rejected and retried)

        # SNES nonlinear solver
        snes = ts.getSNES()
        snes.setTolerances(max_it=10)   # Stop nonlinear solve after 10 iterations (TS will retry with shorter step)

        # KSP linear solver
        ksp = snes.getKSP()
        ksp.setType("preonly")
        pc = ksp.getPC()
        pc.setType("gasm")

        ts.setFromOptions()
        tstart = ts.getTime()
        self._evalSolution(tstart, x)

        # Solve nonlinear equation
        ts.solve(x)
        # te = ts.getTime()
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

        # Clean solver
        pc.destroy()
        ksp.destroy()
        snes.destroy()
        ts.destroy()

        # Get diffused sediment thicknesses
        x.copy(result=self.h)
        self.dh.waxpy(-1.0, self.hGlobal, self.h)
        self.dm.globalToLocal(self.dh, self.tmpL)
        ndepo = self.tmpL.getArray().copy()
        self.tmpL.setArray(ndepo)
        self.dm.localToGlobal(self.tmpL, self.tmp)

        x.destroy()
        f.destroy()
        if self.memclear:
            del ndepo, sedK
            gc.collect()

        return

    def _depMarineSystem(self, sedflux):
        r"""
        Setup matrix for the marine sediment deposition.

        The upstream incoming sediment flux is obtained from the total river sediment flux :math:`\mathrm{Q_{t_i}}` where:

        .. math::

            \mathrm{Q_{t_i}^{t+\Delta t} - \sum_{ups} w_{i,j} Q_{t_u}^{t+\Delta t}}= \mathrm{(\eta_i^{t} - \eta_i^{t+\Delta t}) \frac{\Delta t}{\Omega_i}}

        which gives:

        .. math::

            \mathrm{Q_{s_i}} = \mathrm{Q_{t_i}} - \mathrm{(\eta_i^{t} - \eta_i^{t+\Delta t}) \frac{\Delta t}{\Omega_i}}

        And the evolution of marine elevation is based on incoming sediment flux resulting.

        .. math::

            \mathrm{\frac{\eta_i^{t+\Delta t}-\eta_i^t}{\Delta t}} = \mathrm{G{_m} Q_{s_i} / \Omega_i}

        This system of coupled equations is solved implicitly using PETSc by assembling the matrix and vectors using the nested submatrix and subvectors and by using the ``fieldsplit`` preconditioner combining two separate preconditioners for the collections of variables.

        :arg sedflux: incoming marine sediment volumes

        :return: volDep (the deposited volume of the distributed sediments)
        """

        hl = self.hLocal.getArray()
        fDepm = np.full(self.lpoints, self.Gmar)
        fDepm[fDepm > 0.99] = 0.99
        fDepm[hl > self.sealevel] = 0.

        # Define submatrices
        A00 = self._matrix_build_diag(-fDepm)
        A00.axpy(1.0, self.iMat)
        A01 = self._matrix_build_diag(-fDepm * self.dt / self.larea)
        A10 = self._matrix_build_diag(self.larea / self.dt)

        # Assemble the matrix for the coupled system
        mats = [[A00, A01], [A10, self.dMat2]]
        sysMat = petsc4py.PETSc.Mat().createNest(mats=mats, comm=MPIcomm)
        sysMat.assemblyBegin()
        sysMat.assemblyEnd()

        # Clean up
        A00.destroy()
        A01.destroy()
        A10.destroy()
        self.dMat2.destroy()
        mats[0][0].destroy()
        mats[0][1].destroy()
        mats[1][0].destroy()
        mats[1][1].destroy()

        # Create nested vectors
        self.tmpL.setArray(1. - fDepm)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.tmp.pointwiseMult(self.tmp, self.hGlobal)

        self.tmpL.setArray(sedflux / self.dt)
        self.dm.localToGlobal(self.tmpL, self.tmp1)
        self.h.pointwiseMult(self.hGlobal, self.areaGlobal)
        self.h.scale(1. / self.dt)
        self.tmp1.axpy(1., self.h)

        rhs_vec = petsc4py.PETSc.Vec().createNest([self.tmp, self.tmp1], comm=MPIcomm)
        rhs_vec.setUp()
        hq_vec = rhs_vec.duplicate()

        # Define solver and precondition conditions
        ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
        ksp.setType(petsc4py.PETSc.KSP.Type.TFQMR)
        ksp.setOperators(sysMat)
        ksp.setTolerances(rtol=self.rtol)

        pc = ksp.getPC()
        pc.setType("fieldsplit")
        nested_IS = sysMat.getNestISs()
        pc.setFieldSplitIS(('h', nested_IS[0][0]), ('q', nested_IS[0][1]))

        subksps = pc.getFieldSplitSubKSP()
        subksps[0].setType("preonly")
        subksps[0].getPC().setType("asm")
        subksps[1].setType("preonly")
        subksps[1].getPC().setType("bjacobi")

        ksp.solve(rhs_vec, hq_vec)
        r = ksp.getConvergedReason()
        if r < 0:
            KSPReasons = self._make_reasons(petsc4py.PETSc.KSP.ConvergedReason())
            if MPIrank == 0:
                print(
                    "Linear solver for marine deposition failed to converge after iterations",
                    ksp.getIterationNumber(),
                    flush=True,
                )
                print("with reason: ", KSPReasons[r], flush=True)
        else:
            if MPIrank == 0 and self.verbose:
                print(
                    "Linear solver for marine deposition converge after %d iterations"
                    % ksp.getIterationNumber(),
                    flush=True,
                )

        # Update the solution
        self.newH = hq_vec.getSubVector(nested_IS[0][0])

        # Clean up
        subksps[0].destroy()
        subksps[1].destroy()
        nested_IS[0][0].destroy()
        nested_IS[1][0].destroy()
        nested_IS[0][1].destroy()
        nested_IS[1][1].destroy()
        pc.destroy()
        ksp.destroy()
        sysMat.destroy()
        hq_vec.destroy()
        rhs_vec.destroy()

        # Get the marine deposition volume
        self.tmp.waxpy(-1.0, self.hGlobal, self.newH)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        volDep = self.tmpL.getArray().copy() * self.larea
        volDep[volDep < 0] = 0.

        return volDep

    def _distOcean(self, sedflux):
        """
        Based on the incoming marine volumes of sediment and maximum clinoforms slope we distribute
        locally sediments downslope.

        :arg sedflux: incoming marine sediment volumes

        :return: vdep (the deposited volume of the distributed sediments)
        """

        marVol = self.maxDepQs.copy()
        sinkVol = sedflux.copy()
        self.tmpL.setArray(sinkVol)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        vdep = np.zeros(self.lpoints, dtype=float)

        step = 0
        while self.tmp.sum() > 1.0:

            # Move to downstream nodes
            self.dMat1.mult(self.tmp, self.tmp1)
            self.dm.globalToLocal(self.tmp1, self.tmpL)

            # In case there is too much sediment coming in
            sinkVol = self.tmpL.getArray().copy()
            excess = sinkVol >= marVol
            sinkVol[excess] -= marVol[excess]
            vdep[excess] += marVol[excess]
            marVol[excess] = 0.0

            # In case there is some room to deposit sediments
            noexcess = np.invert(excess)
            marVol[noexcess] -= sinkVol[noexcess]
            vdep[noexcess] += sinkVol[noexcess]
            vdep[self.idBorders] = 0.0
            sinkVol[noexcess] = 0.0
            sinkVol[self.idBorders] = 0.0

            # Find where excess and sink
            self.tmpL.setArray(sinkVol)
            self.dm.localToGlobal(self.tmpL, self.tmp)

            sumExcess = self.tmp.sum()
            if MPIrank == 0 and self.verbose:
                if step % 100 == 0:
                    print(
                        "  --- Marine excess (sum in km3) %0.05f | iter %d"
                        % (sumExcess * 10.e-9, step),
                        flush=True
                    )

            step += 1

        if self.memclear:
            del marVol, sinkVol
            gc.collect()

        return vdep

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
        hl = self.hLocal.getArray().copy()
        if self.clinSlp > 0.0:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                self._distanceCoasts(hl)
            # From the distance to coastline define the upper limit of the shelf to ensure a maximum slope angle
            self.clinoH = self.sealevel - 1.0e-3 - self.coastDist * self.clinSlp
        else:
            self.clinoH = np.full(self.lpoints, self.sealevel - 1.0e-3, dtype=np.float64)
        self.clinoH[hl >= self.sealevel] = hl[hl >= self.sealevel]
        self.clinoH[self.clinoH < hl] = hl[self.clinoH < hl]
        self.maxDepQs = (self.clinoH - hl) * self.larea

        # Get the volumetric marine sediment (m3) to distribute during the time step.
        self.vSedLocal.copy(result=self.QsL)
        sedFlux = self.QsL.getArray().copy() * self.dt
        sedFlux[np.invert(self.sinkIDs)] = 0.0
        sedFlux[sedFlux < 0] = 0.0

        # Downstream direction matrix for ocean distribution
        self._matOcean()
        marDep = self._distOcean(sedFlux)
        if not self.flatModel and self.Gmar > 0.:
            vdep = self._depMarineSystem(marDep)
            marDep = self._distOcean(vdep)
        self.dMat1.destroy()

        # Diffuse downstream
        dh = np.divide(marDep, self.larea, out=np.zeros_like(self.larea), where=self.larea != 0)
        dh[dh < 1.e-3] = 0.
        self.tmpL.setArray(dh)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        dh = self.tmpL.getArray().copy()
        self._diffuseOcean(dh)

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

        if self.memclear:
            del marDep, dh, add_rate, hl, sedFlux
            gc.collect()

        return
