import os
import gc
import sys
import petsc4py
import numpy as np
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import sethillslopecoeff
    from gospl._fortran import hillslp_nl
    from gospl._fortran import jacobiancoeff
    from gospl._fortran import fctcoeff

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD


class hillSLP(object):
    """
    This class encapsulates all the functions related to hillslope (soil creep) processes (both linear and non-linear cases are implemented).
    """

    def __init__(self, *args, **kwargs):
        """
        The initialisation of `hillSLP` class consists in the declaration of several PETSc vectors.
        """

        self.dlim = False
        self.Dlimit = 5.
        self.dexp = 0.05
        self.minDiff = 1.e-4
    
        self.mat = self.dm.createMatrix()
        self.mat.setOption(self.mat.Option.NEW_NONZERO_LOCATIONS, True)

        self.h = self.hGlobal.duplicate()
        self.hl = self.hLocal.duplicate()
        self.dh = self.hGlobal.duplicate()

        return

    def _hillSlope(self, smooth=0):
        r"""
        This function computes hillslope using a **linear** diffusion law commonly referred to as **soil creep**:

        .. math::
          \frac{\partial z}{\partial t}= \kappa_{D} \nabla^2 z

        in which :math:`\kappa_{D}` is the diffusion coefficient and can be defined with different values for the marine and land environments (set with `hillslopeKa` and `hillslopeKm` in the YAML input file). It encapsulates, in a simple formulation, processes operating on superficial sedimentary layers. Main controls on variations of :math:`\kappa_{D}` include substrate, lithology, soil depth, climate and biological activity.

        .. note::
            The hillslope processes in goSPL are considered to be happening at the same rate for coarse and fine sediment sizes.

        :arg smooth: integer specifying if the diffusion equation is used for ice flow (1) and marine deposits (2).
        """

        if smooth == 0:
            if self.Cda == 0.0 and self.Cdm == 0.0:
                return

        t0 = process_time()
        # Diffusion matrix construction
        if smooth == 1:
            Cd = np.full(self.lpoints, self.gaussIce, dtype=np.float64)
            Cd[~self.iceIDs] = 0.0
        elif smooth == 2:
            # Hard-coded coefficients here, used to generate a smooth surface
            # for computing marine flow directions...
            Cd = np.full(self.lpoints, 1.e5, dtype=np.float64)
            Cd[self.seaID] = 5.e6
        else:
            Cd = np.full(self.lpoints, self.Cda, dtype=np.float64)
            Cd[self.seaID] = self.Cdm
        diffCoeffs = sethillslopecoeff(self.lpoints, Cd * self.dt)
        if self.flatModel:
            diffCoeffs[self.idBorders, 1:] = 0.0
            diffCoeffs[self.idBorders, 0] = 1.0

        diffMat = self._matrix_build_diag(diffCoeffs[:, 0])
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)

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
            )
            tmpMat.assemblyEnd()
            diffMat.axpy(1.0, tmpMat)
            tmpMat.destroy()

        # Get elevation values for considered time step
        if smooth == 1:
            if self.tmp1.max()[1] > 0:
                self._solve_KSP(True, diffMat, self.tmp1, self.tmp)
            else:
                self.tmp1.copy(result=self.tmp)
            diffMat.destroy()
            self.dm.globalToLocal(self.tmp, self.tmpL)
            return self.tmpL.getArray().copy()
        elif smooth == 2:
            self._solve_KSP(True, diffMat, self.hGlobal, self.tmp)
            diffMat.destroy()
            self.dm.globalToLocal(self.tmp, self.tmpL)
            return self.tmpL.getArray().copy()
        else:
            self.hGlobal.copy(result=self.hOld)
            self._solve_KSP(True, diffMat, self.hOld, self.hGlobal)
            diffMat.destroy()
            # Update cumulative erosion/deposition and elevation
            self.tmp.waxpy(-1.0, self.hOld, self.hGlobal)
            self.cumED.axpy(1.0, self.tmp)
            self.dm.globalToLocal(self.cumED, self.cumEDLocal)
            self.dm.globalToLocal(self.hGlobal, self.hLocal)

            if self.memclear:
                del ids, indices, indptr, diffCoeffs, Cd
                gc.collect()

            if self.stratNb > 0:
                self.erodeStrat()
                self.deposeStrat()

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Linear Hillslope Processes (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )
        petsc4py.PETSc.garbage_cleanup()

        return

    def _diff_nl_monitor(self, snes, its, norm):
        """
        Non-linear diffusion solver convergence evaluation.
        """

        if MPIrank == 0 and its % 5 == 0:
            print(f"  ---  Non-linear hillslope solver iteration {its}, Residual norm: {norm}", flush=True)

    def _form_residual_nl_hillslope(self, snes, h, F):
        """
        The nonlinear system (SNES) at each time step is solved iteratively by assessing the residual of the hillslope equation.

        Parameters:
        -----------
        snes : PETSc.SNES: The snes object.
        h : PETSc.Vec: The current solution vector (h^{n+1}) at the new time step.
        F : PETSc.Vec: The residual vector to be filled.
        """

        # Current state
        self.dm.globalToLocal(h, self.hl)
        h_array = self.hl.getArray()

        # Compute slope
        if self.K_nb == 0:
            val = hillslp_nl(self.lpoints, h_array, self.Cd_nl, self.K_nl, 0)
        else:
            val = hillslp_nl(self.lpoints, h_array, self.Cd_nl, self.K_sc, self.K_nb)

        if self.flatModel:
            val[self.idBorders] = 0.

        # Residuals
        res = h_array - self.hOldArray - self.dt * val

        # Residual vector
        F.setArray(res[self.glIDs])

        return
    
    def _hillSlopeNL(self):
        r"""
        This function computes hillslope using a **non-linear** diffusion law:

        .. math::
          \frac{\partial h}{\partial t}= \nabla \cdot \left( C_d(h) \nabla h \right)

        Two flavors of non-linear diffusion are possible based on user defined parameters:

        1. a non-critical hillslope model following the work of `Wang et al. (2024) <https://www.sciencedirect.com/science/article/pii/S0169555X24001053>`_.
        2. a non-linear depth-dependent creep law as described in `Barnhart et al. (2019) <https://gmd.copernicus.org/articles/12/1267/2019/gmd-12-1267-2019.pdf>`_ (section 3.4.5).

        .. note::
            In this implementation of the SNES, we do not form the Jacobian and PETSc calculates it based on the residual function. Here, a Nonlinear Generalized Minimum Residual method is used ``ngmres``, a Preconditioned Conjugate Gradient ``cg`` method is defined for the KSP and the preconditioner allows for multi-grid methods based on the HYPRE BoomerAMG approach.
        """

        t0 = process_time()

        self.hGlobal.copy(result=self.hOld)
        self.Cd_nl = np.full(self.lpoints, self.Cda, dtype=np.float64)
        self.Cd_nl[self.seaID] = self.Cdm
        self.hOldArray = self.hLocal.getArray().copy()

        snes = petsc4py.PETSc.SNES().create(comm=petsc4py.PETSc.COMM_WORLD)
        snes.setTolerances(rtol=self.snes_rtol, atol=self.snes_atol,
                           max_it=self.snes_maxit)
        
        # Set a monitor to see residual values
        if self.verbose:
            snes.setMonitor(self._diff_nl_monitor)

        f = self.hGlobal.duplicate()
        snes.setFunction(self._form_residual_nl_hillslope, f)

        # SNES solvers
        snes.setType('ngmres')
        ksp = snes.getKSP()
        ksp.setType("cg")
        pc = ksp.getPC()
        pc.setType("hypre")
        petsc_options = petsc4py.PETSc.Options()
        petsc_options['pc_hypre_type'] = 'boomeramg'
        ksp.getPC().setFromOptions()
        ksp.setTolerances(rtol=1.0e-6)

        b, x = None, f.duplicate()
        self.hGlobal.copy(result=x)
        snes.solve(b, x)

        # Clean solver
        pc.destroy()
        ksp.destroy()
        snes.destroy()

        x.copy(result=self.hGlobal)
        f.destroy()
        x.destroy()

        # Update cumulative erosion/deposition and elevation
        self.tmp.waxpy(-1.0, self.hOld, self.hGlobal)
        self.cumED.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        if self.stratNb > 0:
            self.erodeStrat()
            self.deposeStrat()

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Non-Linear Hillslope Processes (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )
        petsc4py.PETSc.garbage_cleanup()

        return

    def getHillslope(self):
        """
        This function is the main entry point to compute linear and non-linear hillslope processes. 
        """

        h = self.hLocal.getArray().copy()
        self.seaID = np.where(h <= self.sealevel)[0]

        # Specify the hillslope diffusion law to use
        if self.cptSoil:
            self.diffuseSoil()
            return
        
        if self.K_nl == 1.0 and self.K_nb == 0:
            self._hillSlope(smooth=0)
        else:
            self._hillSlopeNL()

        # Update layer elevation
        if self.stratNb > 0:
            self.elevStrat()
        if self.memclear:
            del h
            gc.collect()

        # Update erosion/deposition rates
        self.dm.globalToLocal(self.tmp, self.tmpL)
        add_rate = self.tmpL.getArray() / self.dt
        self.tmpL.setArray(add_rate)
        self.EbLocal.axpy(1.0, self.tmpL)

        if self.memclear:
            del add_rate
            gc.collect()

        return

    def _evalFunctionMardDiff(self, ts, t, x, xdot, f):
        """
        The non-linear system for freshly-deposited marine sediments diffusion at each time step is solved iteratively using PETSc time stepping and SNES solution and is based on Rosenbrock W-scheme (``rosw``).

        Here we evaluate the residual function on a DMPlex for an implicit time-stepping method.

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

    def _evalJacobianMardDiff(self, ts, t, x, xdot, a, A, B):
        """
        The non-linear system for freshly-deposited marine sediments diffusion at each time step is solved iteratively using PETSc time stepping and SNES solution and is based on Rosenbrock W-scheme (``rosw``).

        Here, we define the Jacobian matrix J and the preconditioner matrix P on a DMPlex.

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

    def _evalSolutionMardDiff(self, t, x):
        """
        Evaluate the initial solution of the SNES system.
        """

        assert t == 0.0, "only for t=0.0"
        x.setArray(self.h.getArray())

        return

    def _diffuseOcean(self, dh):
        r"""
        For river-transported sediments reaching the marine realm, this function computes the related marine deposition diffusion. It is based on a non-linear diffusion approach.

        .. math::
          \frac{\partial h}{\partial t}= \nabla \cdot \left( C_d(h) \nabla h \right)

        It calls the following *private functions*:

        - _evalFunctionMardDiff
        - _evalJacobianMardDiff
        - _evalSolutionMardDiff

        .. note::

            PETSc SNES and time stepping TS approaches are used to solve the non-linear equation above over the considered time step.


        :arg dh: Numpy Array of incoming marine depositional thicknesses
        """

        t0 = process_time()

        x = self.tmp1.duplicate()
        f = self.tmp1.duplicate()

        # Get diffusion coefficients based on sediment type
        self.Cd = np.zeros(self.lpoints)
        self.Cd[self.seaID] = self.nlK
        self.hl.setArray(dh)
        self.hl.axpy(1.0, self.hLocal)
        self.dm.localToGlobal(self.hl, self.h)

        # Time stepping definition
        ts = petsc4py.PETSc.TS().create(comm=petsc4py.PETSc.COMM_WORLD)
        # arkimex: IMEX Runge-Kutta schemes | rosw: Rosenbrock W-schemes
        ts.setType("rosw")

        ts.setIFunction(self._evalFunctionMardDiff, self.tmp1)
        ts.setIJacobian(self._evalJacobianMardDiff, self.mat)

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
        self._evalSolutionMardDiff(tstart, x)

        # Solve nonlinear equation
        ts.solve(x)
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
                ),
                flush=True,
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
            del ndepo
            gc.collect()
        petsc4py.PETSc.garbage_cleanup()

        return
