import os
import gc
import sys
import petsc4py
import numpy as np
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import fctcoeff
    from gospl._fortran import local_spl
    from gospl._fortran import jacobiancoeff

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD


class soilSPL(object):
    """
    The class computes river incision expressed using a **stream power formulation** function of river discharge and slope also **accounting for soil production**.

    A non-linear diffusion of soil based on soil thickness is also implemented in this class.

    If the user has turned-on the sedimentation capability, this class will solve implicitly the **stream power formulation** accounting for a sediment transport/deposition term (`Yuan et al, 2019 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JF004867>`_).
    """

    def __init__(self, *args, **kwargs):
        """
        Initialisation of `soilSPL` class.
        """

        self.soil_rtol = 1.0e-6
        self.soil_atol = 1.e-6
        self.soil_maxit = 500

        # self.cptSoil = True
        # self.Ksoil = 2 * self.K # erodibility coefficient for soil
        # self.P0 = 50. * 1.e-6 # soil production maximum rate (50 m/Myr)
        # self.Hs = 0.5 # soil production decay depth
        # self.h_star = 1.0 # roughness length_scale
        # self.H0 = 0.7 # soil transport decay depth for diffusion
        # self.Sperc = 0.0001 # soil / bedrock transition limit ratio factor of production

        if self.Sperc > 0:
            self.soil_transition = -np.log(self.Sperc) * self.Hs
        else:
            self.soil_transition = 100.0
        self.Gsoil = self.hGlobal.duplicate()
        self.Lsoil = self.hLocal.duplicate()
        self.lHbed = self.hLocal.duplicate()
        self.gHbed = self.hGlobal.duplicate()

        # Allocate initial soil thicknesses
        if self.soilFile is not None:
            loadData = np.load(self.soilFile)
            soilH = loadData[self.soilData]
            self.Lsoil.setArray(soilH[self.locIDs])
            self.dm.localToGlobal(self.Lsoil, self.Gsoil)
        else:
            # Create an uniform soil thickness distribution
            self.Gsoil.set(self.cstSoilH)
            self.Lsoil.set(self.cstSoilH)

        return

    def _form_residual_soil(self, snes, h, F):
        """
        The nonlinear system (SNES) at each time step is solved iteratively by assessing the residual of the SPL equation accounting for erosion, deposition (transport-limited) and soil production.

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
        S = local_spl(self.flowDir, h_array, self.rcvIDi,
                      self.distRcvi, self.wghtVali)
        S[S < 0.] = 0.

        # Compute upstream sediment flux
        self.tmp.waxpy(-1.0, h, self.hOld)
        self.tmp1.pointwiseMult(self.tmp, self.areaGlobal)
        self.tmp1.scale(1. / self.dt)
        self._solve_KSP(True, self.fMati, self.tmp1, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        Qt = self.tmpL.getArray()
        Qt[Qt < 0.] = 0.

        # Compute soil thickness
        hSoil = self.soilH + h_array - self.hOldArray
        hSoil[hSoil < 0.] = 0.

        # Residuals based on the equation: 
        # h(t+dt) (1-G) - h(t) (1-G) + dt * K * A^m * S^n - dt * G * Qt / Area = 0
        res = (h_array - self.hOldArray) * (1.0 - self.fDep)  
        res += self.Kbr * (1.0 - np.exp(-hSoil / self.h_star) ) * S**self.spl_n 
        res += self.K_soil * np.exp(-hSoil / self.h_star) * S**self.spl_n
        res -= self.fDep * self.dt * Qt / self.larea

        # Residual vector
        F.setArray(res[self.glIDs])

        return

    def _monitorsoil(self, snes, its, norm):
        """
        Non-linear SPL with soil production solver convergence evaluation.
        """

        if MPIrank == 0 and its % 5 == 0:
            print(f"  ---  Non-linear soil SPL solver iteration {its}, Residual norm: {norm}", flush=True)

    def _solveSoil(self):
        """
        Solves the non-linear stream power law for the transport limited and soil case. This calls the following *private function*:

        - _form_residual_soil

        .. note::

            PETSc SNES approach is used to solve the nonlinear equation. In this implementation of the SNES, we do not form the Jacobian and PETSc calculates it based on the residual function. A Nonlinear Generalized Minimum Residual method is used ``ngmres``, a Preconditioned Conjugate Gradient ``cg`` method is defined for the KSP and the preconditioner allows for multi-grid methods based on the HYPRE BoomerAMG approach.
        """

        self.oldH = self.hGlobal.getArray()
        self.hOldArray = self.hLocal.getArray().copy()

        # Apply production rate to existing soil
        self.soilH = self.Lsoil.getArray().copy()
        self.soilH += self.dt * self.P0 * np.exp(-self.soilH / self.Hs)

        # Upstream-averaged mean annual precipitation rate based on drainage area
        PA = self.FAL.getArray()

        # Define erosion limiter to prevent formation of flat
        dh = []
        for k in range(0, self.flowDir):
            dh.append(self.hOldArray - self.hOldArray[self.rcvIDi[:, k]])
        dh = np.array(dh).max(0)
        elimiter = np.divide(dh, dh + 1.0e-2, out=np.zeros_like(dh),
                             where=dh != 0)

        # Incorporate the effect of local mean annual precipitation rate on erodibility (for soil and bedrock)
        if self.sedfacVal is not None:
            self.Kbr = self.K * self.sedfacVal * (self.rainVal ** self.coeffd)
        else:
            self.Kbr = self.K * (self.rainVal ** self.coeffd)
        self.K_soil = self.Ksoil * (self.rainVal ** self.coeffd)
        self.Kbr *= self.dt * (PA ** self.spl_m) * elimiter
        self.K_soil *= self.dt * (PA ** self.spl_m) * elimiter
        self.Kbr[self.seaID] = 0.0
        self.K_soil[self.seaID] = 0.0

        # Dimensionless depositional coefficient
        self.fDep = np.divide(self.fDepa * self.larea, PA, out=np.zeros_like(PA), where=PA != 0)
        self.fDep[self.seaID] = 0.
        self.fDep[self.fDep > 0.99] = 0.99
        if self.flatModel:
            self.fDep[self.idBorders] = 0.

        snes = petsc4py.PETSc.SNES().create(comm=petsc4py.PETSc.COMM_WORLD)
        snes.setTolerances(rtol=self.soil_rtol, atol=self.soil_atol,
                           max_it=self.soil_maxit)

        # Set a monitor to see residual values
        if self.verbose:
            snes.setMonitor(self._monitorsoil)

        f = self.hGlobal.duplicate()
        snes.setFunction(self._form_residual_soil, f)

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

        # Get eroded sediment thicknesses
        self.tmp.waxpy(-1.0, self.hOld, x)
        f.destroy()
        x.destroy()

        # Update soil thicknesses
        self.dm.globalToLocal(self.tmp, self.tmpL)
        nHsoil = self.soilH +  self.tmpL.getArray()
        nHsoil[nHsoil < 0.] = 0.
        # Limit soil thickness
        nHsoil[nHsoil > self.soil_transition] = self.soil_transition
        self.Lsoil.setArray(nHsoil)
        self.dm.localToGlobal(self.Lsoil, self.Gsoil)

        petsc4py.PETSc.garbage_cleanup()

        return
    
    def _getEroDepRateSoil(self):
        """
        This function computes erosion deposition rates in metres per year and associated soil evolution. This is done on the filled elevation. 

        The approach is based on **BasicHySa** governing equations from Terrainbento (as described in Appendix B20 from `Barnhart et al. (2019) <https://gmd.copernicus.org/articles/12/1267/2019/gmd-12-1267-2019.pdf>`_).

        .. note::

           The approach uses a continuous layer of soil–alluvium, which influences both hillslope and river-induced erosion. It relies on the SPACE algorithm of `Shobe et al. (2017) <https://gmd.copernicus.org/articles/10/4577/2017/>`_.
        """

        t0 = process_time()

        # Build the SPL erosion arrays
        if self.flexOn:
            self.hLocal.copy(result=self.hOldFlex)

        self._solveSoil()

        # Update erosion rate (positive for incision)
        E = -self.tmp.getArray().copy()
        E = np.divide(E, self.dt)
        self.Eb.setArray(E)
        self.dm.globalToLocal(self.Eb, self.EbLocal)
        E = self.EbLocal.getArray().copy()
        if self.flatModel:
            E[self.idBorders] = 0.0
        E[self.lsink] = 0.0
        self.EbLocal.setArray(E)
        self.dm.localToGlobal(self.EbLocal, self.Eb)

        if MPIrank == 0 and self.verbose:
            print(
                "Finalise erosion deposition rates (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def updateSoilThickness(self):
        """
        Updates soil thickness through time.
        """

        self.dm.globalToLocal(self.tmp, self.tmpL)
        nHsoil = self.Lsoil.getArray() + self.tmpL.getArray()

        # Limit soil thickness
        nHsoil[nHsoil < 0.] = 0.
        nHsoil[nHsoil > self.soil_transition] = self.soil_transition

        self.Lsoil.setArray(nHsoil)
        self.dm.localToGlobal(self.Lsoil, self.Gsoil)

        return

    def erodepSPLsoil(self):
        """
        Modified **stream power law** model used to represent erosion by rivers also taking into account the role played by sediments in modulating erosion and deposition rate, considering **non-linear slope dependency** and accounting for soil production.

        It calls the private function `_getEroDepRateSoil` described above. Once erosion/deposition rates have been calculated, the function computes local thicknesses and soil evolution for the considered time step and update local elevation and cumulative erosion, deposition values.
        """

        t0 = process_time()

        # Computes the erosion deposition rates based on flow accumulation
        self.Eb.set(0.0)
        self.hGlobal.copy(result=self.hOld)
        self.dm.globalToLocal(self.hOld, self.hOldLocal)
        self._getEroDepRateSoil()

        # Get erosion / deposition thicknesses
        Eb = self.Eb.getArray().copy()
        self.tmp.setArray(-Eb * self.dt)
        self.cumED.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)
        self.tmp1.pointwiseMult(self.tmp, self.areaGlobal)

        # Update stratigraphic layers
        if self.stratNb > 0:
            self.erodeStrat()
            self.deposeStrat()

        # Update erosion/deposition rates
        self.dm.globalToLocal(self.tmp, self.tmpL)
        add_rate = self.tmpL.getArray() / self.dt
        self.EbLocal.setArray(add_rate)

        # Destroy flow matrices
        self.fMati.destroy()
        self.fMat.destroy()

        if MPIrank == 0 and self.verbose:
            print(
                "Get Erosion Deposition values (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        if self.memclear:
            del Eb
            gc.collect()

        return

    def _evalFunctionSoil(self, ts, t, x, xdot, f):
        """
        The non-linear system for soil diffusion is solved iteratively using PETSc time stepping and SNES solution and is based on Rosenbrock W-scheme (``rosw``).

        Here again, we evaluate the residual function on a DMPlex for an implicit time-stepping method.

        Parameters:
        -----------
        ts : PETSc.TS: The time-stepper object.
        t : float: The current time.
        x : PETSc.Vec: The current solution vector (h^{n+1}) at the new time step.
        xdot : PETSc.Vec: The time derivative approximation (h^{n+1} - h^n) / dt.
        f : PETSc.Vec: The residual vector to be filled.
        """

        self.dm.globalToLocal(x, self.hl)
        with self.hl as hl, self.lHbed as zb, xdot as hdot:
            dh = hl - zb
            dh[dh < 0.1] = 0.0
            Cd = self.minDiff + np.multiply(self.Cd, (1.0 - np.exp(-dh/self.H0)))
            nlvec = fctcoeff(hl, Cd)
            f.setArray(hdot + nlvec[self.glIDs])

        return

    def _evalJacobianSoil(self, ts, t, x, xdot, a, A, B):
        """
        The non-linear system for soil diffusion is solved iteratively using PETSc time stepping and SNES solution and is based on Rosenbrock W-scheme (``rosw``).

        Here, we define the Jacobian matrix A and the preconditioner matrix B on a DMPlex.

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

        with self.hl as hl, self.lHbed as zb:
            dh = hl - zb
            dh[dh < 0.1] = 0.0
            Cd = self.minDiff + np.multiply(self.Cd, (1.0 - np.exp(-dh/self.H0)))

            # Coefficient derivatives
            Cp = np.multiply(self.Cd, np.exp(-dh/self.H0)/self.H0)
            nlC = jacobiancoeff(hl, Cd, Cp)

            for row in range(self.lpoints):
                B.setValuesLocal(row, row, a + nlC[row, 0])
                cols = self.FVmesh_ngbID[row, :]
                B.setValuesLocal(row, cols, nlC[row, 1:])
            B.assemble()

            if A != B:
                A.assemble()

        return True

    def _evalSolutionSoil(self, t, x):
        """
        Evaluate the initial solution of the SNES system.
        """

        assert t == 0.0, "only for t=0.0"
        x.setArray(self.h.getArray())

        return

    def diffuseSoil(self):
        r"""
        For river-transported sediments reaching the marine realm, this function computes the related marine deposition diffusion. It is based on a non-linear diffusion approach.

        .. math::
          \frac{\partial h}{\partial t}= \nabla \cdot \left( C_d \times (1.0 - exp^{h_s/H_0} \nabla h \right)

        It calls the following *private functions*:

        - _evalFunctionSoil
        - _evalJacobianSoil
        - _evalSolutionSoil

        .. note::

            PETSc SNES and time stepping TS approaches are used to solve the non-linear equation above over the considered time step.
        """

        t0 = process_time()

        x = self.tmp1.duplicate()
        f = self.tmp1.duplicate()

        # Get diffusion soil coefficient
        self.Cd = np.full(self.lpoints, self.Cda, dtype=np.float64)
        self.Cd[self.seaID] = self.Cdm

        # Remove the soil thickness from the elevation
        self.hLocal.copy(result=self.hl)
        self.dm.localToGlobal(self.hl, self.h)
        self.gHbed.waxpy(-1.0, self.Gsoil, self.hGlobal) 
        self.lHbed.waxpy(-1.0, self.Lsoil, self.hLocal) 

        # Time stepping definition
        ts = petsc4py.PETSc.TS().create(comm=petsc4py.PETSc.COMM_WORLD)
        # arkimex: IMEX Runge-Kutta schemes | rosw: Rosenbrock W-schemes
        ts.setType("rosw")

        ts.setIFunction(self._evalFunctionSoil, self.tmp1)
        ts.setIJacobian(self._evalJacobianSoil, self.mat)

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
        self._evalSolutionSoil(tstart, x)

        # Solve nonlinear equation
        ts.solve(x)

        if MPIrank == 0 and self.verbose:
            print(
                "Nonlinear soil diffusion solution (%0.02f seconds)" % (process_time() - t0),
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
        self.dh.waxpy(-1.0, self.hGlobal, x)
        self.dm.globalToLocal(self.dh, self.tmpL)
        chgSoil = self.tmpL.getArray().copy()
        self.tmpL.setArray(chgSoil)
        self.dm.localToGlobal(self.tmpL, self.tmp)

        x.destroy()
        f.destroy()
        if self.memclear:
            del depSoil
            gc.collect()
        petsc4py.PETSc.garbage_cleanup()

        # Update cumulative erosion and deposition as well as elevation
        self.cumED.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        # Update soil thickness
        self.updateSoilThickness()

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
                "Diffuse Soil Sediments (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        return