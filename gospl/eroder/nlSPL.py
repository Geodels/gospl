import os
import gc
import sys
import petsc4py
import numpy as np
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import local_spl
    from gospl._fortran import local_spl_coeff

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD


class nlSPL(object):
    """
    The class computes river incision expressed using a **stream power formulation** function of river discharge and slope with non-linear slope dependency.
    If the user has turned-on the sedimentation capability, this class also solves implicitly the **stream power formulation** accounting for a sediment transport/deposition term (`Yuan et al, 2019 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JF004867>`_).
    """

    def __init__(self, *args, **kwargs):
        """
        Initialisation of `nlSPL` class.
        """

        self.snes_rtol = 1.0e-6
        self.snes_atol = 1.e-6
        self.snes_maxit = 500

        return

    def _form_residual(self, snes, h, F):
        """
        The nonlinear system (SNES) at each time step is solved iteratively by assessing the residual of the detachment-limited SPL equation.

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
        S[S < 0] = 0.

        # Residuals based on the equation: h(t+dt) - h(t) + dt * K * A^m * S^n = 0
        res = h_array - self.hOldArray + self.Kbr * S**self.spl_n
        if self.iceOn:
            res += self.Kbi * S**self.spl_n

        # Residual vector
        F.setArray(res[self.glIDs])

        return

    def _form_residual_ed(self, snes, h, F):
        """
        The nonlinear system (SNES) at each time step is solved iteratively by assessing the residual of the SPL equation accounting for erosion and deposition (transport-limited).

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

        # Residuals based on the equation
        # h(t+dt) (1-G) - h(t) (1-G) + dt * K * A^m * S^n - dt * G * Qt / Area = 0
        res = (h_array - self.hOldArray) * (1.0 - self.fDep)
        res += self.Kbr * S**self.spl_n
        if self.iceOn:
            res += self.Kbi * S**self.spl_n
        res -= self.fDep * self.dt * Qt / self.larea

        # Residual vector
        F.setArray(res[self.glIDs])

        return

    def _monitor(self, snes, its, norm):
        """
        Non-linear SPL solver convergence evaluation.
        """

        if MPIrank == 0 and its % 50 == 0:
            print(f"  ---  Non-linear SPL solver iteration {its}, Residual norm: {norm}", flush=True)

    def _form_jacobian(self, snes, h, J, P):
        """
        For the detachment-limited case, the Jacobian is calculated and provided to the nonlinear system (SNES).

        Parameters:
        -----------
        snes : PETSc.SNES: The snes object.
        h : PETSc.Vec: The current solution vector (h^{n+1}) at the new time step.
        J : PETSc.Mat: The Jacobian matrix to be filled.
        P : PETSc.Mat: The preconditioner matrix to be filled.
        """

        # Current state
        self.dm.globalToLocal(h, self.hl)
        h_array = self.hl.getArray()

        # Compute slope and coefficients
        S, coeffs = local_spl_coeff(self.flowDir, h_array, self.rcvIDi,
                                    self.distRcvi, self.wghtVali)
        derivatives = np.zeros(self.lpoints)
        S[S < 0] = 0.
        coeffs[S == 0, :] = 0.
        ids = np.where(S > 0)[0]
        derivatives[ids] = self.spl_n * self.Kbr * np.power(S[ids], self.spl_n - 1.0)
        if self.iceOn:
            derivatives[ids] += self.spl_n * self.Kbi * np.power(S[ids], self.spl_n - 1.0)

        # Build Jacobian matrix
        for row in range(self.lpoints):
            P.setValuesLocal(row, row, 1.0 + derivatives[row] * coeffs[row, 0])
            cols = self.rcvIDi[row, :]
            P.setValuesLocal(row, cols, derivatives[row] * coeffs[row, 1:])
        P.assemble()

        if J != P:
            J.assemble()

        return True

    def _solveNL_ed(self):
        """
        Solves the non-linear stream power law for the transport limited case. This calls the following *private function*:

        - _form_residual_ed

        .. note::

            PETSc SNES approach is used to solve the nonlinear equation above over the considered time step. In this implementation of the SNES, we do not form the Jacobian and PETSc calculates it based on the residual function. A Nonlinear Generalized Minimum Residual method is used ``ngmres``, a Preconditioned Conjugate Gradient ``cg`` method is defined for the KSP and the preconditioner allows for multi-grid methods based on the HYPRE BoomerAMG approach.
        """

        self.oldH = self.hGlobal.getArray()
        self.hOldArray = self.hLocal.getArray().copy()

        # Upstream-averaged mean annual precipitation rate based on drainage area
        PA = self.FAL.getArray()

        # Define erosion limiter to prevent formation of flat
        dh = []
        for k in range(0, self.flowDir):
            dh.append(self.hOldArray - self.hOldArray[self.rcvIDi[:, k]])
        dh = np.array(dh).max(0)
        elimiter = np.divide(dh, dh + 1.0e-2, out=np.zeros_like(dh),
                             where=dh != 0)

        # Incorporate the effect of local mean annual precipitation rate on erodibility
        if self.sedfacVal is not None:
            self.Kbr = self.K * self.sedfacVal * (self.rainVal ** self.coeffd)
        else:
            self.Kbr = self.K * (self.rainVal ** self.coeffd)
        self.Kbr *= self.dt * (PA ** self.spl_m) * elimiter
        self.Kbr[self.seaID] = 0.0

        # In case glacial erosion is accounted for
        if self.iceOn:
            Ai = self.iceFAL.getArray()
            self.Kbi = self.Kice * self.dt * (Ai ** self.spl_m) * elimiter
            self.Kbi[self.seaID] = 0.0

        # Dimensionless depositional coefficient
        self.fDep = np.divide(self.fDepa * self.larea, PA, out=np.zeros_like(PA), where=PA != 0)
        self.fDep[self.seaID] = 0.
        self.fDep[self.fDep > 0.99] = 0.99
        if self.flatModel:
            self.fDep[self.idBorders] = 0.

        snes = petsc4py.PETSc.SNES().create(comm=petsc4py.PETSc.COMM_WORLD)
        snes.setTolerances(rtol=self.snes_rtol, atol=self.snes_atol,
                           max_it=self.snes_maxit)

        # Set a monitor to see residual values
        if self.verbose:
            snes.setMonitor(self._monitor)

        f = self.hGlobal.duplicate()
        snes.setFunction(self._form_residual_ed, f)

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

        petsc4py.PETSc.garbage_cleanup()

        return

    def _solveNL(self):
        """
        Solves the non-linear stream power law for the detachment limited case. This calls the following *private functions*:

        - _form_residual
        - _form_jacobian

        .. note::

            PETSc SNES approach is used to solve the nonlinear equation above over the considered time step. In this implementation of the SNES, we provide the Jacobian. A Nonlinear Generalized Minimum Residual method is used ``nrichardson``, a Generalized Minimum Residual method is used ``gmres`` for the KSP and the preconditioner allows for multi-grid methods based on the HYPRE BoomerAMG approach.
        """

        self.oldH = self.hGlobal.getArray()
        self.hOldArray = self.hLocal.getArray().copy()

        # Upstream-averaged mean annual precipitation rate based on drainage area
        PA = self.FAL.getArray()

        # Define erosion limiter to prevent formation of flat
        dh = []
        for k in range(0, self.flowDir):
            dh.append(self.hOldArray - self.hOldArray[self.rcvIDi[:, k]])
        dh = np.array(dh).max(0)
        elimiter = np.divide(dh, dh + 1.0e-2, out=np.zeros_like(dh),
                             where=dh != 0)

        # Incorporate the effect of local mean annual precipitation rate on erodibility
        if self.sedfacVal is not None:
            self.Kbr = self.K * self.sedfacVal * (self.rainVal ** self.coeffd)
        else:
            self.Kbr = self.K * (self.rainVal ** self.coeffd)
        self.Kbr *= self.dt * (PA ** self.spl_m) * elimiter
        self.Kbr[self.seaID] = 0.0

        # In case glacial erosion is accounted for
        if self.iceOn:
            Ai = self.iceFAL.getArray()
            self.Kbi = self.Kice * self.dt * (Ai ** self.spl_m) * elimiter
            self.Kbi[self.seaID] = 0.0

        snes = petsc4py.PETSc.SNES().create(comm=petsc4py.PETSc.COMM_WORLD)
        snes.setTolerances(rtol=self.snes_rtol, atol=self.snes_atol,
                           max_it=self.snes_maxit)

        # Set a monitor to see residual values
        if self.verbose:
            snes.setMonitor(self._monitor)

        f = self.hGlobal.duplicate()
        snes.setFunction(self._form_residual, f)

        # Setting the Jacobian function
        J = self.dm.createMatrix()
        J.setOption(J.Option.NEW_NONZERO_LOCATIONS, True)
        snes.setJacobian(self._form_jacobian, J)

        # SNES solvers
        snes.setType('nrichardson')
        ksp = snes.getKSP()
        ksp.setType("gmres")
        pc = ksp.getPC()
        pc.setType("hypre")  # hypre or gamg
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
        J.destroy()

        petsc4py.PETSc.garbage_cleanup()

        return

    def _getEroDepRateNL(self):
        r"""
        This function computes erosion deposition rates in metres per year. This is done on the filled elevation. We use the filled-limited elevation to ensure that erosion/deposition is not going to be underestimated by small depressions which are likely to be filled (either by sediments or water) during a single time step.

        The simplest law to simulate fluvial incision is based on the detachment-limited stream power law, in which erosion rate  depends on drainage area :math:`A`, net precipitation :math:`P` and local slope :math:`S` and takes the form:

        .. math::

          E = − \kappa P^d (PA)^m S^n

        :math:`\kappa` is a dimensional coefficient describing the erodibility of the channel bed as a function of rock strength, bed roughness and climate, :math:`d`, :math:`m` and :math:`n` are dimensionless positive constants.

        A similar approach is used to compute ice induced erosion where the ice flow accumulation is defined based on downstream nodes and is smoothed to better represent the erosion induced by glaciers. The ice-induced erosion uses the stream power law equation with a eordibility coefficient which is user defined. Under glacier terminus point, melted glacier flow is added to river flow accumulation likewise is the glacier-induced transported sediment flux.

        Default formulation assumes :math:`d = 0`, :math:`m = 0.5` and :math:`n = 1`. The precipitation exponent :math:`d` allows for representation of climate-dependent chemical weathering of river bed across non-uniform rainfall.

        .. important::

            Here, the coefficient `n` can be fixed by the user to value different than 1.0 and the equation is also dependent on `m`, `d` and the erodibility :math:`\kappa`.

        In addition, an alternative method to the purely detachment-limited approach consists in accounting for the role played by sediment in modulating erosion and deposition rates. It follows the model of `Yuan et al, 2019 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JF004867>`_, whereby the deposition flux depends on a deposition coefficient :math:`G` and is proportional to the ratio between cell area :math:`\mathrm{\Omega}` and water discharge :math:`\mathrm{Q}=\bar{P}A`.
        """

        t0 = process_time()

        # Build the SPL erosion arrays
        if self.flexOn:
            self.hLocal.copy(result=self.hOldFlex)

        if self.fDepa == 0:
            self._solveNL()
        else:
            self._solveNL_ed()

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

    def erodepSPLnl(self):
        """
        Modified **stream power law** model used to represent erosion by rivers also taking into account the role played by sediments in modulating erosion and deposition rate and considering **non-linear slope dependency**.

        It calls the private function `_getEroDepRateNL` described above. Once erosion/deposition rates have been calculated, the function computes local thicknesses for the considered time step and update local elevation and cumulative erosion, deposition values.
        """

        t0 = process_time()

        # Computes the erosion deposition rates based on flow accumulation
        self.Eb.set(0.0)
        self.hGlobal.copy(result=self.hOld)
        self.dm.globalToLocal(self.hOld, self.hOldLocal)
        self._getEroDepRateNL()

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
