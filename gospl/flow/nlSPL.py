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


class NLSPL(object):
    """
    This class calculates
    """

    def __init__(self, *args, **kwargs):
        """
        The initialisation of `NLSPL` class consists in the declaration of PETSc vectors and matrices.
        """

        self.snes_rtol = 1.0e-6
        self.snes_atol = 1.e-6
        self.snes_max_it = 500

        return

    def form_residual(self, snes, h, F):

        # Current state
        self.dm.globalToLocal(h, self.hl)
        h_array = self.hl.getArray()

        # Compute slope
        S = local_spl(self.flowDir, h_array, self.rcvIDi,
                      self.distRcvi, self.wghtVali)
        S[S < 0] = 0.

        # Residuals based on the equation: h(t+dt) - h(t) + dt * K * A^m * S^n = 0
        res = h_array - self.hOldArray + self.Kbr * S**self.spl_n

        # Residual vector
        F.setArray(res[self.glIDs])

        return

    def monitor(self, snes, its, norm):

        if MPIrank == 0 and its % 10 == 0:
            print(f"  ---  Non-linear SPL solver iteration {its}, Residual norm: {norm}", flush=True)

    def form_jacobian(self, snes, h, J, P):

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

        # Build Jacobian matrix
        for row in range(self.lpoints):
            P.setValuesLocal(row, row, 1.0 + derivatives[row] * coeffs[row, 0])
            cols = self.rcvIDi[row, :]
            P.setValuesLocal(row, cols, derivatives[row] * coeffs[row, 1:])
        P.assemble()

        if J != P:
            J.assemble()

        return True

    def _solveNL(self):
        r"""
        It calls the following *private functions*:

        - _evalFunction
        - _evalJacobian

        .. note::

            PETSc SNES and time stepping TS approaches are used to solve the nonlinear equation above over the considered time step.

        :return: ndepo (updated deposition numpy arrays)
        """

        self.oldH = self.hGlobal.getArray()
        self.hOldArray = self.hLocal.getArray().copy()

        # Upstream-averaged mean annual precipitation rate based on drainage area
        PA = self.FAL.getArray()

        # Incorporate the effect of local mean annual precipitation rate on erodibility
        if self.sedfacVal is not None:
            self.Kbr = self.K * self.sedfacVal * (self.rainVal ** self.coeffd)
        else:
            self.Kbr = self.K * (self.rainVal ** self.coeffd)
        self.Kbr *= self.dt * (PA ** self.spl_m)
        self.Kbr[self.seaID] = 0.0

        snes = petsc4py.PETSc.SNES().create(comm=petsc4py.PETSc.COMM_WORLD)

        snes.setTolerances(rtol=self.snes_rtol, atol=self.snes_atol,
                           max_it=self.snes_max_it)

        # Set a monitor to see residual values
        if self.verbose:
            snes.setMonitor(self.monitor)

        f = self.hGlobal.duplicate()
        snes.setFunction(self.form_residual, f)

        # Setting the Jacobian function
        J = self.dm.createMatrix()
        J.setOption(J.Option.NEW_NONZERO_LOCATIONS, True)
        snes.setJacobian(self.form_jacobian, J)

        # SNES solvers
        snes.setType('nrichardson')
        ksp = snes.getKSP()
        ksp.setType("gmres")
        pc = ksp.getPC()
        pc.setType("hypre")  # hypre or gamg
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
        """

        t0 = process_time()

        # Build the SPL erosion arrays
        if self.flexOn:
            self.hLocal.copy(result=self.hOldFlex)

        self._solveNL()

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

    def erodepSPLNL(self):
        """
        Modified **stream power law** model used to represent erosion by rivers also taking into account the role played by sediment in modulating erosion and deposition rate.

        It calls the private function `_getEroDepRate` described above. Once erosion/deposition rates have been calculated, the function computes local thicknesses for the considered time step and update local elevation and cumulative erosion, deposition values.

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
