import os
import gc
import sys
import petsc4py
import numpy as np
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import mfdreceivers

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD


class SPL(object):
    """
    The class computes river incision expressed using a **stream power formulation** function of river discharge and slope.

    .. note::
        This class assumes that the stream power law is defined such that it varies linearly with slope (*i.e.*, ``n`` exponent set to 1.). In such case, a linear expression is used and is most simpler to compute than for the non-linear case.

    If the user has turned-on the sedimentation capability, this class solves implicitly the **stream power formulation** accounting for a sediment transport/deposition term (`Yuan et al, 2019 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JF004867>`_).
    """

    def __init__(self, *args, **kwargs):
        """
        The initialisation of `SPL` class consists in the declaration of PETSc vectors.
        """

        # Petsc vectors
        self.hOld = self.hGlobal.duplicate()
        self.hOldLocal = self.hLocal.duplicate()
        self.hOldFlex = self.hLocal.duplicate()
        self.Eb = self.hGlobal.duplicate()
        self.stepED = self.hGlobal.duplicate()
        self.EbLocal = self.hLocal.duplicate()
        self.newH = self.hGlobal.duplicate()
        self.EbLocal.set(0.)

        return

    def _eroMats(self, hOldArray):
        """
        Builds the erosion matrices used to solve implicitly the stream power equations for the river and ice processes.

        :arg hOldArray: local elevation array from previous time step

        :return: eMat, PA where the first is a sparse PETSc matrices related to river and glacial erosion and PA is the flow/ice accumulation rate.
        """

        # Upstream-averaged mean annual precipitation rate based on drainage area
        PA = self.FAL.getArray()

        # Incorporate the effect of local mean annual precipitation rate on erodibility
        if self.sedfacVal is not None:
            Kbr = self.K * self.sedfacVal * (self.rainVal ** self.coeffd)
        else:
            Kbr = self.K * (self.rainVal ** self.coeffd)
        Kbr *= self.dt * (PA ** self.spl_m)
        Kbr[self.seaID] = 0.0

        # In case glacial erosion is accounted for
        if self.iceOn:
            GA = self.iceFAL.getArray()
            Kbi = self.dt * self.Kice * (GA ** self.ice_m)
            PA += GA

        # Initialise matrices...
        eMat = self.iMat.copy()
        wght = self.wghtVali.copy()
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        nodes = indptr[:-1]
        wght[self.seaID, :] = 0.0

        # Define erosion coefficients
        for k in range(0, self.flowDir):

            # Define erosion limiter to prevent formation of flat
            dh = hOldArray - hOldArray[self.rcvIDi[:, k]]
            limiter = np.divide(dh, dh + 1.0e-2, out=np.zeros_like(dh), where=dh != 0)

            # Bedrock erosion processes SPL computation (maximum bedrock incision)
            data = np.divide(
                Kbr * limiter,
                self.distRcvi[:, k],
                out=np.zeros_like(PA),
                where=self.distRcvi[:, k] != 0,
            )
            tmpMat = self._matrix_build()
            data = np.multiply(data, -wght[:, k])
            data[self.rcvIDi[:, k].astype(petsc4py.PETSc.IntType) == nodes] = 0.0
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                self.rcvIDi[:, k].astype(petsc4py.PETSc.IntType),
                data,
            )
            tmpMat.assemblyEnd()
            eMat.axpy(1.0, tmpMat)
            tmpMat.destroy()
            tmpMat = self._matrix_build_diag(data)
            eMat.axpy(-1.0, tmpMat)
            tmpMat.destroy()

            if self.iceOn:
                data = np.divide(
                    Kbi * limiter,
                    self.distRcvi[:, k],
                    out=np.zeros_like(GA),
                    where=self.distRcvi[:, k] != 0,
                )
                tmpMat = self._matrix_build()
                data = np.multiply(data, -wght[:, k])
                data[self.rcvIDi[:, k].astype(petsc4py.PETSc.IntType) == nodes] = 0.0
                tmpMat.assemblyBegin()
                tmpMat.setValuesLocalCSR(
                    indptr,
                    self.rcvIDi[:, k].astype(petsc4py.PETSc.IntType),
                    data,
                )
                tmpMat.assemblyEnd()
                eMat.axpy(1.0, tmpMat)
                tmpMat.destroy()
                tmpMat = self._matrix_build_diag(data)
                eMat.axpy(-1.0, tmpMat)
                tmpMat.destroy()

        if self.memclear:
            del dh, limiter, wght, data
            gc.collect()

        return eMat, PA

    def _coupledEDSystem(self, eMat):
        r"""
        Setup matrix for the coupled linear system in which the SPL model takes into account sediment deposition.

        .. note::

            The approach follows `Yuan et al, 2019 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JF004867>`_, where the deposition flux depends on a deposition coefficient :math:`G` and is proportional to the ratio between cell area :math:`A` and flow accumulation :math:`FA`.

        The approach considers the local balance between erosion and deposition and is based on sediment flux resulting from net upstream erosion.

        .. math::

            \mathrm{\frac{\eta_i^{t+\Delta t}-\eta_i^t}{\Delta t}} =  \mathrm{-\kappa P^d_i \sqrt{Q_i} \frac{\eta_i^{t+\Delta t} - \eta_{rcv}^{t+\Delta t}}{\lambda_{i,rcv}}} + \mathrm{G' Q_{s_i} / \Omega_i}

        where :math:`\mathrm{\lambda_{i,rcv}}` is the length of the edges connecting the considered vertex to its receiver and :math:`\mathrm{\Omega_i}` is the area (voronoi) of the node :math:`i`.

        :math:`\mathrm{Q_{s_i}}` is the upstream incoming sediment flux in m3/yr and :math:`\mathrm{G'}` is equal to :math:`\mathrm{G \Omega_i / \bar{P}A}`.

        The upstream incoming sediment flux is obtained from the total sediment flux :math:`\mathrm{Q_{t_i}}` where:

        .. math::

            \mathrm{Q_{t_i}^{t+\Delta t} - \sum_{ups} w_{i,j} Q_{t_u}^{t+\Delta t}}= \mathrm{(\eta_i^{t} - \eta_i^{t+\Delta t}) \frac{\Delta t}{\Omega_i}}

        which gives:

        .. math::

            \mathrm{Q_{s_i}} = \mathrm{Q_{t_i}} - \mathrm{(\eta_i^{t} - \eta_i^{t+\Delta t}) \frac{\Delta t}{\Omega_i}}

        This system of coupled equations is solved implicitly using PETSc by assembling the matrix and vectors using the nested submatrix and subvectors and by using the ``fieldsplit`` preconditioner combining two separate preconditioners for the collections of variables.

        :arg eMat: erosion matrix (from the purely-erosive SPL model)
        """

        # Define submatrices
        A00 = self._matrix_build_diag(-self.fDep)
        A01 = self._matrix_build_diag(-self.fDep * self.dt / self.larea)
        A10 = self._matrix_build_diag(self.larea / self.dt)

        # Assemble the matrix for the coupled system
        A00.axpy(1.0, eMat)
        mats = [[A00, A01], [A10, self.fMati]]
        sysMat = petsc4py.PETSc.Mat().createNest(mats=mats, comm=MPIcomm)
        sysMat.assemblyBegin()
        sysMat.assemblyEnd()

        # Clean up
        A00.destroy()
        A01.destroy()
        A10.destroy()
        eMat.destroy()
        mats[0][0].destroy()
        mats[0][1].destroy()
        mats[1][0].destroy()
        mats[1][1].destroy()

        # Create nested vectors
        self.tmpL.setArray(1. - self.fDep)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.tmp.pointwiseMult(self.tmp, self.hOld)
        self.tmp1.pointwiseMult(self.hOld, self.areaGlobal)
        self.tmp1.scale(1. / self.dt)
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
        subksps[0].getPC().setType("gasm")
        subksps[1].setType("preonly")
        subksps[1].getPC().setType("gasm")

        ksp.solve(rhs_vec, hq_vec)
        r = ksp.getConvergedReason()
        if r < 0:
            KSPReasons = self._make_reasons(petsc4py.PETSc.KSP.ConvergedReason())
            if MPIrank == 0:
                print(
                    "LinearSolver failed to converge after iterations",
                    ksp.getIterationNumber(),
                    flush=True,
                )
                print("with reason: ", KSPReasons[r], flush=True)
        else:
            if MPIrank == 0 and self.verbose:
                print(
                    "LinearSolver converge after %d iterations"
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

        return

    def _getEroDepRate(self):
        r"""
        This function computes erosion deposition rates in metres per year. This is done on the filled elevation. We use the filled-limited elevation to ensure that erosion/deposition is not going to be underestimated by small depressions which are likely to be filled (either by sediments or water) during a single time step.

        The simplest law to simulate fluvial incision is based on the detachment-limited stream power law, in which erosion rate  depends on drainage area :math:`A`, net precipitation :math:`P` and local slope :math:`S` and takes the form:

        .. math::

          E = − \kappa P^d (PA)^m S^n

        :math:`\kappa` is a dimensional coefficient describing the erodibility of the channel bed as a function of rock strength, bed roughness and climate, :math:`d`, :math:`m` and :math:`n` are dimensionless positive constants.

        A similar approach is used to compute ice induced erosion where the ice flow accumulation is defined based on downstream nodes and is smoothed to better represent the erosion induced by glaciers. The ice-induced erosion uses the stream power law equation with a eordibility coefficient which is user defined. Under glacier terminus point, melted glacier flow is added to river flow accumulation likewise is the glacier-induced transported sediment flux.

        Default formulation assumes :math:`d = 0`, :math:`m = 0.5` and :math:`n = 1`. The precipitation exponent :math:`d` allows for representation of climate-dependent chemical weathering of river bed across non-uniform rainfall.

        .. important::

            Here, the coefficient `n` is fixed and the equation is tuned based on `m`, `d` and the erodibility :math:`\kappa`.

        The erosion rate is solved by an implicit time integration method, the matrix system is based on the receiver distributions and is assembled from local Compressed Sparse Row (**CSR**) matrices into a global PETSc matrix. The PETSc *scalable linear equations solvers* (**KSP**) is used with both an iterative method and a preconditioner and erosion rate solution is obtained using PETSc Richardson solver (`richardson`) with block Jacobian preconditioning (`bjacobi`).

        In addition, an alternative method to the purely detachment-limited approach consists in accounting for the role played by sediment in modulating erosion and deposition rates. It follows the model of `Yuan et al, 2019 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JF004867>`_, whereby the deposition flux depends on a deposition coefficient :math:`G` and is proportional to the ratio between cell area :math:`\mathrm{\Omega}` and water discharge :math:`\mathrm{Q}=\bar{P}A`.
        """

        t0 = process_time()
        # Build the SPL erosion arrays
        hOldArray = self.hLocal.getArray().copy()
        if self.flexOn:
            self.hLocal.copy(result=self.hOldFlex)
        eMat, PA = self._eroMats(hOldArray)

        # Solve SPL erosion implicitly for fluvial and glacial erosion
        if self.fDepa == 0:
            t1 = process_time()
            self._solve_KSP(True, eMat, self.hOld, self.stepED)
            self.tmp.waxpy(-1.0, self.hOld, self.stepED)
            eMat.destroy()
            if MPIrank == 0 and self.verbose:
                print(
                    "Solve SPL erosion (%0.02f seconds)" % (process_time() - t1),
                    flush=True,
                )
        # Accounting for continental sediment deposition
        else:
            t1 = process_time()
            # Dimensionless depositional coefficient
            self.fDep = np.divide(self.fDepa * self.larea, PA, out=np.zeros_like(PA), where=PA != 0)
            self.fDep[self.seaID] = 0.
            self.fDep[self.fDep > 0.99] = 0.99
            if self.flatModel:
                self.fDep[self.idBorders] = 0.
            self._coupledEDSystem(eMat)
            eMat.destroy()
            if MPIrank == 0 and self.verbose:
                print(
                    "Solve SPL accounting for sediment deposition (%0.02f seconds)" % (process_time() - t1),
                    flush=True,
                )
            self.tmp.waxpy(-1.0, self.hOld, self.newH)

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

        if self.memclear:
            del PA, hOldArray, E
            gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "Finalise erosion deposition rates (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def erodepSPL(self):
        """
        Modified **stream power law** model used to represent erosion by rivers also taking into account the role played by sediment in modulating erosion and deposition rate.

        It calls the private function `_getEroDepRate` described above. Once erosion/deposition rates have been calculated, the function computes local thicknesses for the considered time step and update local elevation and cumulative erosion, deposition values.

        """

        t0 = process_time()

        # Computes the erosion deposition rates based on flow accumulation
        self.Eb.set(0.0)
        self.hGlobal.copy(result=self.hOld)
        self.dm.globalToLocal(self.hOld, self.hOldLocal)
        self._getEroDepRate()

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
