import os
import gc
import sys
import petsc4py
import numpy as np
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import setmaxnb
    from gospl._fortran import sethillslopecoeff
    from gospl._fortran import scale_volume

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD


class SEDMesh(object):
    """
    This class encapsulates all the functions related to sediment transport, production and deposition in the continental domain. `gospl` has the ability to track three types of clastic sediment size and one type of carbonate (still under development). The following processes are considered:

    - inland river deposition in depressions and enclosed sea
    - hillslope processes in both marine and inland areas

    .. note::
        All these functions are ran in parallel using the underlying PETSc library.

    """

    def __init__(self, *args, **kwargs):
        """
        The initialisation of `SEDMesh` class consists in the declaration of several PETSc vectors.
        """

        # Petsc vectors
        self.tmp = self.hGlobal.duplicate()
        self.tmpL = self.hLocal.duplicate()
        self.tmp1 = self.hGlobal.duplicate()
        self.Qs = self.hGlobal.duplicate()
        self.QsL = self.hLocal.duplicate()

        self.vSed = self.hGlobal.duplicate()
        self.vSedLocal = self.hLocal.duplicate()
        if self.stratNb > 0:
            if self.stratF is not None:
                self.vSedf = self.hGlobal.duplicate()
                self.vSedfLocal = self.hLocal.duplicate()
            if self.stratW is not None:
                self.vSedw = self.hGlobal.duplicate()
                self.vSedwLocal = self.hLocal.duplicate()
            if self.carbOn:
                self.vSedc = self.hGlobal.duplicate()
                self.vSedcLocal = self.hLocal.duplicate()

        # Get the maximum number of neighbours on the mesh
        maxnb = np.zeros(1, dtype=np.int64)
        maxnb[0] = setmaxnb(self.lpoints)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, maxnb, op=MPI.MAX)
        self.maxnb = maxnb[0]

        return

    def getSedFlux(self):
        """
        This function computes sediment flux in cubic metres per year from incoming rivers. Like for the computation of the flow discharge and erosion rates, the sediment flux is solved by an implicit time integration method, the matrix system is the one obtained from the receiver distributions over the unfilled elevation mesh for the flow discharge (`fMat`).

        The PETSc *scalable linear equations solvers* (**KSP**) is used here again with an iterative method obtained from PETSc Richardson solver (`richardson`) with block Jacobian preconditioning (`bjacobi`).
        """
        t0 = process_time()

        # Multi-lithology case
        if self.stratNb > 0:
            # Coarse sediment
            # Get erosion rate (m/yr) to volume
            self.tmpL.setArray(self.thCoarse)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.tmp.pointwiseMult(self.tmp, self.areaGlobal)
            # Get the volume of sediment transported in m3 per year
            self._solve_KSP(False, self.fMati, self.tmp, self.vSed)
            # Update local vector
            self.dm.globalToLocal(self.vSed, self.vSedLocal)

            # Fine sediment
            if self.stratF is not None:
                # Get erosion rate (m/yr) to volume
                self.tmpL.setArray(self.thFine)
                self.dm.localToGlobal(self.tmpL, self.tmp)
                self.tmp.pointwiseMult(self.tmp, self.areaGlobal)
                # Get the volume of sediment transported in m3 per year
                self._solve_KSP(False, self.fMati, self.tmp, self.vSedf)
                # Update local vector
                self.dm.globalToLocal(self.vSedf, self.vSedfLocal)

            # Weathered sediment
            if self.stratW is not None:
                # Get erosion rate (m/yr) to volume
                self.tmpL.setArray(self.thClay)
                self.dm.localToGlobal(self.tmpL, self.tmp)
                self.tmp.pointwiseMult(self.tmp, self.areaGlobal)
                # Get the volume of sediment transported in m3 per year
                self._solve_KSP(False, self.fMati, self.tmp, self.vSedw)
                # Update local vector
                self.dm.globalToLocal(self.vSedw, self.vSedwLocal)

            # Carbonate sediment
            if self.carbOn:
                # Get erosion rate (m/yr) to volume
                self.tmpL.setArray(self.thCarb)
                self.dm.localToGlobal(self.tmpL, self.tmp)
                self.tmp.pointwiseMult(self.tmp, self.areaGlobal)
                # Get the volume of sediment transported in m3 per year
                self._solve_KSP(False, self.fMati, self.tmp, self.vSedc)
                # Update local vector
                self.dm.globalToLocal(self.vSedc, self.vSedcLocal)

        else:
            # Get erosion rate (m/yr) to volume
            self.Eb.copy(result=self.tmp)
            self.tmp.pointwiseMult(self.tmp, self.areaGlobal)
            # Get the volume of sediment transported in m3 per year
            self._solve_KSP(False, self.fMati, self.tmp, self.vSed)
            # Update local vector
            self.dm.globalToLocal(self.vSed, self.vSedLocal)

        if MPIrank == 0 and self.verbose:
            print(
                "Update Sediment Load (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def _getSedVol(self, stype):
        """
        Pick the relevant PETSc array for the specified sediment type.
        """

        if stype == 0:
            self.vSedLocal.copy(result=self.QsL)
        elif stype == 1:
            self.vSedfLocal.copy(result=self.QsL)
        elif stype == 2:
            self.vSedcLocal.copy(result=self.QsL)
        elif stype == 3:
            self.vSedwLocal.copy(result=self.QsL)
        self.dm.localToGlobal(self.QsL, self.Qs)

        return

    def _moveDownstream(self, vSed, step):
        """
        In cases where river sediment fluxes drain into depressions, they might fill the sink completely and overspill or be deposited in it. This function computes the excess of sediment (if any) able to flow dowstream.

        .. important::

            The excess sediment volume is then added to the downstream sediment flux (`vSed`).

        :arg vSed: excess sediment volume array
        :arg step: downstream distribution step

        :return: excess (boolean set to True is excess sediment fluxes remain to be distributed)
        """
        excess = False

        # Remove points belonging to other processors
        vSed = np.multiply(vSed, self.inIDs)

        # Get volume incoming in each depression
        grp = npi.group_by(self.pitIDs[self.lsink])
        uID = grp.unique
        _, vol = grp.sum(vSed[self.lsink])
        inV = np.zeros(len(self.pitParams), dtype=np.float64)
        ids = uID > -1
        inV[uID[ids]] = vol[ids]

        # Combine incoming volume globally
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, inV, op=MPI.SUM)

        # Get excess volume to distribute downstream
        eV = inV - self.pitVol
        if (eV > 0.0).any():
            eIDs = eV > 0.0
            self.pitVol[eIDs] = 0.0
            spillIDs = self.pitInfo[eIDs, 0]
            localSpill = np.where(self.pitInfo[eIDs, 1] == MPIrank)[0]
            localPts = spillIDs[localSpill]
            nvSed = np.zeros(self.lpoints, dtype=np.float64)
            nvSed[localPts] = eV[eIDs][localSpill]
            ids = np.in1d(self.pitIDs, np.where(eV > 0.0)[0])
            self.sedFilled[ids] = self.lFill[ids]

        # Update unfilled depressions volumes
        eIDs = eV < 0
        self.pitVol[eIDs] -= inV[eIDs]
        self.pitVol[self.pitVol < 0] = 0.0

        # In case there is still remaining water flux to distribute downstream
        if (eV > 1.0e-3).any():
            if step == 100:
                self._buildFlowDirection(self.lFill)
            else:
                self._buildFlowDirection(self.sedFilled)
            self.tmpL.setArray(nvSed)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            if self.tmp.sum() > 0.5 * self.maxarea[0]:
                excess = True
                self._solve_KSP(True, self.fMat, self.tmp, self.tmp1)
                self.dm.globalToLocal(self.tmp1, self.tmpL)

        return excess

    def _distributeSediment(self, hl, stype=0):
        """
        This function finds the continental sediment volumes to distribute downstream (down to the ocean or sinks) for a given iteration.

        """

        t0 = process_time()

        # Define the volume to distribute
        self._getSedVol(stype)

        # Get the volumetric sediment rate (m3/yr) to distribute during the time step and convert it in volume (m3)
        vSed = self.QsL.getArray().copy() * self.dt

        step = 0
        excess = True
        while excess:
            t1 = process_time()
            excess = self._moveDownstream(vSed, step)
            if excess:
                vSed = self.tmpL.getArray().copy()
                nvSed = vSed.copy()
                nvSed[hl < self.sedFilled] = 0.0
                self.tmpL.setArray(nvSed / self.dt)
                vSed[self.idBorders] = 0.0
                if stype == 0:
                    self.vSedLocal.axpy(1.0, self.tmpL)
                if stype == 1:
                    self.vSedfLocal.axpy(1.0, self.tmpL)
                if stype == 3:
                    self.vSedwLocal.axpy(1.0, self.tmpL)
                if stype == 2:
                    self.vSedcLocal.axpy(1.0, self.tmpL)
            if MPIrank == 0 and self.verbose:
                print(
                    "Downstream sediment flux computation step %d (%0.02f seconds)"
                    % (step, process_time() - t1),
                    flush=True,
                )
            step += 1
        if stype == 0:
            self.dm.localToGlobal(self.vSedLocal, self.vSed)
        if stype == 1:
            self.dm.localToGlobal(self.vSedfLocal, self.vSedf)
        if stype == 3:
            self.dm.localToGlobal(self.vSedwLocal, self.vSedw)
        if stype == 2:
            self.dm.localToGlobal(self.vSedcLocal, self.vSedc)

        if MPIrank == 0 and self.verbose:
            print(
                "Distribute Continental Sediments (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        return

    def _updateSinks(self, hl, stype=0):
        """
        This function updates depressions elevations based on incoming sediment volumes. The function runs until all inland sediments have reached either a sink which is not completely filled or the open ocean.

        """

        # Difference between intial volume and remaining one
        depo = self.pitParams[:, 0] - self.pitVol
        depo[depo < 0] = 0.0
        with np.errstate(divide="ignore", over="ignore"):
            scaleV = np.divide(
                depo,
                self.pitParams[:, 0],
                out=np.zeros_like(self.pitParams[:, 0]),
                where=self.pitParams[:, 0] != 0,
            )
        scaleV[scaleV > 1.0] = 1.0
        scale = scale_volume(self.pitIDs, scaleV)

        # Update cumulative erosion and deposition as well as elevation
        self.tmpL.setArray((self.lFill - hl) * scale)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.cumED.axpy(1.0, self.tmp)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        # Update stratigraphic layer parameters
        if self.stratNb > 0:
            self.deposeStrat(stype)
            self.elevStrat()

        # In case there is other sediment type
        self.pitParams[:, 0] = depo.copy()

        return

    def sedChange(self):
        """
        This function is the main entry point to perform both continental and marine river-induced deposition. It calls the private function:

        - _distributeSediment
        - _updateSinks

        """

        # Find Continental Sediment Fluxes
        self.getSedFlux()

        # Compute depressions information as elevations changed due to erosion
        self.fillElevation(sed=True)
        hl = self.hLocal.getArray().copy()
        self.lsink = self.lsinki.copy()
        self.pitVol = self.pitParams[:, 0].copy()

        # Distribute inland sediments
        self.sedFilled = hl.copy()
        if self.stratNb > 0:
            self._distributeSediment(hl, stype=0)
            self._updateSinks(hl, stype=0)
            hl = self.hLocal.getArray().copy()
            if self.carbOn:
                self._distributeSediment(hl, stype=2)
                self._updateSinks(hl, stype=2)
                hl = self.hLocal.getArray().copy()
            if self.stratF is not None:
                self._distributeSediment(hl, stype=1)
                self._updateSinks(hl, stype=1)
                hl = self.hLocal.getArray().copy()
            if self.stratW is not None:
                self._distributeSediment(hl, stype=3)
                self._updateSinks(hl, stype=3)
        else:
            self._distributeSediment(hl, stype=0)
            self._updateSinks(hl, stype=0)

        return

    def getHillslope(self):
        r"""
        This function computes hillslope processes. The code assumes that gravity is the main driver for transport and states that the flux of sediment is proportional to the gradient of topography.

        As a result, we use a linear diffusion law commonly referred to as **soil creep**:

        .. math::
          \frac{\partial z}{\partial t}= \kappa_{D} \nabla^2 z

        in which :math:`\kappa_{D}` is the diffusion coefficient and can be defined with different values for the marine and land environments (set with `hillslopeKa` and `hillslopeKm` in the YAML input file). It encapsulates, in a simple formulation, processes operating on superficial sedimentary layers. Main controls on variations of :math:`\kappa_{D}` include substrate, lithology, soil depth, climate and biological activity.
        """

        # Compute Hillslope Diffusion Law
        h = self.hLocal.getArray().copy()
        self.seaID = np.where(h <= self.sealevel)[0]
        self._hillSlope()

        # Update layer elevation
        if self.stratNb > 0:
            self.elevStrat()
        if self.memclear:
            del h
            gc.collect()

        return

    def _hillSlope(self, smooth=0):
        r"""
        This function computes hillslope using a linear diffusion law commonly referred to as **soil creep**:

        .. math::
          \frac{\partial z}{\partial t}= \kappa_{D} \nabla^2 z

        in which :math:`\kappa_{D}` is the diffusion coefficient and can be defined with different values for the marine and land environments (set with `hillslopeKa` and `hillslopeKm` in the YAML input file).

        .. note::
            The hillslope processes in `gospl` are considered to be happening at the same rate for coarse and fine sediment sizes.

        :arg smooth: integer specifying if the diffusion equation is used for smoothing (1) or for marine deposits (2).
        """

        if not smooth:
            if self.Cda == 0.0 and self.Cdm == 0.0:
                return

        t0 = process_time()

        # Diffusion matrix construction
        if smooth == 1:
            Cd = np.full(self.lpoints, self.smthK, dtype=np.float64)
        elif smooth == 2:
            Cd = np.full(self.lpoints, self.smthD, dtype=np.float64)
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
            diffMat += tmpMat
            tmpMat.destroy()

        # Get elevation values for considered time step
        if smooth == 1:
            self._solve_KSP(True, diffMat, self.hGlobal, self.tmp)
            diffMat.destroy()
            self.dm.globalToLocal(self.tmp, self.tmpL)
            return self.tmpL.getArray().copy()
        elif smooth == 2:
            self._solve_KSP(True, diffMat, self.tmp1, self.tmp)
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
                if self.stratW is not None:
                    self.deposeStrat(3)
                else:
                    self.deposeStrat(0)

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Hillslope Processes (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return
