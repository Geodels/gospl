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
    from gospl._fortran import scale_volume

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD


class SEDMesh(object):
    """
    This class encapsulates all the functions related to sediment transport and deposition in the continental domain.
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
        self.nQs = self.hLocal.duplicate()

        self.vSed = self.hGlobal.duplicate()
        self.vSedLocal = self.hLocal.duplicate()

        # Get the maximum number of neighbours on the mesh
        maxnb = np.zeros(1, dtype=np.int64)
        maxnb[0] = setmaxnb(self.lpoints)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, maxnb, op=MPI.MAX)
        self.maxnb = maxnb[0]

        return

    def _getSedFlux(self):
        """
        This function computes sediment flux in cubic metres per year from incoming rivers. Like for the computation of the flow discharge and erosion rates, the sediment flux is solved by an implicit time integration method, the matrix system is the one obtained from the receiver distributions over the unfilled elevation mesh for the flow discharge (`fMat`).

        The PETSc *scalable linear equations solvers* (**KSP**) is used here again with an iterative method obtained from PETSc Richardson solver (`richardson`) with block Jacobian preconditioning (`bjacobi`).
        """
        t0 = process_time()

        # Stratigraphic layers exist
        if self.stratNb > 0:
            # Get erosion rate (m/yr) to volume
            self.tmpL.setArray(self.thCoarse)
            self.dm.localToGlobal(self.tmpL, self.tmp)
        else:
            # Get erosion rate (m/yr) to volume
            self.Eb.copy(result=self.tmp)

        # Get the volume of sediment transported in m3/yr
        self.tmp.pointwiseMult(self.tmp, self.areaGlobal)
        self._solve_KSP(False, self.fMati, self.tmp, self.vSed)
        self.fMati.destroy()

        # Update local vector
        self.dm.globalToLocal(self.vSed, self.vSedLocal)
        if MPIrank == 0 and self.verbose:
            print(
                "Update Sediment Load (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )
        petsc4py.PETSc.garbage_cleanup()

        return

    def _moveDownstream(self, vSed, step):
        """
        In cases where river sediment fluxes drain into depressions, they might fill the sink completely and overspill or be deposited in it. This function computes the excess of sediment (if any) able to flow dowstream.

        .. important::

            The excess sediment volume is then added to the downstream sediment flux (`vSed`).

        :arg vSed: excess sediment volume array
        :arg step: downstream distribution step

        :return: excess (boolean set to True if excess sediment fluxes remain to be distributed)
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

        # In case there is still remaining sediment flux to distribute downstream
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
            self.fMat.destroy()

        return excess

    def _distributeSediment(self, hl):
        """
        This function finds the continental sediment volumes to distribute downstream (down to the ocean or sinks) for a given iteration.

        :arg hl: local elevation prior deposition
        """

        t0 = process_time()

        # Define the volume to distribute
        self.vSedLocal.copy(result=self.QsL)
        self.dm.localToGlobal(self.QsL, self.Qs)

        # Get the volumetric sediment rate (m3/yr) to distribute during the time step and convert it in volume (m3)
        vSed = self.QsL.getArray().copy() * self.dt

        step = 0
        excess = True
        self.fMat.destroy()
        while excess:
            t1 = process_time()
            excess = self._moveDownstream(vSed, step)
            if excess:
                vSed = self.tmpL.getArray().copy()
                nvSed = vSed.copy()
                nvSed[hl < self.sedFilled] = 0.0
                self.tmpL.setArray(nvSed / self.dt)
                vSed[self.idBorders] = 0.0
                self.vSedLocal.axpy(1.0, self.tmpL)
            if MPIrank == 0 and self.verbose:
                print(
                    "Downstream sediment flux computation step %d (%0.02f seconds)"
                    % (step, process_time() - t1),
                    flush=True,
                )
            step += 1
        self.dm.localToGlobal(self.vSedLocal, self.vSed)

        if MPIrank == 0 and self.verbose:
            print(
                "Distribute Continental Sediments (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        return

    def _updateSinks(self, hl):
        """
        This function updates depressions elevations based on incoming sediment volumes. It runs until all inland sediments have reached either a sink which is not completely filled or the open ocean.

        :arg hl: local elevation prior deposition
        """

        # Difference between initial volume and remaining one
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

        # Update soil thickness
        if self.cptSoil:
            self.updateSoilThickness()

        # Update stratigraphic layer parameters
        if self.stratNb > 0:
            self.deposeStrat()
            self.elevStrat()

        # In case there is other sediment type
        self.pitParams[:, 0] = self.pitVol.copy()

        return

    def sedChange(self):
        """
        This function is the main entry point to perform continental river-induced deposition. It calls the private functions:

        - _getSedFlux
        - _distributeSediment
        - _updateSinks

        """

        # Find Continental Sediment Fluxes
        self._getSedFlux()

        # Compute depressions information as elevations changed due to erosion
        self.fillElevation(sed=True)
        hl = self.hLocal.getArray().copy()
        self.lsink = self.lsinki.copy()
        self.pitVol = self.pitParams[:, 0].copy()

        # Distribute inland sediments
        self.sedFilled = hl.copy()
        self._distributeSediment(hl)

        self._updateSinks(hl)

        return
