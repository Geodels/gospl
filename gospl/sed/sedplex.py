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

        :return: (excess, sedFilled_changed) where excess is True when sediment
            still needs to be redistributed, and sedFilled_changed is True when
            at least one pit saturated this iteration (signalling that the
            flow direction matrix should be rebuilt next call).
        """
        excess = False
        sedFilled_changed = False

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
            ids = np.isin(self.pitIDs, np.where(eV > 0.0)[0])
            self.sedFilled[ids] = self.lFill[ids]
            sedFilled_changed = True

        # Update unfilled depressions volumes
        eIDs = eV < 0
        self.pitVol[eIDs] -= inV[eIDs]
        self.pitVol[self.pitVol < 0] = 0.0

        # In case there is still remaining sediment flux to distribute downstream
        if (eV > 1.0e-3).any():
            # Only rebuild the flow direction matrix when the topography has
            # actually changed (a pit saturated this iteration), or on the
            # very first iteration when no matrix exists yet, or on the
            # special step==100 fallback that switches to lFill. Otherwise
            # reuse the matrix built on the previous iteration.
            need_rebuild = (
                sedFilled_changed or step == 0 or step == 100
            )
            if need_rebuild:
                if step > 0 and self.fMat is not None:
                    self.fMat.destroy()
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
            # Note: matrix lifecycle moved to _distributeSediment so we don't
            # destroy/rebuild on every iteration when the topography is stable.

        return excess, sedFilled_changed

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

        # Drop any matrix from the prior context; _moveDownstream will rebuild
        # on iteration 0 and only rebuild thereafter when sedFilled changes.
        self.fMat.destroy()
        self.fMat = None

        step = 0
        excess = True
        max_iters = 5000
        while excess:
            if step >= max_iters:
                if MPIrank == 0:
                    print(
                        "Continental sediment routing did not converge after "
                        "%d iterations; continuing." % max_iters,
                        flush=True,
                    )
                break
            t1 = process_time()
            excess, _ = self._moveDownstream(vSed, step)
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

        # Final cleanup of the cached flow matrix
        if self.fMat is not None:
            self.fMat.destroy()
            self.fMat = None

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

        Sediment is deposited as a bottom-up fill: for each pit we find the
        lake-surface elevation `target_lvl` such that the volume below that
        level (above the bathymetry) equals the deposited volume. Nodes
        already above `target_lvl` get no deposit; nodes below are raised to
        `target_lvl`. This produces a horizontal lake surface and a deposit
        thickness that varies with the bathymetry, which is far more
        physical than the previous uniform "(lFill-hl) * fraction" fill.

        :arg hl: local elevation prior deposition
        """

        # Sediment volume deposited per pit (m^3)
        depo = self.pitParams[:, 0] - self.pitVol
        depo[depo < 0] = 0.0

        # Find target lake surface per pit by interpolating the discrete
        # depth-volume curves built in fillElevation. Prepend a (0, pitMin)
        # point so very small deposits also interpolate correctly (the curve
        # otherwise starts at filled_lvl[p, 0] which is 20% above pit floor).
        num_pits = len(self.pitParams)
        target_lvl = np.full(num_pits, -np.inf, dtype=np.float64)
        active = np.where(depo > 0)[0]
        for p in active:
            pit_min = self.pitParams[p, 1] - self.pitParams[p, 2]
            vols = np.concatenate(([0.0], self.filled_vol[p]))
            lvls = np.concatenate(([pit_min], self.filled_lvl[p]))
            target_lvl[p] = np.interp(depo[p], vols, lvls)

        # Broadcast per-pit target level to nodes; zero deposit outside pits
        node_target_lvl = np.full(self.lpoints, -np.inf, dtype=np.float64)
        in_pit = self.pitIDs >= 0
        node_target_lvl[in_pit] = target_lvl[self.pitIDs[in_pit]]

        # Cap at the spill elevation (rim) so over-fills cannot push the
        # surface above the pit's outlet. Then bottom-up: each node rises to
        # max(hl, surface).
        surface = np.minimum(node_target_lvl, self.lFill)
        delta = np.maximum(0.0, surface - hl)
        delta[~in_pit] = 0.0

        # Update cumulative erosion and deposition as well as elevation
        self.tmpL.setArray(delta)
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
