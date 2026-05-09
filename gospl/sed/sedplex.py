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

    def _bottomUpDelta(self, hl, depo, pit_select):
        """
        Compute per-node deposit thickness for the pits flagged in
        `pit_select` using the bottom-up rule: lake surface stays horizontal
        at the level matching the deposited volume.

        :arg hl: local elevation prior deposition
        :arg depo: per-pit deposited volume (m^3)
        :arg pit_select: boolean array (num_pits) selecting which pits to fill

        :return: per-node delta thickness (m); zero outside selected pits
        """

        num_pits = len(self.pitParams)
        target_lvl = np.full(num_pits, -np.inf, dtype=np.float64)
        active = np.where(pit_select & (depo > 0))[0]
        for p in active:
            pit_min = self.pitParams[p, 1] - self.pitParams[p, 2]
            vols = np.concatenate(([0.0], self.filled_vol[p]))
            lvls = np.concatenate(([pit_min], self.filled_lvl[p]))
            target_lvl[p] = np.interp(depo[p], vols, lvls)

        node_target_lvl = np.full(self.lpoints, -np.inf, dtype=np.float64)
        in_pit = self.pitIDs >= 0
        node_target_lvl[in_pit] = target_lvl[self.pitIDs[in_pit]]

        surface = np.minimum(node_target_lvl, self.lFill)
        delta = np.maximum(0.0, surface - hl)

        # Mask out non-selected pits
        if active.size > 0:
            sel_node = np.isin(self.pitIDs, active)
        else:
            sel_node = np.zeros(self.lpoints, dtype=bool)
        delta[~sel_node] = 0.0
        return delta

    def _spillCoords(self):
        """
        Globally-known XYZ coordinates of each pit's spill point.

        Each rank knows its own spill points (those it owns); we publish the
        full table by zero-filling the others and Allreducing with SUM.

        :return: (num_pits, 3) numpy array of spill coordinates
        """
        num_pits = len(self.pitParams)
        spill_xyz = np.zeros((num_pits, 3), dtype=np.float64)
        spill_idx = self.pitInfo[:, 0].astype(int)
        spill_owner = self.pitInfo[:, 1].astype(int)
        my_pits = np.where(spill_owner == MPIrank)[0]
        if my_pits.size > 0:
            valid = spill_idx[my_pits] >= 0
            mp = my_pits[valid]
            spill_xyz[mp] = self.lcoords[spill_idx[mp]]
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, spill_xyz, op=MPI.SUM)
        return spill_xyz

    def _addPitMicroTilt(self, delta, pit_select):
        """
        Tilt each lake's deposited surface very slightly toward its spill
        point so the next flow-routing step doesn't have to deal with large
        flat areas. The tilt is mass-conservative per pit (mean removed) and
        small enough (~1e-6 m/m) to be physically negligible.

        :arg delta: per-node deposit thickness from bottom-up fill
        :arg pit_select: boolean array (num_pits) selecting tilted pits

        :return: tilted delta (same shape, mass preserved per pit)
        """
        num_pits = len(self.pitParams)
        # pit_select is globally consistent (built from pitParams / pitVol
        # which are MPI-reduced), so this early-exit is collective-safe.
        # We must NOT add an early-exit on sel_node.any() below: that mask is
        # local-only, so ranks with zero locally-selected cells would skip
        # the Allreduce inside _spillCoords and the per-pit reduce while
        # other ranks wait, causing a hang.
        if not pit_select.any():
            return delta

        active_pit = np.where(pit_select)[0]
        sel_node = np.isin(self.pitIDs, active_pit) & (delta > 0)

        spill_xyz = self._spillCoords()

        # Distance from each selected node to its pit's spill
        in_pit = self.pitIDs >= 0
        pit_safe = np.where(in_pit, self.pitIDs, 0)
        diff_xyz = self.lcoords - spill_xyz[pit_safe]
        dist = np.linalg.norm(diff_xyz, axis=1)

        # Per-pit mean distance over selected nodes (MPI reduced)
        # Build local per-pit sums and counts then reduce.
        local_sum = np.zeros(num_pits, dtype=np.float64)
        local_cnt = np.zeros(num_pits, dtype=np.float64)
        sel_local_owned = sel_node & (self.inIDs == 1)
        if sel_local_owned.any():
            ids = self.pitIDs[sel_local_owned]
            grp = npi.group_by(ids)
            up, sums = grp.sum(dist[sel_local_owned])
            _, cnts = grp.sum(np.ones_like(dist[sel_local_owned]))
            local_sum[up] = sums
            local_cnt[up] = cnts
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, local_sum, op=MPI.SUM)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, local_cnt, op=MPI.SUM)
        with np.errstate(divide="ignore", invalid="ignore"):
            mean_dist = np.where(
                local_cnt > 0, local_sum / np.maximum(local_cnt, 1.0), 0.0
            )

        # Slope coefficient: ~1e-6 m/m. Lake surface is higher far from the
        # spill, lower near it -- gentle gradient pointing toward outlet.
        slope_coeff = 1.0e-6
        tilt = slope_coeff * (dist - mean_dist[pit_safe])
        tilt[~sel_node] = 0.0

        # Final delta with tilt; clip negatives in case tilt would push the
        # surface below bedrock (shouldn't happen with such a small slope).
        tilted = delta + tilt
        tilted[tilted < 0] = 0.0
        return tilted

    def _identifyPitInlets(self):
        """
        Identify cells that receive flow from a donor outside their own pit.
        These are the natural inlets where upstream sediment enters a
        depression, used as the initial-pile location for the non-linear
        diffusion path in large/deep partially-filled pits.

        :return: boolean array (lpoints) marking inlet cells (owned cells only)
        """

        if not hasattr(self, "rcvIDi"):
            return np.zeros(self.lpoints, dtype=bool)

        flowDir = self.rcvIDi.shape[1]
        src = np.repeat(np.arange(self.lpoints), flowDir)
        dst = self.rcvIDi.flatten().astype(int)
        # Drop self-loops (cells with no defined receiver)
        good = src != dst
        src = src[good]
        dst = dst[good]

        src_pit = self.pitIDs[src]
        dst_pit = self.pitIDs[dst]
        dst_in_pit = dst_pit >= 0
        cross_pit = src_pit != dst_pit
        owned_dst = self.inIDs[dst] == 1
        inlet_edges = dst_in_pit & cross_pit & owned_dst

        inlet_mask = np.zeros(self.lpoints, dtype=bool)
        if inlet_edges.any():
            inlet_mask[dst[inlet_edges]] = True
        return inlet_mask

    def _diffuseLargePit(self, hl, depo, pit_select):
        """
        Deposit sediment in large/deep partially-filled pits via inlet pile
        plus marine-style non-linear diffusion. Mass is approximately
        conserved by the diffusion solver; a per-pit rescale at the end
        absorbs any boundary drift.

        :arg hl: local elevation prior deposition
        :arg depo: per-pit deposited volume
        :arg pit_select: boolean array (num_pits) selecting diffused pits

        :return: per-node delta thickness (m)
        """

        num_pits = len(self.pitParams)
        active = np.where(pit_select & (depo > 0))[0]
        if active.size == 0:
            return np.zeros(self.lpoints, dtype=np.float64)

        # Build node-level mask for selected pits
        sel_node = np.isin(self.pitIDs, active)

        # Identify inlets among selected pits
        inlet_mask = self._identifyPitInlets() & sel_node

        # Per-pit inlet count (globally)
        inlet_count = np.zeros(num_pits, dtype=np.float64)
        if inlet_mask.any():
            ids = self.pitIDs[inlet_mask]
            grp = npi.group_by(ids)
            up, cnts = grp.sum(np.ones_like(ids, dtype=np.float64))
            inlet_count[up] = cnts
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, inlet_count, op=MPI.SUM)

        # Pits in the active set that have no detected inlet on any rank
        # (rare; can happen for orphan basins with all-flat boundaries) fall
        # back to the highest in-pit cell as a one-cell pile.
        no_inlet = active[inlet_count[active] == 0]
        if no_inlet.size > 0:
            # Batched two-step reduction so the cost is independent of
            # len(no_inlet). Step 1 finds the global max elevation per pit;
            # step 2 finds which rank held that max (lowest rank wins ties).
            # The previous element-wise MAX over [elev, rank, idx] was
            # incorrect because the rank field maxes independently of elev,
            # so only the highest MPIrank ever wrote, regardless of which
            # rank actually owned the highest cell.
            mpisize = MPI.COMM_WORLD.Get_size()
            local_top_elev = np.full(num_pits, -np.inf, dtype=np.float64)
            local_top_idx = np.full(num_pits, -1, dtype=np.int64)
            for p in no_inlet:
                pit_local = np.where(
                    (self.pitIDs == p) & (self.inIDs == 1)
                )[0]
                if pit_local.size > 0:
                    k = np.argmax(hl[pit_local])
                    local_top_elev[p] = hl[pit_local[k]]
                    local_top_idx[p] = pit_local[k]

            # Step 1: max elevation per pit, globally
            global_top_elev = local_top_elev.copy()
            MPI.COMM_WORLD.Allreduce(
                MPI.IN_PLACE, global_top_elev, op=MPI.MAX
            )

            # Step 2: each rank publishes its rank if it owns a cell whose
            # elevation matches the global max for that pit; sentinel mpisize
            # otherwise. MIN reduction picks the lowest-ranked match (stable
            # tie-break).
            winner = np.full(num_pits, mpisize, dtype=np.int64)
            for p in no_inlet:
                if (
                    local_top_idx[p] >= 0
                    and abs(local_top_elev[p] - global_top_elev[p]) < 1.0e-12
                ):
                    winner[p] = MPIrank
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, winner, op=MPI.MIN)

            # Only the winner rank writes the inlet mask + bumps inlet_count
            for p in no_inlet:
                if (
                    int(winner[p]) == MPIrank
                    and local_top_idx[p] >= 0
                ):
                    inlet_mask[int(local_top_idx[p])] = True
                    inlet_count[p] = 1.0
            # Step 3: sync inlet_count so all ranks see the bumps
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, inlet_count, op=MPI.MAX)

        # Initial deposit thickness: depo[pit] / (n_inlets * inlet_area)
        delta_init = np.zeros(self.lpoints, dtype=np.float64)
        if inlet_mask.any():
            iIDs = np.where(inlet_mask)[0]
            for c in iIDs:
                p = self.pitIDs[c]
                if inlet_count[p] > 0:
                    delta_init[c] = depo[p] / (inlet_count[p] * self.larea[c])

        # Run non-linear diffusion over the selected pit cells. Diffusion is
        # zero outside this mask so deposit cannot spread out of the basin.
        delta_smooth = self._diffuseImplicit(
            delta_init, sel_node, self.nl_pit_K, label="lake"
        )

        # Per-pit mass rescale: enforce sum(delta_smooth * larea) == depo[p]
        # Compute global sum per pit
        actual_vol = np.zeros(num_pits, dtype=np.float64)
        owned = (self.inIDs == 1) & sel_node
        if owned.any():
            ids = self.pitIDs[owned]
            vol_per_node = (delta_smooth * self.larea)[owned]
            grp = npi.group_by(ids)
            up, sums = grp.sum(vol_per_node)
            actual_vol[up] = sums
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, actual_vol, op=MPI.SUM)

        with np.errstate(divide="ignore", invalid="ignore"):
            scale = np.where(
                actual_vol > 0, depo / np.maximum(actual_vol, 1e-30), 1.0
            )
        in_pit = self.pitIDs >= 0
        pit_safe = np.where(in_pit, self.pitIDs, 0)
        node_scale = scale[pit_safe]
        node_scale[~sel_node] = 0.0
        delta_smooth = delta_smooth * node_scale

        return delta_smooth

    def _updateSinks(self, hl):
        """
        Update depression elevations based on incoming sediment volumes.

        Three-path algorithm depending on the pit's fill state and size:
          1. Pit fully filled (depo >= pitVolume) -- bottom-up to lFill (rim).
          2. Partial fill, pit small/shallow (volume below `nl_pit_volume` OR
             depth below `nl_pit_depth`) -- bottom-up + per-pit micro-tilt
             toward the spill point so the lake surface is not perfectly flat
             (avoids large flats in the next flow-routing step).
          3. Partial fill, pit large AND deep -- inlet-pile initial condition
             + marine-style non-linear diffusion. Produces emergent
             clinoform-like deposit progradation from the inlets.

        :arg hl: local elevation prior deposition
        """

        depo = self.pitParams[:, 0] - self.pitVol
        depo[depo < 0] = 0.0
        pit_volume = self.pitParams[:, 0]
        pit_depth = self.pitParams[:, 2]

        # Path classification
        active = depo > 0
        full_fill = active & (depo >= pit_volume - 1.0e-6)
        partial = active & ~full_fill
        is_large = (pit_volume >= self.nl_pit_volume) & (
            pit_depth >= self.nl_pit_depth
        )
        diffuse_path = partial & is_large
        bottomup_path = (full_fill | (partial & ~is_large))

        # Path 1+2: bottom-up
        delta = self._bottomUpDelta(hl, depo, bottomup_path)

        if bottomup_path.any():
            delta = self._addPitMicroTilt(delta, bottomup_path)

        # Path 3: large/deep partial fill via diffusion
        if diffuse_path.any():
            delta_diff = self._diffuseLargePit(hl, depo, diffuse_path)
            delta = delta + delta_diff

        # Apply deposit
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
