import os
import gc
import sys
import petsc4py
from gospl.tools.petscgc import safe_garbage_cleanup
import numpy as np
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import setmaxnb

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()


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

        # Dual-lithology fine sediment flux (m3/yr). vSed carries the TOTAL
        # (coarse+fine) sediment; vSedF carries the fine sub-flux, both routed
        # by the same upstream-integration operator. fineFrac is the per-node
        # fine fraction of the routed flux (= arriving-sediment composition at
        # sinks), used by deposeStrat. Allocated unconditionally (negligible)
        # and registered in destroy_DMPlex; only populated when stratLith.
        self.vSedF = self.hGlobal.duplicate()
        self.vSedFLocal = self.hLocal.duplicate()
        self.fineFrac = np.zeros(self.lpoints, dtype=np.float64)
        # depoFineFrac is the per-node fine fraction actually used by
        # deposeStrat. It starts each step as the routed arriving composition
        # (fineFrac) and is overridden inside pits by _pitFineFraction (3b:
        # fine biased to the depocenter, coarse to the margins/inlet).
        # _totFlux snapshots the pre-cascade total routed flux (vSedLocal is
        # later mutated by the pit cascade) so per-pit fine fractions stay valid.
        self.depoFineFrac = np.zeros(self.lpoints, dtype=np.float64)
        # Per-pit fine volume RETAINED through the cascade under coarse-settles-
        # first (filled pits keep coarse, fine overspills); reset each step in
        # _distributeSediment and used by _pitFineFraction. _routedFine carries
        # the fine overspill routed downstream by _moveDownstream.
        self._pitRetFine = np.zeros(1, dtype=np.float64)
        self._routedFine = np.zeros(self.lpoints, dtype=np.float64)

        # In-model provenance tracers (opt-in `provenance:`; DESIGN_PROVENANCE.md
        # §6). One routed sub-flux per source class (vSedP[c], mirroring vSedF),
        # the arriving/deposited per-node composition (provFrac/depoProvFrac),
        # and per-class mass-balance diagnostics. Allocated only when provOn.
        if self.provOn:
            self.vSedP = [self.hGlobal.duplicate() for _ in range(self.provNb)]
            self.vSedPLocal = self.hLocal.duplicate()
            self.provFrac = np.zeros((self.lpoints, self.provNb), dtype=np.float64)
            self.depoProvFrac = np.zeros((self.lpoints, self.provNb), dtype=np.float64)
            self._provEroded = np.zeros(self.provNb, dtype=np.float64)
            self._provDeposited = np.zeros(self.provNb, dtype=np.float64)
            # B2b-pit: per-pit retained provenance (accumulated through the
            # cascade in _moveDownstream, used by _pitProvFraction) and the
            # provenance overspill routed downstream each iteration.
            self._pitRetProv = np.zeros((1, self.provNb), dtype=np.float64)
            self._routedProv = np.zeros((self.lpoints, self.provNb), dtype=np.float64)

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
            # thCoarse is already an erosion-positive rate (m/yr) from
            # stratplex.erodeStrat; use as-is. In dual-lithology mode the
            # transported sediment is the TOTAL (coarse + fine) so the
            # deposition machinery handles all eroded volume.
            if self.stratLith:
                self.tmpL.setArray(self.thCoarse + self.thFine)
            else:
                self.tmpL.setArray(self.thCoarse)
            self.dm.localToGlobal(self.tmpL, self.tmp)
        else:
            # Eb is in the thickness-rate convention (positive deposition,
            # negative incision); the upstream-integration solve below
            # wants an erosion-positive source so that vSed (m^3/yr)
            # accumulates as positive sediment flux downstream. Negate.
            self.Eb.copy(result=self.tmp)
            self.tmp.scale(-1.0)

        # Get the volume of sediment transported in m3/yr
        self.tmp.pointwiseMult(self.tmp, self.areaGlobal)
        self._solve_KSP(False, self.fMati, self.tmp, self.vSed)

        # Dual lithology: route the fine sub-flux through the SAME operator
        # (linear, so this is exact), then snapshot the per-node fine fraction
        # of the routed (upstream-integrated) flux — i.e. the composition of
        # sediment arriving at each node. Done before fMati is destroyed.
        if self.stratNb > 0 and self.stratLith:
            self.tmpL.setArray(self.thFine)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.tmp.pointwiseMult(self.tmp, self.areaGlobal)
            self._solve_KSP(False, self.fMati, self.tmp, self.vSedF)
            self.dm.globalToLocal(self.vSed, self.vSedLocal)
            self.dm.globalToLocal(self.vSedF, self.vSedFLocal)
            vtot = self.vSedLocal.getArray()
            vfin = self.vSedFLocal.getArray()
            self.fineFrac = np.divide(
                vfin, vtot, out=np.zeros_like(vtot), where=vtot > 0.0
            )
            np.clip(self.fineFrac, 0.0, 1.0, out=self.fineFrac)
            # Seed the deposit composition with the arriving (routed) fine
            # fraction; _updateSinks (pits, coarse-settles-first) and seaChange
            # (marine) refine it.
            self.depoFineFrac = self.fineFrac.copy()

        # Provenance: route each source class's eroded sub-flux through the SAME
        # operator (linear, exact, Σ over classes == total). provFrac is the
        # arriving per-node composition; depoProvFrac (no settling bias — a
        # passive tracer) is what deposeStrat lays down. (Through-flux routing;
        # the pit-cascade redistribution refinement is B2b.)
        if self.provOn:
            self.dm.globalToLocal(self.vSed, self.vSedLocal)
            vtot = self.vSedLocal.getArray()
            for c in range(self.provNb):
                self.tmpL.setArray(self.provEro[:, c])
                self.dm.localToGlobal(self.tmpL, self.tmp)
                self.tmp.pointwiseMult(self.tmp, self.areaGlobal)
                self._solve_KSP(False, self.fMati, self.tmp, self.vSedP[c])
                self.dm.globalToLocal(self.vSedP[c], self.vSedPLocal)
                self.provFrac[:, c] = np.divide(
                    self.vSedPLocal.getArray(), vtot,
                    out=np.zeros(self.lpoints), where=vtot > 0.0,
                )
            np.clip(self.provFrac, 0.0, 1.0, out=self.provFrac)
            self.depoProvFrac = self.provFrac.copy()

        self.fMati.destroy()

        # Update local vector
        self.dm.globalToLocal(self.vSed, self.vSedLocal)
        if MPIrank == 0 and self.verbose:
            print(
                "Update Sediment Load (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )
        safe_garbage_cleanup()

        return

    def _moveDownstream(self, vSed, vSedF, step, vSedP=None):
        """
        In cases where river sediment fluxes drain into depressions, they might fill the sink completely and overspill or be deposited in it. This function computes the excess of sediment (if any) able to flow dowstream.

        .. important::

            The excess sediment volume is then added to the downstream sediment flux (`vSed`).

        Dual lithology: the fine sub-volume ``vSedF`` is threaded in lockstep so
        overspill is **fine-enriched** — coarse settles first (retained up to the
        pit capacity), so a filled pit keeps a coarse-enriched deposit and the
        excess that overspills carries the remaining fine. The fine overspill is
        routed through the same flow matrix as the total (linear); the per-pit
        retained fine is accumulated in ``self._pitRetFine`` (used by
        ``_pitFineFraction``) and the routed-downstream fine is returned in
        ``self._routedFine``.

        :arg vSed: excess sediment volume array (total)
        :arg vSedF: excess fine sediment volume array (ignored when single)
        :arg step: downstream distribution step

        :return: (excess, sedFilled_changed) where excess is True when sediment
            still needs to be redistributed, and sedFilled_changed is True when
            at least one pit saturated this iteration (signalling that the
            flow direction matrix should be rebuilt next call).
        """
        excess = False
        sedFilled_changed = False
        dual = self.stratLith
        prov = getattr(self, "provOn", False)

        # Remove points belonging to other processors
        vSed = np.multiply(vSed, self.inIDs)
        if dual:
            vSedF = np.multiply(vSedF, self.inIDs)
        if prov:
            vSedP = vSedP * self.inIDs[:, None]

        # CLOSED terminal sinks: flow-terminal nodes that are not a labelled
        # depression (and, via `lsink`, not an outlet or below sea level). On a
        # closed (wall) domain the sediment arriving here can neither drain nor
        # fill a pit, so the per-pit accounting below (which excludes pit-id -1)
        # would silently drop it. Deposit it in place instead — accumulate it
        # and remove it from the routed flux. Empty set on open/marine domains.
        closed = getattr(self, "_closedDepo", None)
        if closed is not None:
            csink = self.lsink & (self.pitIDs < 0)
            if csink.any():
                self._closedDepo[csink] += vSed[csink]
                vSed = vSed.copy()
                vSed[csink] = 0.0
                if dual:
                    vSedF = vSedF.copy()
                    vSedF[csink] = 0.0
                if prov:
                    vSedP = vSedP.copy()
                    vSedP[csink, :] = 0.0

        # Get volume incoming in each depression
        grp = npi.group_by(self.pitIDs[self.lsink])
        uID = grp.unique
        _, vol = grp.sum(vSed[self.lsink])
        inV = np.zeros(len(self.pitParams), dtype=np.float64)
        ids = uID > -1
        inV[uID[ids]] = vol[ids]
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, inV, op=MPI.SUM)

        inVf = np.zeros(len(self.pitParams), dtype=np.float64)
        if dual:
            _, volf = grp.sum(vSedF[self.lsink])
            inVf[uID[ids]] = volf[ids]
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, inVf, op=MPI.SUM)
            inVf = np.minimum(inVf, inV)

        # Provenance entering each pit (per class). Proportional tracer: the pit
        # retains and overspills each class in proportion to its incoming mix.
        inVp = np.zeros((len(self.pitParams), self.provNb), dtype=np.float64)
        if prov:
            for c in range(self.provNb):
                _, volp = grp.sum(vSedP[self.lsink, c])
                inVp[uID[ids], c] = volp[ids]
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, inVp, op=MPI.SUM)

        nvSed = np.zeros(self.lpoints, dtype=np.float64)
        nvSedF = np.zeros(self.lpoints, dtype=np.float64)
        nvSedP = np.zeros((self.lpoints, self.provNb), dtype=np.float64)

        # Get excess volume to distribute downstream
        eV = inV - self.pitVol
        if (eV > 0.0).any():
            eIDs = eV > 0.0
            spillIDs = self.pitInfo[eIDs, 0]
            localSpill = np.where(self.pitInfo[eIDs, 1] == MPIrank)[0]
            localPts = spillIDs[localSpill]
            nvSed[localPts] = eV[eIDs][localSpill]
            if dual:
                # Coarse settles first: the pit retains its capacity, coarse
                # filling it before fine, so the retained deposit is coarse-
                # enriched and the overspill (eV) is fine-enriched.
                coarse_in = inV - inVf
                ret_coarse = np.minimum(coarse_in, self.pitVol)
                ret_fine = np.clip(self.pitVol - ret_coarse, 0.0, inVf)
                ovf_fine = inVf - ret_fine
                self._pitRetFine[eIDs] += ret_fine[eIDs]
                nvSedF[localPts] = ovf_fine[eIDs][localSpill]
            if prov:
                # Proportional: retained = pitVol fraction of the incoming mix,
                # overspill = the rest (each class in the pit's composition).
                ret_frac = np.divide(
                    self.pitVol, inV, out=np.zeros_like(inV), where=inV > 0.0
                )
                ret_prov = inVp * ret_frac[:, None]
                ovf_prov = inVp - ret_prov
                self._pitRetProv[eIDs] += ret_prov[eIDs]
                nvSedP[localPts] = ovf_prov[eIDs][localSpill]
            self.pitVol[eIDs] = 0.0
            ids = np.isin(self.pitIDs, np.where(eV > 0.0)[0])
            self.sedFilled[ids] = self.lFill[ids]
            sedFilled_changed = True

        # Update unfilled depressions volumes (retain all incoming)
        eIDs = eV < 0
        self.pitVol[eIDs] -= inV[eIDs]
        self.pitVol[self.pitVol < 0] = 0.0
        if dual:
            self._pitRetFine[eIDs] += inVf[eIDs]
        if prov:
            self._pitRetProv[eIDs] += inVp[eIDs]      # unfilled pits retain all

        self._routedFine = np.zeros(self.lpoints, dtype=np.float64)
        self._routedProv = np.zeros((self.lpoints, self.provNb), dtype=np.float64)
        # In case there is still remaining sediment flux to distribute downstream
        if (eV > 1.0e-3).any():  # TODO-REFACTOR: value matches DEPOSIT_FLOOR but distinct role (sediment-routing convergence threshold); do not replace
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
                # seed=True: the sediment cascade solves the SAME substochastic
                # (I - W^T) system as flow accumulation, so on a cold start
                # (first step, no warm tmp1) the runoff/sediment RHS is a valid
                # lower-bound guess. Without it the step-1 solve propagates one
                # hop/iteration over the full network and grinds to max_it on a
                # high-resolution global mesh (~10x slower first sediment step).
                self._solve_KSP(True, self.fMat, self.tmp, self.tmp1, seed=True)
                self.dm.globalToLocal(self.tmp1, self.tmpL)
                if dual:
                    # Route the fine overspill through the SAME matrix (linear).
                    self.nQs.setArray(nvSedF)
                    self.dm.localToGlobal(self.nQs, self.tmp)
                    self._solve_KSP(True, self.fMat, self.tmp, self.tmp1, seed=True)
                    self.dm.globalToLocal(self.tmp1, self.nQs)
                    self._routedFine = self.nQs.getArray().copy()
                if prov:
                    # Route each class's overspill through the SAME matrix.
                    for c in range(self.provNb):
                        self.nQs.setArray(nvSedP[:, c])
                        self.dm.localToGlobal(self.nQs, self.tmp)
                        self._solve_KSP(True, self.fMat, self.tmp, self.tmp1, seed=True)
                        self.dm.globalToLocal(self.tmp1, self.nQs)
                        self._routedProv[:, c] = self.nQs.getArray().copy()
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

        # Dual lithology: distribute the fine sub-volume in lockstep so the
        # coarse-settles-first overspill in _moveDownstream is conservative.
        # vSedFLocal mirrors vSedLocal EXACTLY (it is NOT zeroed): it keeps the
        # initial fine flux and accumulates the residual overspill, so at the
        # marine sinks it ends up as the post-cascade fine reaching the sea.
        if self.stratLith:
            vSedF = self.vSedFLocal.getArray().copy() * self.dt
            self._pitRetFine = np.zeros(len(self.pitParams), dtype=np.float64)
        else:
            vSedF = None

        # Provenance: thread each class's sub-flux through the cascade so a pit's
        # retained mix (and its overspill into downstream pits) is exact.
        # _pitRetProv accumulates the per-pit retained provenance for
        # _pitProvFraction; no separate accumulator is needed.
        if self.provOn:
            vSedP = np.zeros((self.lpoints, self.provNb), dtype=np.float64)
            for c in range(self.provNb):
                self.dm.globalToLocal(self.vSedP[c], self.vSedPLocal)
                vSedP[:, c] = self.vSedPLocal.getArray() * self.dt
            self._pitRetProv = np.zeros(
                (len(self.pitParams), self.provNb), dtype=np.float64
            )
        else:
            vSedP = None

        # Drop any matrix from the prior context; _moveDownstream will rebuild
        # on iteration 0 and only rebuild thereafter when sedFilled changes.
        self.fMat.destroy()
        self.fMat = None

        # Accumulator for sediment that reaches CLOSED terminal sinks (no pit, no
        # outlet, above sea) during the cascade — deposited after the loop so a
        # closed (wall) domain conserves mass. Empty on open/marine domains.
        self._closedDepo = np.zeros(self.lpoints, dtype=np.float64)

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
            excess, _ = self._moveDownstream(vSed, vSedF, step, vSedP=vSedP)
            if excess:
                vSed = self.tmpL.getArray().copy()
                nvSed = vSed.copy()
                nvSed[hl < self.sedFilled] = 0.0
                self.tmpL.setArray(nvSed / self.dt)
                vSed[self.outletIDs] = 0.0
                self.vSedLocal.axpy(1.0, self.tmpL)
                if self.stratLith:
                    # Mirror the residual-flux accumulation for the fine pile.
                    vSedF = self._routedFine.copy()
                    nvSedF = vSedF.copy()
                    nvSedF[hl < self.sedFilled] = 0.0
                    self.tmpL.setArray(nvSedF / self.dt)
                    vSedF[self.outletIDs] = 0.0
                    self.vSedFLocal.axpy(1.0, self.tmpL)
                if self.provOn:
                    # Carry the routed provenance overspill to the next pit
                    # (downstream-lake chains), matching the total's treatment.
                    vSedP = self._routedProv.copy()
                    vSedP[self.outletIDs, :] = 0.0
            if MPIrank == 0 and self.verbose:
                print(
                    "Downstream sediment flux computation step %d (%0.02f seconds)"
                    % (step, process_time() - t1),
                    flush=True,
                )
            step += 1

        # Diagnostic (verbose): continental-routing cascade iteration count
        # (number of spillover rounds), surfaced only when notably high.
        if MPIrank == 0 and self.verbose and step > 20:
            print(
                "[sed] continental routing cascade ran %d iterations" % step,
                flush=True,
            )

        # Final cleanup of the cached flow matrix
        if self.fMat is not None:
            self.fMat.destroy()
            self.fMat = None

        # Conservation closure for CLOSED (wall) sinks: deposit the sediment the
        # cascade routed to closed terminal sinks that are not a labelled pit
        # (collected in `_closedDepo` inside _moveDownstream — flow-terminal
        # nodes that are not a pit, an outlet, or below sea level, e.g. the low
        # corner of an all-wall box). No-op on open / marine domains. NOTE: this
        # does NOT make a fully-closed multi-step domain conserve once its basins
        # saturate — see the "wall" boundary note in the docs / AGENTS.
        # Whether ANY rank has closed-sink sediment to deposit must be decided
        # collectively: the block below runs DM scatters (localToGlobal /
        # globalToLocal), which are collective over the comm. On an open domain
        # a closed flow-terminal node (interior pit-less low point) can land on
        # one rank but not another, so a rank-local `_closedDepo.any()` guard
        # would let one rank enter the scatters while another skips ahead to
        # `_updateSinks` -> deadlock. Reduce the flag with LOR so every rank
        # takes the same branch; ranks with no closed sink contribute zeros.
        has_closed = self._closedDepo is not None and bool(self._closedDepo.any())
        has_closed = MPI.COMM_WORLD.allreduce(has_closed, op=MPI.LOR)
        if has_closed:
            delta = np.zeros(self.lpoints)
            nz = self._closedDepo > 0.0
            delta[nz] = self._closedDepo[nz] / self.larea[nz]   # volume -> thickness
            self.tmpL.setArray(delta)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.cumED.axpy(1.0, self.tmp)
            self.hGlobal.axpy(1.0, self.tmp)
            self.dm.globalToLocal(self.cumED, self.cumEDLocal)
            self.dm.globalToLocal(self.hGlobal, self.hLocal)

        self.dm.localToGlobal(self.vSedLocal, self.vSed)
        if self.stratLith:
            self.dm.localToGlobal(self.vSedFLocal, self.vSedF)

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
        Deposit sediment in large/deep partially-filled pits.

        The initial deposit shape is a bathymetric bottom-up fill (most of the
        volume distributed proportional to water depth below the target lake
        level) combined with a small inlet bias that seeds delta progradation.
        Marine-style non-linear diffusion then refines the surface, and a
        per-pit rescale at the end enforces exact volume conservation.

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

        # Initial deposit shape: (1 - bias_frac) of each pit's volume is
        # distributed proportional to local water depth below the target lake
        # level (a bottom-up bathymetric fill), and bias_frac is concentrated
        # at the inlets to seed delta progradation. Starting the diffusion
        # solver from a near-final lake-floor geometry avoids the rim-only
        # blocky artefact produced by a pure inlet-point pile. The fraction
        # is configurable via the `nlPitInletBias` YAML key (clamped [0, 1]).
        bias_frac = self.nl_pit_inlet_bias

        # Target lake water level per active pit (same interp the bottom-up
        # path uses; volume -> level mapping built when the pit was labelled).
        target_lvl = np.full(num_pits, -np.inf, dtype=np.float64)
        for p in active:
            pit_min = self.pitParams[p, 1] - self.pitParams[p, 2]
            vols = np.concatenate(([0.0], self.filled_vol[p]))
            lvls = np.concatenate(([pit_min], self.filled_lvl[p]))
            target_lvl[p] = np.interp(depo[p], vols, lvls)

        in_pit = self.pitIDs >= 0
        pit_safe = np.where(in_pit, self.pitIDs, 0)
        node_lvl = np.where(in_pit, target_lvl[pit_safe], -np.inf)
        water_depth = np.maximum(0.0, node_lvl - hl)
        water_depth[~sel_node] = 0.0

        # Per-pit sum of water_depth * area over owned cells, MPI-reduced.
        denom = np.zeros(num_pits, dtype=np.float64)
        owned = (self.inIDs == 1) & sel_node
        if owned.any():
            ids = self.pitIDs[owned]
            grp = npi.group_by(ids)
            up, sums = grp.sum((water_depth * self.larea)[owned])
            denom[up] = sums
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, denom, op=MPI.SUM)

        # Bathymetric baseline: thick where the lake is deep, thin near rim.
        delta_init = np.zeros(self.lpoints, dtype=np.float64)
        denom_node = denom[pit_safe]
        ok = (denom_node > 0) & sel_node
        delta_init[ok] = (
            (1.0 - bias_frac)
            * depo[pit_safe[ok]]
            * water_depth[ok]
            / np.maximum(denom_node[ok], 1e-30)
        )

        # Inlet bias: bias_frac * depo[p] / (n_inlets * area). For pits where
        # the bathymetric baseline is degenerate (denom == 0) fall back to
        # the original inlet-only behaviour so all the volume is still placed.
        bias_at_inlet = np.where(denom > 0, bias_frac, 1.0)
        if inlet_mask.any():
            iIDs = np.where(inlet_mask)[0]
            for c in iIDs:
                p = self.pitIDs[c]
                if inlet_count[p] > 0:
                    delta_init[c] += (
                        bias_at_inlet[p]
                        * depo[p]
                        / (inlet_count[p] * self.larea[c])
                    )

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

    def _pitFineFraction(self, hl, delta, depo):
        """
        Dual-lithology (3b + fine-enriched overspill): set the per-node fine
        fraction of the continental pit/lake deposit so fine settles toward the
        depocenter and coarse builds the inlet/margins, conserving each pit's
        RETAINED fine volume.

        Composition only — the deposit geometry ``delta`` is NOT modified, so
        model dynamics (elevations, flow) are untouched; the worst case of a
        bug here is a wrong recorded composition, caught by the per-fraction
        invariant + fine-conservation tests.

        Mechanism:
          - Per-pit retained fine fraction ``ff_pit = _pitRetFine / depo`` —
            the fine actually kept in the pit under coarse-settles-first
            (accumulated through the cascade by ``_moveDownstream``); a pit that
            overspilled keeps a coarse-enriched deposit, the excess fine having
            travelled downstream.
          - Bathymetric depth ``d = max(lFill − hl, 0)`` is the depocenter
            proxy (deep at the centre, shallow at the rim/inlet — independent
            of where the deposit piled).
          - Fine fraction is biased ∝ d and renormalised to the deposit-weighted
            mean depth, so ``Σ(delta·larea·ffrac) == ff_pit·Σ(delta·larea) ==``
            the fine volume retained in the pit.

        :arg hl: pre-deposition local elevation
        :arg delta: per-node continental deposit thickness (local)
        :arg depo: per-pit deposited (retained) volume (m^3)
        """
        num_pits = len(self.pitParams)
        in_pit = self.pitIDs >= 0
        owned = (self.inIDs == 1) & in_pit

        # ---- per-pit retained fine fraction (coarse-settles-first) ----
        ff_pit = np.divide(
            self._pitRetFine, depo, out=np.zeros(num_pits), where=depo > 0
        )
        np.clip(ff_pit, 0.0, 1.0, out=ff_pit)

        # ---- deposit-weighted mean bathymetric depth per pit ----
        depth = np.maximum(self.lFill - hl, 0.0)
        dw = np.zeros(num_pits, dtype=np.float64)   # Σ delta*larea
        dwd = np.zeros(num_pits, dtype=np.float64)  # Σ delta*larea*depth
        if owned.any():
            ids = self.pitIDs[owned]
            grp = npi.group_by(ids)
            up = grp.unique
            dl = (delta * self.larea)[owned]
            dw[up] = grp.sum(dl)[1]
            dwd[up] = grp.sum(dl * depth[owned])[1]
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, dw, op=MPI.SUM)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, dwd, op=MPI.SUM)

        depthbar = np.divide(dwd, dw, out=np.zeros(num_pits), where=dw > 0)

        # ---- per-node biased fine fraction (depocenter-weighted) ----
        pit_safe = np.where(in_pit, self.pitIDs, 0)
        ffp = ff_pit[pit_safe]
        db = depthbar[pit_safe]
        ffrac = np.where(
            in_pit & (db > 0.0), ffp * depth / np.maximum(db, 1.0e-30), 0.0
        )
        np.clip(ffrac, 0.0, 1.0, out=ffrac)
        self.depoFineFrac[in_pit] = ffrac[in_pit]

        return

    def _pitProvFraction(self, depo):
        """
        Provenance (B2b-pit): set the per-node source composition of the
        continental pit/lake deposit to each pit's RETAINED provenance mix,
        tracked through the overspill cascade by ``_moveDownstream``.

        Unlike the dual-lithology fine fraction there is **no depocenter bias**
        — provenance is a passive label, so every node in a pit records the same
        retained mix. ``_pitRetProv[pit, :]`` sums over classes to the pit's
        retained volume ``depo[pit]`` (the cascade conserves each class), so the
        per-pit fractions ``_pitRetProv / depo`` sum to 1 and ``depoProvFrac``
        stays summed-to-1 — ``stratP`` partitions ``stratH`` machine-exactly.

        Composition only — the deposit geometry is untouched.

        :arg depo: per-pit deposited (retained) volume (m^3)
        """
        num_pits = len(self.pitParams)
        in_pit = self.pitIDs >= 0

        # Per-pit retained provenance fraction (overspill-chain aware).
        ffp = np.divide(
            self._pitRetProv,
            depo[:, None],
            out=np.zeros((num_pits, self.provNb)),
            where=depo[:, None] > 0,
        )
        # Guard rounding: renormalise rows with deposit to sum to 1.
        rsum = ffp.sum(axis=1)
        good = rsum > 0
        ffp[good] /= rsum[good][:, None]

        pit_safe = np.where(in_pit, self.pitIDs, 0)
        self.depoProvFrac[in_pit, :] = ffp[pit_safe][in_pit, :]

        return

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

        # Dual lithology (3b): bias the deposit composition within pits — fine
        # to the depocenter, coarse to the inlet/margins, using the cascade-
        # tracked retained fine (coarse-settles-first). Geometry (delta) is
        # unchanged; only self.depoFineFrac (read by deposeStrat) is set.
        if self.stratLith:
            self._pitFineFraction(hl, delta, depo)

        # Provenance (B2b-pit): record each pit's retained source mix (uniform
        # within the pit; cascade-tracked overspill keeps downstream-lake chains
        # exact). Geometry unchanged; only self.depoProvFrac (read by
        # deposeStrat) is set.
        if self.provOn:
            self._pitProvFraction(depo)

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
