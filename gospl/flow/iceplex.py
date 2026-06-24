import os
import gc
import petsc4py
import numpy as np

from mpi4py import MPI
from time import process_time

from gospl.tools.constants import (
    BOUNDARY_FLOW_SENTINEL,
    MISSING_DATA_SENTINEL,
)

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import ice_velocity
    from gospl._fortran import ice_lateral_erosion
    from gospl._fortran import mfdreceivers
    from gospl._fortran import epsfill

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()


class IceMesh(object):
    r"""
    Diagnostic glacial-erosion model.

    Rather than time-integrating an ice-dynamics PDE, this class derives a cheap,
    stable diagnostic of the glacial state each step and uses it to drive glacial
    erosion, sediment transport and ice loading. For the equilibrium-line-altitude
    (ELA) surface mass balance :math:`\dot{m}`:

    1. The accumulation :math:`\dot{m}^+` is routed downhill on the
       (epsilon-filled) bed using a multiple-flow-direction algorithm — the same
       machinery as the river flow accumulation — into an **ice discharge**
       :math:`Q` (m\ :sup:`3`/yr); one linear solve, no time integration.
    2. The **ice thickness** follows a Bahr width–area scaling of the discharge,
       :math:`H = e_h\, f_w\, Q^{0.3}`.
    3. The **basal sliding velocity** comes from Glen's sliding law on that
       thickness and the bed-surface slope (the ``ice_velocity`` kernel,
       :math:`u_b \propto H^{n-1}|\nabla s|^{n-1}\nabla s`) — physically bounded,
       so erosion does not spike at flow-convergence cells.

    From these the model computes velocity-based glacial abrasion
    (:math:`E_g = K_g |u_b|^l`, with an optional lateral valley-wall term for
    U-shaping), ice-transported till deposited as moraine where the ice melts out
    (conserving rock volume), the discharge-conserving glacial meltwater
    re-injected into the river flow accumulation, and the ice load feeding the
    flexural isostasy.

    The diagnostic is robust and physical at any resolution and over goSPL's long
    (:math:`10^2`–:math:`10^4` yr) timesteps, where a true ice-dynamics solve is
    stiff and over-thickens km-scale continental ice (see
    ``docs/DESIGN_ICE_SHEET.md`` for the history).

    Useful links:
    - `Braun et al. (1999) <https://doi.org/10.3189/172756499781821797>`_
    - `Herman & Braun (2008) <https://doi.org/10.1029/2007JF000807>`_
    - `Egholm (2012) <https://doi.org/10.1016/j.geomorph.2011.12.019>`_
    - `Deal & Prasicek (2020) <https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020GL089263>`_
    - `Hergarten (2021) <http://dx.doi.org/10.5194/esurf-2021-1>`_
    - `Liebl et al. (2023) <https://gmd.copernicus.org/articles/16/1315/2023/>`_
    """

    def __init__(self, *args, **kwargs):
        """
        Initialisation of the `IceMesh` class.

        This method initializes the ice-related fields (ice thickness, basal
        velocity, meltwater and flexural response) based on the configuration
        flags `iceOn` and `flexOn`. Memory is allocated for these fields.
        """

        if self.iceOn:
            # Diagnostic ice thickness (Bahr width-area scaling of the discharge).
            self.iceHL = self.hLocal.duplicate()
            # Precip-scaled ablation pattern used as the till melt-out / moraine
            # deposition weight.
            self.iceMeltL = self.hLocal.duplicate()
            # Glacial meltwater delivered to the rivers (discharge-conserving),
            # re-injected into the river FA source in flowplex so glacial
            # discharge is not lost downstream. Distinct from iceMeltL.
            self.iceMeltRiverL = self.hLocal.duplicate()
            if self.flexOn:
                self.iceFlex = self.hLocal.duplicate()
            # Basal sliding speed (m/yr); the driver for the velocity-based
            # glacial abrasion law. Registered in destroy_DMPlex.
            self.iceUbL = self.hLocal.duplicate()
            # Glacial abrasion rate E_g = Kg|u_b|^l (m/yr); a diagnostic output
            # field, populated every step from the basal velocity (zero where
            # abrasion is off, i.e. Kg = 0).
            self.iceAbrL = self.hLocal.duplicate()
            # Ice discharge (m^3/yr): the ELA accumulation routed downhill;
            # drives the Bahr thickness (and hence the sliding velocity).
            self.iceFAL = self.hLocal.duplicate()
            self.iceFAG = self.hGlobal.duplicate()
            self.iceHL.set(0.0)
            self.iceMeltL.set(0.0)
            self.iceMeltRiverL.set(0.0)
            self.iceUbL.set(0.0)
            self.iceAbrL.set(0.0)
            self.iceFAL.set(0.0)
            self.iceFAG.set(0.0)
            # Glacial-till mass-balance diagnostics (m^3, owned-node running
            # totals reduced by the till-conservation test).
            self._tillEroded = 0.0
            self._tillDeposited = 0.0

        return

    def _iceMassBalance(self, elaH, iceH):
        """
        Bed elevation `zbed` and the ELA surface mass balance `mdot`
        (m ice/yr, no `meltfac`): accumulation above the ELA, ablation below.
        """
        zbed = self.hLocal.getArray().copy()              # bed elevation (fixed this step)
        # ELA mass-balance ramp (1 above the ice cap, 0 at the ELA, negative
        # below). elaH/iceH may be scalars (uniform) or per-vertex arrays
        # (spatial maps). Where iceH <= elaH (degenerate / no-ice-band cell) the
        # ramp is zeroed so the division never blows up; this is exact for the
        # uniform case, where the caller's guard already ensures iceH > elaH.
        denom = iceH - elaH
        safe = np.where(denom > 0.0, denom, 1.0)
        ramp = np.where(denom > 0.0, (zbed - elaH) / safe, 0.0)
        ramp = np.minimum(ramp, 1.0)
        mdot = self.rainVal * ramp                        # surface mass balance (m ice/yr)
        # Scale / cap the ACCUMULATION (positive mdot) to realistic ice rates:
        # full precipitation is rarely all snow/ice, so `accum_factor` converts
        # precipitation to ice accumulation and `accum_max` caps it. Ablation
        # (negative mdot) is left untouched. Defaults (1.0, None) are a no-op.
        if self.ice_accum_factor != 1.0 or self.ice_accum_max is not None:
            acc = np.maximum(mdot, 0.0) * self.ice_accum_factor
            if self.ice_accum_max is not None:
                acc = np.minimum(acc, self.ice_accum_max)
            mdot = np.where(mdot > 0.0, acc, mdot)
        return zbed, mdot

    def _matrixIceFlow(self, dir_ice=1, terminus=None):
        """
        Build the multiple-flow-direction (MFD) routing matrix for ice on the
        (epsilon-filled) bed surface — used by the diagnostic ``mfd`` flow model
        to accumulate the ELA mass balance into an ice discharge. Mirrors the
        water flow-direction matrix; deterministic slope tie-break via ``gid``.

        The eps-fill is seeded from the glacier **terminus** (``terminus``, the
        physical ice outlet — ice melts out at/below it), and every cell at or
        below the terminus is made an absorbing sink in the routing matrix. This
        is essential for a well-posed solve: seeding/draining the whole bed down
        to sea level instead leaves the ice-free sub-terminus region as interior
        drainage targets that are NOT absorbed, so ``(I − Wᵀ)`` is singular, the
        discharge KSP diverges, and the bounded fallback zeroes it — producing
        no ice anywhere. Anchoring the outlet at the terminus drains the
        above-terminus ice region to a well-posed boundary.
        """
        hl = self.hLocal.getArray().copy()
        if terminus is None:
            terminus = max(self.hGlobal.min()[1] + 0.1, self.sealevel)
        # Absorbing outlet = the ice-free melt-out region (bed at/below the
        # terminus) plus the domain borders. Captured on the BED elevation
        # (matches the `fa[zbed < terminus] = 0` discharge mask downstream).
        sinks = hl <= terminus
        if self.flatModel:
            hl[self.idBorders] = BOUNDARY_FLOW_SENTINEL
            sinks[self.idBorders] = True
        fillz = np.zeros(self.mpoints, dtype=np.float64) + MISSING_DATA_SENTINEL
        fillz[self.locIDs] = hl
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, fillz, op=MPI.MAX)
        if MPIrank == 0:
            fillz = epsfill(terminus, fillz)
        fillEPS = MPI.COMM_WORLD.bcast(fillz, root=0)
        fillz = fillEPS[self.locIDs]

        rcv, _, wght = mfdreceivers(
            dir_ice, 1.0, fillz, BOUNDARY_FLOW_SENTINEL, self.gid,
        )
        sink_ids = np.where(sinks)[0]
        rcv[sink_ids, :] = sink_ids[:, None]
        wght[sink_ids, :] = 0.0

        self.iceMat = self.iMat.copy()
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        nodes = indptr[:-1]
        for k in range(dir_ice):
            tmpMat = self._matrix_build()
            data = wght[:, k].copy()
            data[rcv[:, k].astype(petsc4py.PETSc.IntType) == nodes] = 0.0
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr, rcv[:, k].astype(petsc4py.PETSc.IntType), data,
            )
            tmpMat.assemblyEnd()
            self.iceMat.axpy(-1.0, tmpMat)
            tmpMat.destroy()
        if self.memclear:
            del data, indptr, nodes, hl, fillz, fillEPS, rcv, wght
            gc.collect()
        self.iceMat.transpose()
        return

    def _iceFlowMFD(self, elaH, iceH, iceT):
        r"""
        Diagnostic glacial driver (no ice-dynamics solve).

        1. ELA surface mass balance (``_iceMassBalance``: accumulation above the
           ELA, ablation below — including the ``accum_factor``/``accum_max``
           controls).
        2. Route the accumulation downhill through the MFD matrix into an **ice
           discharge** ``iceFA`` (one linear solve; no time integration).
        3. **Thickness** ``iceHL`` from a Bahr width–area scaling of the discharge.
        4. **Basal sliding velocity** ``iceUbL`` from Glen's sliding law on that
           thickness + bed-surface slope (``ice_velocity``), the abrasion driver
           ``E_g = K_g |u_b|^l``. Physically bounded (``∝ H^{n-1}|∇s|^{n-1}∇s``),
           unlike a raw balance velocity ``Q/(H·W)`` which blows up with catchment
           size and spikes at flow-convergence cells.

        No stiffness, no CFL — one routing solve per step.
        """
        zbed, mdot = self._iceMassBalance(elaH, iceH)
        terminus = max(iceT, self.sealevel)

        # (2) Route the positive mass balance (accumulation) into a discharge.
        # Source is a VOLUME rate (m^3/yr) = accumulation rate x cell area, so
        # the routed `iceFA` is an ice discharge (m^3/yr).
        self._matrixIceFlow(self.iceDir, terminus)
        iceA = np.maximum(mdot, 0.0) * self.larea
        self.tmpL.setArray(iceA)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        # seed=True: cold-start the ice discharge from the accumulation source
        # (a valid lower bound for the substochastic (I - W^T) routing system,
        # exactly as the water flow-accumulation solve does). Without it the
        # first-step solve starts from iceFAG == 0, cannot propagate across the
        # network within the iteration budget, fails, and is zeroed by the
        # bounded fallback -- leaving NO ice discharge (hence no thickness) for
        # the whole run.
        self._solve_KSP(True, self.iceMat, self.tmp, self.iceFAG, seed=True)
        self.dm.globalToLocal(self.iceFAG, self.iceFAL)
        fa = self.iceFAL.getArray()
        fa[fa < 0.0] = 0.0
        fa[zbed < terminus] = 0.0
        self.iceFAL.setArray(fa)
        self.dm.localToGlobal(self.iceFAL, self.iceFAG)

        # Ablation pattern for the till melt-out / deposition weight (precip-
        # scaled, gated to ice presence — what the till machinery expects).
        self.iceMeltL.setArray(np.maximum(-mdot, 0.0) * self.larea * (fa > 1.0e-8))
        # Glacial meltwater delivered to the rivers (discharge-conserving by
        # default: the routed accumulation released where the ice melts out, so
        # Σ river-melt == Σ accumulation; closes the glacial water budget).
        self.iceMeltRiverL.setArray(self._glacialMeltwater(zbed, mdot))

        # Smooth the discharge (robustness / morphology), then Bahr thickness.
        self.iceFAG.copy(result=self.tmp1)
        smth = self._hillSlope(smooth=1)
        smth[smth < 0.0] = 0.0
        self.iceFAL.setArray(smth)
        self.dm.localToGlobal(self.iceFAL, self.iceFAG)

        # (3) Ice thickness: Bahr width–area scaling of the discharge.
        H = self.icewe * self.icewf * np.power(smth, 0.3)
        H[zbed < terminus] = 0.0
        H[H < 1.0e-1] = 0.0
        self.iceHL.setArray(H)

        # (4) Basal sliding velocity from the diagnostic thickness + bed-surface
        # slope, via Glen's sliding law (``ice_velocity``: u_s ∝
        # H^(n-1)|∇s|^(n-1)∇s). This is PHYSICALLY BOUNDED by thickness and slope
        # — unlike a raw balance velocity Q/(H·W), which blows up with catchment
        # size (Q grows with the upstream area while the cell width does not) and
        # spikes at MFD flow-convergence cells, driving runaway abrasion and
        # tens-of-km erosion/deposition spikes in downstream sinks.
        ub = ice_velocity(self.lpoints, H, zbed, self.ice_slide, self.ice_glen)
        ub[H <= 1.0] = 0.0
        self.iceUbL.setArray(ub)
        abr = self.ice_Kg * np.power(np.maximum(ub, 0.0), self.ice_abr_l)
        abr[self.seaID] = 0.0
        self.iceAbrL.setArray(abr)

        return

    def iceAccumulation(self):
        """
        Main glacial update.

        Computes the diagnostic glacial state (:meth:`_iceFlowMFD`) from the
        ice-cap altitude, the equilibrium-line altitude (ELA) and the glacier
        terminus. Each of these may be a uniform scalar (time-varying) or a
        per-vertex map — the latter for global models where the ELA varies
        strongly with latitude. It also snapshots the ice load for the flexural
        isostasy when enabled.
        """

        ti = process_time()

        # Snapshot the current ice load for flexure (applyFlexure loads the
        # increment iceHL - iceFlex), taken BEFORE the thickness is updated.
        if self.flexOn:
            self.iceHL.copy(result=self.iceFlex)

        # Resolve the glacier-geometry altitudes: per-vertex map where supplied
        # (spatial ELA for global tropical-vs-polar runs), otherwise the uniform
        # time-varying scalar. `spatial` flags that any field is a map.
        iceH = self.iceMesh[self.locIDs] if self.iceMesh is not None else self.iceH(self.tNow)
        elaH = self.elaMesh[self.locIDs] if self.elaMesh is not None else self.elaH(self.tNow)
        iceT = self.termMesh[self.locIDs] if self.termMesh is not None else self.iceT(self.tNow)
        spatial = (
            self.iceMesh is not None
            or self.elaMesh is not None
            or self.termMesh is not None
        )

        # Uniform case keeps the cheap global short-circuit (byte-identical):
        # no ice when the whole surface is below the ELA, or for a degenerate
        # config where the ice-cap altitude is not above the ELA. With spatial
        # maps these are per-node and handled robustly in _iceMassBalance.
        if not spatial:
            max_elev = self.hGlobal.max()[1]
            if max_elev < elaH or iceH <= elaH:
                self.iceHL.set(0.)
                self.iceUbL.set(0.)
                self.iceAbrL.set(0.)
                self.iceMeltL.set(0.)
                self.iceMeltRiverL.set(0.)
                self.iceFAL.set(0.)
                self.iceFAG.set(0.)
                if self.flexOn:
                    self.iceFlex.set(0.)
                if MPIrank == 0 and self.verbose:
                    why = ("max elevation %.0f m below ELA %.0f m" % (max_elev, elaH)
                           if max_elev < elaH
                           else "ice-cap altitude %.0f m not above ELA %.0f m" % (iceH, elaH))
                    print("Glaciers Accumulation (%0.02f seconds) — no ice (%s)"
                          % (process_time() - ti, why), flush=True)
                return

        # Diagnostic glacial driver (routing proxy — glacial erosion morphology
        # without an ice-dynamics solve).
        self._iceFlowMFD(elaH, iceH, iceT)

        # At the first step, seed the flexure reference with the ice just
        # diagnosed so the initial load is not applied as a transient.
        if self.flexOn and self.tNow == self.tStart:
            self.iceHL.copy(result=self.iceFlex)

        # Instructive diagnostics (only with -v / verbose). The reductions are
        # COLLECTIVE so they run on every rank (gated only by the uniform
        # `self.verbose`, never by MPIrank); the formatted line prints on rank 0.
        if self.verbose:
            comm = MPI.COMM_WORLD
            own = self.inIDs == 1
            H = self.iceHL.getArray()
            ub = self.iceUbL.getArray()
            abr = self.iceAbrL.getArray()
            vol = comm.allreduce(float(np.sum(H[own] * self.larea[own])), op=MPI.SUM)
            hmax = comm.allreduce(float(H[own].max()) if own.any() else 0.0, op=MPI.MAX)
            ncov = comm.allreduce(int((H[own] > 0.1).sum()), op=MPI.SUM)
            ntot = comm.allreduce(int(own.sum()), op=MPI.SUM)
            ubmax = comm.allreduce(float(ub[own].max()) if own.any() else 0.0, op=MPI.MAX)
            abrmax = comm.allreduce(float(abr[own].max()) if own.any() else 0.0, op=MPI.MAX)
            if MPIrank == 0:
                cov = 100.0 * ncov / max(ntot, 1)
                geom = (
                    "ELA %.0f m, terminus %.0f m, ice-cap %.0f m"
                    % (elaH, max(iceT, self.sealevel), iceH)
                    if not spatial else "spatial ELA/terminus/ice-cap maps"
                )
                print(
                    "Glaciers Accumulation (%0.02f seconds) — %s"
                    % (process_time() - ti, geom),
                    flush=True,
                )
                print(
                    "   ice volume %.4g km3 | max thickness %.1f m | cover %.1f%% "
                    "(%d cells) | max sliding %.3g m/yr | max abrasion %.3g m/yr"
                    % (vol / 1.0e9, hmax, cov, ncov, ubmax, abrmax),
                    flush=True,
                )

        return

    def _glacialLateralErosion(self):
        r"""
        Lateral glacial erosion rate (m/yr) — valley-wall abrasion by adjacent
        fast ice, the explicit term that widens glaciated valleys toward a
        U-profile (vertical abrasion alone only deepens the trough).

        Returns a per-cell rate (≥0), nonzero only on subaerial 'wall' cells
        (little ice of their own) standing within the ice column of a
        faster-sliding neighbour — see the ``ice_lateral_erosion`` kernel. No-op
        (zeros) unless ``ice.abrasion.Kl > 0``.
        """
        if self.ice_Kl <= 0.0:
            return np.zeros(self.lpoints, dtype=np.float64)
        H = self.iceHL.getArray()
        ub = self.iceUbL.getArray()
        zbed = self.hLocal.getArray()
        elat = ice_lateral_erosion(
            self.lpoints, H, zbed, np.maximum(ub, 0.0),
            self.ice_Kl, self.ice_lat_l, 1.0,   # heps = 1 m ice-presence threshold
        )
        elat[self.seaID] = 0.0                  # no marine erosion
        return elat

    def glacialTill(self):
        r"""
        Glacial till — production, transport and deposition.

        Glacial abrasion (:math:`E_g = K_g |u_b|^l`) erodes rock under sliding
        ice; the abraded material becomes **till** carried by the ice and
        deposited where the ice **melts out** — the ablation zone / terminal
        moraine. The bed is therefore **lowered** under fast ice and **raised**
        in the ablation zone: a transport, not a source/sink.

        Two modes, selected by whether stratigraphy is on:

        - **Bulk bed** (no stratigraphy): the abraded volume ``Vtot`` is
          redistributed to the ablation cells weighted by the meltwater rate
          (``iceMeltL``) as a bed-to-bed transport, so ``Σ deposited == Vtot``
          and the net bed-volume change is zero by construction.
        - **Stratigraphic** (:meth:`_glacialTillStrata`): the till is removed
          from the layers it was abraded from and re-deposited as a fresh
          moraine layer in the ablation zone — split into the coarse/fine
          lithology fractions under dual lithology. Conservation is on the
          *solid* phase (per fraction), so the bed bulks up by the porosity
          contrast (uncompacted till vs compacted source rock).

        Both modes are guarded by a dedicated till-conservation test (the
        dual-lithology lesson: a volume the total budget can't see needs its
        own guard).

        Active only with ``till.on`` and ``Kg > 0``; a no-op otherwise.

        **Deposition distribution.** By default the till is spread across the
        whole ablation zone weighted by the meltwater rate. With
        ``till.route: True`` it is instead **routed down the ice-surface flow
        network** and melts out toward each catchment's terminus
        (:meth:`_routeTill`) — appropriate for high-resolution regional runs
        where individual glacier catchments and termini are resolved. Both
        produce a per-cell deposition weight (summing to one) consumed
        identically by the bulk and stratigraphic paths, so conservation is
        unchanged.
        """
        if not (self.iceOn and self.ice_till_on) or (
            self.ice_Kg <= 0.0 and self.ice_Kl <= 0.0
        ):
            return

        ti = process_time()
        ub = self.iceUbL.getArray()
        # Vertical abrasion under sliding ice + lateral wall abrasion (U-shaping);
        # both are bed-lowering routed into the same conserved till -> moraine.
        Eg = self.ice_Kg * np.power(np.maximum(ub, 0.0), self.ice_abr_l)  # m/yr
        Eg = Eg + self._glacialLateralErosion()
        Eg[self.seaID] = 0.0
        Vero = Eg * self.dt * self.larea                  # abraded volume (m^3) per cell
        owned = self.inIDs == 1
        Vtot = MPI.COMM_WORLD.allreduce(float(np.sum(Vero[owned])), op=MPI.SUM)

        # No abrasion -> nothing happens (skipping the erosion too keeps it
        # conservative: no orphaned till).
        if Vtot <= 0.0:
            return

        # Per-cell deposition weight (local array, Σ over owned nodes = 1):
        # melt-weighted spreading, or catchment-routed melt-out to the termini.
        if self.ice_till_route:
            dep_w = self._routeTill(Vero, Vtot, owned)
        else:
            melt = self.iceMeltL.getArray()               # ablation volume (m^3/yr)
            Wtot = MPI.COMM_WORLD.allreduce(
                float(np.sum(melt[owned])), op=MPI.SUM
            )
            # No ablation zone to receive the melt-spread till -> no-op.
            if Wtot <= 0.0:
                return
            dep_w = melt / Wtot

        # Routing can leave nothing depositable (e.g. no resolved outlet); guard.
        Wdep = MPI.COMM_WORLD.allreduce(float(np.sum(dep_w[owned])), op=MPI.SUM)
        if Wdep <= 0.0:
            return

        dz_ero = -Eg * self.dt                            # bed lowering (m, ≤0)

        if self.stratNb > 0:
            # Stratigraphy on: route the till through the stratigraphic pile so
            # the abraded material leaves the layers it came from and the
            # moraine is recorded as a fresh deposit (split coarse/fine under
            # dual lithology). See _glacialTillStrata.
            self._glacialTillStrata(dz_ero, dep_w, owned)
        else:
            # Bulk bed transport: erode at abrasion cells, deposit the till per
            # the deposition weight. Σ(dz·area) == 0 (rock moved, not created).
            dz_dep = (Vtot * dep_w) / self.larea
            dz = dz_ero + dz_dep
            self.tmpL.setArray(dz)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.cumED.axpy(1.0, self.tmp)
            self.hGlobal.axpy(1.0, self.tmp)
            self.dm.globalToLocal(self.cumED, self.cumEDLocal)
            self.dm.globalToLocal(self.hGlobal, self.hLocal)
            self._tillEroded += Vtot
            self._tillDeposited += MPI.COMM_WORLD.allreduce(
                float(np.sum((dz_dep * self.larea)[owned])), op=MPI.SUM
            )

        if MPIrank == 0 and self.verbose:
            print(
                "Glacial Till (%0.02f seconds)" % (process_time() - ti),
                flush=True,
            )

        return

    def _glacialTillStrata(self, dz_ero, dep_w, owned):
        r"""
        Stratigraphic / dual-lithology coupling for glacial till.

        The abraded bed lowering ``dz_ero`` (m, ≤0) is removed from the
        stratigraphic pile with :meth:`erodeStrat`, which returns the
        uncompacted solid removed per fraction (``thCoarse`` / ``thFine``,
        m/yr). That solid is the till: its total uncompacted volume is
        redistributed to the ablation cells weighted by the meltwater rate and
        laid down as a fresh moraine layer with :meth:`deposeStrat`, carrying
        the **abraded fine fraction** (a single, ice-mixed composition).

        Solid mass is conserved per fraction by construction: the fine volume
        deposited equals the fine volume eroded (so the dual-lithology
        ``_fineEroded`` / ``_fineDeposited`` budget stays balanced), and the
        coarse volume likewise. The bed bulks up by the porosity contrast
        between the uncompacted till and the compacted source rock.

        ``erodeStrat`` / ``deposeStrat`` overwrite the transient routing arrays
        (``thCoarse`` / ``thFine`` / ``depoFineFrac``) that the *fluvial*
        ``sedChange`` consumes later this step, so they are saved and restored
        around the till calls.
        """
        # Preserve the fluvial routing state (set by the eroder's erodeStrat,
        # consumed by the later sedChange/_getSedFlux). glacialTill runs
        # between the two.
        saved_thCoarse = self.thCoarse.copy() if hasattr(self, "thCoarse") else None
        saved_thFine = (
            self.thFine.copy()
            if self.stratLith and hasattr(self, "thFine")
            else None
        )
        saved_depoFineFrac = self.depoFineFrac.copy()

        # --- Erosion: remove the abraded bulk from the stratigraphic pile. ---
        self.tmpL.setArray(dz_ero)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.erodeStrat()                                 # sets thCoarse/thFine
        # Apply the bed lowering to the elevation / cumulative-erosion fields.
        self.cumED.axpy(1.0, self.tmp)
        self.hGlobal.axpy(1.0, self.tmp)

        # Uncompacted solid abraded (m^3): the till volume to redeposit.
        coarseV = MPI.COMM_WORLD.allreduce(
            float(np.sum((self.thCoarse * self.dt * self.larea)[owned])), op=MPI.SUM
        )
        if self.stratLith:
            fineV = MPI.COMM_WORLD.allreduce(
                float(np.sum((self.thFine * self.dt * self.larea)[owned])), op=MPI.SUM
            )
        else:
            fineV = 0.0
        Vsolid = coarseV + fineV
        if Vsolid <= 0.0:
            # Nothing actually came out of the pile (e.g. emptied to bedrock
            # floor): no till to deposit. Restore and return.
            self.thCoarse = saved_thCoarse
            if saved_thFine is not None:
                self.thFine = saved_thFine
            self.depoFineFrac = saved_depoFineFrac
            self.dm.globalToLocal(self.cumED, self.cumEDLocal)
            self.dm.globalToLocal(self.hGlobal, self.hLocal)
            self._tillEroded += 0.0
            return

        # Ice-mixed till composition (single fine fraction for the moraine).
        ffrac = fineV / Vsolid if self.stratLith else 0.0

        # --- Deposition: lay the till down per the deposition weight. ---
        dz_dep = (Vsolid * dep_w) / self.larea            # bulk till thickness (m)
        self.depoFineFrac = np.zeros(self.lpoints, dtype=np.float64)
        self.depoFineFrac[dz_dep > 0.0] = ffrac
        # Snapshot the current layer so the bed/cumED update uses the thickness
        # deposeStrat ACTUALLY laid down (it floors sub-1e-4 m deposits), keeping
        # elevation, the stratigraphic pile and the diagnostics consistent.
        h_before = self.stratH[:, self.stratStep].copy()
        self.tmpL.setArray(dz_dep)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.deposeStrat()                                # adds to stratH/stratHf
        dep = self.stratH[:, self.stratStep] - h_before   # actual deposit (m)
        self.tmpL.setArray(dep)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.cumED.axpy(1.0, self.tmp)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        # Diagnostics: solid abraded vs solid actually re-deposited (the two
        # match to the deposeStrat thickness floor; any sub-floor remainder is
        # dropped, as for every other deposit in the model).
        self._tillEroded += Vsolid
        self._tillDeposited += MPI.COMM_WORLD.allreduce(
            float(np.sum((dep * self.larea)[owned])), op=MPI.SUM
        )

        # Restore the fluvial routing state for sedChange.
        if saved_thCoarse is not None:
            self.thCoarse = saved_thCoarse
        if saved_thFine is not None:
            self.thFine = saved_thFine
        self.depoFineFrac = saved_depoFineFrac

        return

    def _glacialMeltwater(self, zbed, mdot):
        r"""
        Glacial meltwater volume (m^3/yr) delivered to the rivers.

        Two models, selected by ``ice.melt_conserve``:

        - **Discharge-conserving** (default): the accumulation (the water that
          fell as ice above the ELA, :math:`\dot{m}^+\,A`) is routed down the
          ice-surface flow network and **released as meltwater where the ice
          melts out** (fraction :math:`f`, forced to 1 at the margin/terminus),
          so :math:`\sum \mathrm{melt} = \sum \mathrm{accumulation}`. Over goSPL's
          long timesteps a land-terminating glacier is ~steady, so all the ice
          that accumulates leaves as meltwater at the snout — this closes the
          glacial water budget (the precipitation removed from runoff above the
          ELA returns downstream). Transport-with-loss on the ice network, same
          MPI-correct flow-matrix / KSP machinery as ``_routeTill``.
        - **Precip-scaled** (``melt_conserve: False``): the local ablation rate
          :math:`\max(-\dot{m},0)\,A` where ice is present — cheaper, but it
          loses water (the ablation generally under-returns the accumulation).
        """
        H = self.iceHL.getArray()
        if not self.ice_melt_conserve:
            return np.maximum(-mdot, 0.0) * self.larea * (H > 1.0e-2)

        # Halo-sync the thickness (ghosts needed for receivers / margin test).
        self.dm.localToGlobal(self.iceHL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        H = self.tmpL.getArray().copy()
        s = zbed + H                                       # ice surface
        rcv, _, wght = mfdreceivers(1, self.flowExp, s, self.sealevel, self.gid)
        rcv0 = rcv[:, 0].astype(petsc4py.PETSc.IntType)
        w0 = wght[:, 0].copy()

        thr = 1.0e-2
        ice = H > thr
        A = np.maximum(mdot, 0.0) * self.larea             # accumulation volume (m^3/yr)
        # Melt-out fraction: ablation share of the ice column this step, =1 off-ice
        # and at the margin (receiver ice-free) so all remaining ice melts out and
        # the transport is exactly conservative (Σ melt = Σ accumulation).
        abl = np.maximum(-mdot, 0.0) * self.dt
        f = np.clip(abl / np.maximum(H, thr), 0.0, 1.0)
        f[~ice] = 1.0
        f[ice & (H[rcv0] <= thr)] = 1.0

        outw = w0 * (1.0 - f)
        outw[self.ghostIDs] = 0.0
        if self.flatModel:
            outw[self.idBorders] = 0.0

        mat = self.iMat.copy()
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        nodes = indptr[:-1]
        data = -outw
        data[rcv0 == nodes] = 0.0
        tmpMat = self._matrix_build()
        tmpMat.assemblyBegin()
        tmpMat.setValuesLocalCSR(indptr, rcv0, data)
        tmpMat.assemblyEnd()
        mat.axpy(1.0, tmpMat)
        tmpMat.destroy()
        mat.transpose()

        self.tmpL.setArray(A)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self._solve_KSP(True, mat, self.tmp, self.tmp1)
        mat.destroy()
        self.dm.globalToLocal(self.tmp1, self.tmpL)
        L = self.tmpL.getArray().copy()                    # ice flux through each cell
        melt = f * L                                       # melt-out (m^3/yr)
        melt[melt < 0.0] = 0.0

        # The transport conserves analytically; the KSP solve is only accurate to
        # its tolerance, so renormalise so the released meltwater equals the
        # accumulation exactly (closing the water budget bit-for-bit).
        owned = self.inIDs == 1
        Mtot = MPI.COMM_WORLD.allreduce(float(np.sum(melt[owned])), op=MPI.SUM)
        Atot = MPI.COMM_WORLD.allreduce(float(np.sum(A[owned])), op=MPI.SUM)
        if Mtot > 0.0:
            melt *= Atot / Mtot
        return melt

    def _routeTill(self, Vero, Vtot, owned):
        r"""
        Catchment-aware till routing on the ice-surface flow network
        (``till.route: True``).

        The abraded till is transported down the **ice surface**
        :math:`s = z_\mathrm{bed} + H` along steepest descent (the same
        direction the ice flux follows) and **melts out** progressively toward
        each catchment's terminus, building moraine where the ice actually ends
        rather than smearing it across the whole ablation zone. This matters at
        high (sub-km) resolution where individual glacier catchments and termini
        are resolved.

        It is a transport-with-loss on the ice network, solved with the same
        MPI-correct flow-matrix / KSP machinery as the river accumulation:

        .. math::
            L_i = A_i + \sum_{u \to i} w_{ui}\,(1 - f_u)\,L_u,
            \qquad D_i = f_i\,L_i

        where :math:`A_i` is the local abraded volume, :math:`L_i` the till
        load passing through cell :math:`i`, and :math:`f_i \in [0,1]` the
        melt-out fraction — the share of the local ice column lost to ablation
        this step, :math:`f_i = \min(1, \dot{a}_i\,\Delta t / H_i)`. At the ice
        margin (a cell whose steepest-descent receiver is ice-free) :math:`f`
        is forced to 1 so no till leaks onto bare ground, which makes the
        transport exactly conservative: :math:`\sum_i D_i = \sum_i A_i`.

        :arg Vero: per-cell abraded volume (m^3, local array).
        :arg Vtot: total abraded volume over owned nodes (m^3).
        :arg owned: boolean mask of owned (non-ghost) nodes.

        :return: per-cell deposition weight (local array, Σ over owned = 1).
        """
        # Halo-synced ice thickness and meltwater (ghost values are needed for
        # the steepest-descent receivers and the melt-out fraction).
        self.dm.localToGlobal(self.iceHL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        H = self.tmpL.getArray().copy()
        self.dm.localToGlobal(self.iceMeltL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        melt = self.tmpL.getArray().copy()

        zbed = self.hLocal.getArray()
        s = zbed + H                                       # ice surface

        # Steepest-descent (single-flow-direction) receivers on the ice surface.
        rcv, _, wght = mfdreceivers(1, self.flowExp, s, self.sealevel, self.gid)
        rcv0 = rcv[:, 0].astype(petsc4py.PETSc.IntType)
        w0 = wght[:, 0].copy()

        thr = 1.0e-2                                       # ice-presence threshold (m)
        ice = H > thr
        # Melt-out fraction f = (ablation thickness this step) / H, clamped.
        meltThick = (melt * self.dt) / np.maximum(self.larea, 1.0e-12)
        f = np.clip(meltThick / np.maximum(H, thr), 0.0, 1.0)
        # Deposit fully (no downstream carry) off-ice and at the ice margin
        # (receiver ice-free), so nothing leaks onto bare ground -> conservative.
        f[~ice] = 1.0
        f[ice & (H[rcv0] <= thr)] = 1.0

        # Downstream-carried weight (1 - f); zero on ghosts/borders.
        outw = w0 * (1.0 - f)
        outw[self.ghostIDs] = 0.0
        if self.flatModel:
            outw[self.idBorders] = 0.0

        # Accumulation matrix (I - W^T) on the ice network, then L = solve(·, A).
        mat = self.iMat.copy()
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        nodes = indptr[:-1]
        data = -outw
        data[rcv0 == nodes] = 0.0
        tmpMat = self._matrix_build()
        tmpMat.assemblyBegin()
        tmpMat.setValuesLocalCSR(indptr, rcv0, data)
        tmpMat.assemblyEnd()
        mat.axpy(1.0, tmpMat)
        tmpMat.destroy()
        mat.transpose()

        self.tmpL.setArray(Vero)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self._solve_KSP(True, mat, self.tmp, self.tmp1)
        mat.destroy()
        self.dm.globalToLocal(self.tmp1, self.tmpL)
        L = self.tmpL.getArray().copy()

        # Deposited volume per cell. The transport conserves analytically, but
        # the KSP solve is only accurate to its tolerance, so normalise by the
        # ACTUAL routed total (not Vtot) to make the weight sum to one exactly —
        # mass conservation then matches the melt-weighted path bit-for-bit.
        # The analytical load is non-negative; clamp away sub-tolerance KSP
        # residual noise so the deposition weight can never go slightly negative.
        dep_vol = np.maximum(f * L, 0.0)
        total = MPI.COMM_WORLD.allreduce(
            float(np.sum(dep_vol[owned])), op=MPI.SUM
        )
        if total <= 0.0:
            return np.zeros(self.lpoints, dtype=np.float64)
        return dep_vol / total
