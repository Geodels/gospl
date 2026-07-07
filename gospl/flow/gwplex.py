import os
import petsc4py
import numpy as np

from time import process_time
from mpi4py import MPI

from gospl.tools.constants import ICE_COVER_MIN

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import jacobiancoeff

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()


class GWMesh(object):
    r"""
    Water table (groundwater) + generic duricrust — opt-in near-surface hydrology
    and chemical armoring. See ``docs/DESIGN_WATERTABLE_DURICRUST.md`` and the
    technical guide ``docs/tech_guide/groundwater.rst``.

    Enabled by the YAML ``groundwater:`` block (``self.gwOn`` / ``self.duriOn``);
    when absent, every path here is a no-op and goSPL is **byte-identical** to a
    run without it. The whole update runs once per step in ``updateGroundwater``
    (after flow accumulation, before erosion), so each process reads the previous
    process's state — a slow, stable, explicit/sequential coupling.

    **Water table.** An implicit (backward-Euler) Dupuit-Boussinesq head solve on
    the DMPlex (``_solveHead``): ``(I + (Δt/S)·L(T))·h = h_old + (Δt/S)·R`` with
    transmissivity ``T = Ksat·max(h − z_bed, b_min)``, a seepage free boundary
    (``h ≤ z`` at sea / lakes / open outlets) and a dry-aquifer floor
    (``h ≥ z_bed``). ``z_bed`` is prescribed (``aquifer_base`` scalar/map) or tied
    to the bedrock (``aquifer_base: from_soil`` ⇒ ``lHbed − bedrock_depth``).
    Recharge ``R = f_infil·max(0, rain − evap)`` is zeroed under water and ice.
    Opt-in baseflow closure returns the seepage discharge to the rivers
    (``conserve_baseflow``).

    **Duricrust.** A generic capillary-fringe crust (``_updateDuricrust``) grows
    where the water-table depth sits in the fringe band (``Φ`` Gaussian on
    ``wt = z − h``), fed by a weathering supply ``Ψ`` (climate proxy, an explicit
    Maher-Chamberlain rate, or the soil-production congruency; regolith-limited
    when soil is tracked). Its induration ``duriF ∈ [0,1]`` **armors erodibility**
    through the single ``_surfaceArmoringK`` hook (``1 − armor_max·duriF``), which
    slows erosion of crusted cells and drives relief inversion. When stratigraphy
    is on, the induration is archived per layer (``stratDuri``) so buried crusts
    are preserved and **re-arm on exhumation** (stacked-duricrust / cratonic
    laterite behaviour).

    **State & lifecycle.** Persistent Vecs (``headL/headG``, ``duriHL/duriHG``,
    ``rechargeL``, ``baseflowL``) and the cached elliptic solver (``_ksp_gw``,
    fgmres + hypre BoomerAMG) are registered in ``destroy_DMPlex``; ``head`` and
    ``duriH`` are model memory written to / restored from the output HDF5 on
    restart. Compatible with dual lithology and provenance (independent per-layer
    fields; armoring composes at the shared K hook).
    """

    def __init__(self, *args, **kwargs):
        """
        Initialise the groundwater / duricrust state. Allocated only when
        ``self.gwOn`` (set by ``inputparser._readGroundwater``); otherwise the
        cached-solver handles are still set to ``None`` so ``destroy_DMPlex`` and
        any ``getattr`` guards are safe.
        """

        # Cached elliptic head solver + operator (built lazily on first solve).
        # Set unconditionally so destroy_DMPlex / guards never hit a missing attr.
        self._gwMat = None
        self._ksp_gw = None
        # Cached Level-B solute-transport solver + operator (built lazily; G1+).
        self._soluteMat = None
        self._ksp_solute = None

        if getattr(self, "gwOn", False):
            # --- PETSc state (persistent, halo-synced; in destroy_DMPlex) ---
            # Water-table head h (elevation of the saturated surface, m).
            self.headL = self.hLocal.duplicate()
            self.headG = self.hGlobal.duplicate()
            # Duricrust thickness duriH (m).
            self.duriHL = self.hLocal.duplicate()
            self.duriHG = self.hGlobal.duplicate()
            # Net recharge R this step (m/yr) — diagnostic/output.
            self.rechargeL = self.hLocal.duplicate()
            # Seepage return to rivers (m^3/yr) — only used when conserve_baseflow.
            self.baseflowL = self.hLocal.duplicate()

            # Resolve a per-vertex aquifer_base map (done here — needs locIDs).
            # `[file, key]` -> file + ".npz" subset to the local partition; the
            # z_bed depth below the surface (m) then varies in space (regolith /
            # weathering-front depth). Scalar / 'from_soil' left untouched.
            basemap = getattr(self, "_gwAquiferBaseMap", None)
            if basemap is not None:
                data = np.load(basemap[0] + ".npz")
                self.gwAquiferBase = data[basemap[1]][self.locIDs].astype(
                    np.float64
                )

            # Seed the head at the aquifer BASE (a dry start) so recharge fills
            # it UP to the steady water table. Seeding at the surface is wrong:
            # recharge would push h above z, the seepage clip would pin it there,
            # and it could never drain below the surface (spurious full
            # saturation). For a scalar / map aquifer_base, z_bed = z − base;
            # for `from_soil` (Phase 5, a string) fall back to the surface.
            base = self.gwAquiferBase
            if isinstance(base, str):
                self.hGlobal.copy(result=self.headG)
                self.hLocal.copy(result=self.headL)
            else:
                self.headL.setArray(self.hLocal.getArray() - base)
                self.dm.localToGlobal(self.headL, self.headG)
                self.dm.globalToLocal(self.headG, self.headL)
            self.duriHL.set(0.0)
            self.duriHG.set(0.0)
            self.rechargeL.set(0.0)
            self.baseflowL.set(0.0)

            # --- numpy state (rank-local, no halo) ---
            self.wtDepth = np.zeros(self.lpoints, dtype=np.float64)     # z - h (m)
            self.duriF = np.zeros(self.lpoints, dtype=np.float64)       # induration 0..1
            self.duriKarmor = np.ones(self.lpoints, dtype=np.float64)   # K multiplier (<=1)
            self.gwSeepIDs = np.zeros(0, dtype=np.int64)                # Dirichlet seepage nodes
            # Surface at the previous groundwater update — the duricrust breakdown
            # term reads the per-step incision (z_last − z > 0 = surface lowered).
            self.gwZlast = self.hLocal.getArray().copy()
            # Signed across-bed groundwater flux (m³/yr, >0 = into the surface),
            # for the opt-in lake ↔ aquifer volume coupling. 0 until the first
            # head solve, so the cascade's first read is a no-op.
            self.gwLakeFlux = np.zeros(self.lpoints, dtype=np.float64)
            # Running total of net groundwater volume exchanged with lakes (m³,
            # owned nodes) — diagnostic for the lake-exchange budget.
            self.gwLakeInflow = 0.0
            # Cached local temperature (K) for the optional Arrhenius weathering
            # term (loaded lazily from the soil tempMap when weather_Ea > 0), and
            # a lazily-resolved per-vertex `weatherability` map (rate mode).
            self._gwTempK = None
            self._duriWeatherArr = None

            # Resolve a per-vertex infiltration map (done here — needs locIDs,
            # unavailable at parse time). `[file, key]` -> file + ".npz", subset
            # to the local partition. Scalar `gwInfiltration` is left untouched.
            infilmap = getattr(self, "_gwInfilMap", None)
            if infilmap is not None:
                data = np.load(infilmap[0] + ".npz")
                self.gwInfiltration = data[infilmap[1]][self.locIDs].astype(
                    np.float64
                )

            # --- Level-B geochemistry state (opt-in; DESIGN_WATERTABLE_GEOCHEM.md
            # G0). Allocated only when `gwGeochemOn`; nothing is solved yet (G1+).
            # The solute is an (lpoints, n_species) array — single-tracer when
            # n_species=1, multi-tracer otherwise — sharing one transport template.
            if getattr(self, "gwGeochemOn", False):
                nsp = int(self.gwNspecies)
                if getattr(self, "duriOn", False) and MPIrank == 0 and self.verbose:
                    print(
                        "[gw] geochem on: the duricrust source is the transported "
                        "precipitated solute — the per-species `precip_rate` (a "
                        "1/yr rate) governs crust growth; the Level-A `form_rate` "
                        "is NOT used. Scale `precip_rate` down for long timesteps.",
                        flush=True,
                    )
                # Per-species parameters (parser lists -> numpy, len n_species).
                # `gwGeoWeather` stays a RAW list (each entry a scalar, a
                # per-vertex `[file, key]` map, or a table used with a lithology
                # label) — resolved lazily into `_gwGeoWeatherArr` on first use.
                self._gwGeoWeatherArr = None
                self.gwGeoRiverDecay = np.asarray(self.gwGeoRiverDecay, dtype=np.float64)
                self.gwGeoCsat = np.asarray(self.gwGeoCsat, dtype=np.float64)
                self.gwGeoPrecip = np.asarray(self.gwGeoPrecip, dtype=np.float64)
                self.gwGeoVsolid = np.asarray(self.gwGeoVsolid, dtype=np.float64)
                # Groundwater solute concentration and the dissolvable source pool
                # (rank-local, per tracer). The pool is seeded per species as
                # `source_pool · cell area` (weatherable rock mass, default 1e6);
                # dissolution debits it, so a pool too small for the run's
                # weathering rate exhausts and the tracer goes inert. Set
                # `source_pool` large (a thick saprolite source) to keep weathering
                # rate-limited over long, fast-weathering runs. `_solutePoolWarned`
                # gates the one-time exhaustion warning.
                self.gwGeoSourcePool = np.asarray(self.gwGeoSourcePool, dtype=np.float64)
                self.gwSolute = np.zeros((self.lpoints, nsp), dtype=np.float64)
                self.gwSourcePool = self.larea[:, None] * self.gwGeoSourcePool[None, :]
                self._solutePoolWarned = False
                # Cumulative mass budget per tracer (m³-equiv, owned nodes) for the
                # conservation guard / diagnostics: dissolved = precipitated +
                # exported (ocean) + currently in solution.
                self.gwDissolved = np.zeros(nsp, dtype=np.float64)
                self.gwPrecip = np.zeros(nsp, dtype=np.float64)
                self.gwOceanFlux = np.zeros(nsp, dtype=np.float64)
                # Per-node baseflow-carried solute export (m³/yr, summed over
                # tracers) — the spatial output field (G3) — and its per-species
                # split (`soluteflux_<name>` output + the river-routing source).
                self.gwSoluteFlux = np.zeros(self.lpoints, dtype=np.float64)
                self.gwSoluteFluxSp = np.zeros((self.lpoints, nsp), dtype=np.float64)
                # Per-node cumulative crust precipitated by each tracer, and the
                # dominant crust-forming tracer (G4 typing: -1 = no crust).
                self.gwCrustBySpecies = np.zeros((self.lpoints, nsp), dtype=np.float64)
                self.gwCrustType = np.full(self.lpoints, -1, dtype=np.int32)
                # Solute-source provenance (G5, only with in-model provenance on):
                # cumulative crust attributed to each source-rock class (where the
                # solute dissolved), and the dominant source class per node.
                if getattr(self, "provOn", False):
                    self.gwCrustProv = np.zeros(
                        (self.lpoints, int(self.provNb)), dtype=np.float64
                    )
                    self.gwCrustSource = np.full(self.lpoints, -1, dtype=np.int32)
                # Scratch Vec pair for the per-species transport solve (G1+).
                self.soluteL = self.hLocal.duplicate()
                self.soluteG = self.hGlobal.duplicate()
                self.soluteL.set(0.0)
                self.soluteG.set(0.0)
                # ext 2: river dissolved-load routing — accumulated solute flux
                # (m³/yr) down the surface network + its ocean-delivered total.
                if getattr(self, "gwRiverLoad", False):
                    self.riverSoluteG = self.hGlobal.duplicate()
                    self.riverSoluteL = self.hLocal.duplicate()
                    self.riverSoluteG.set(0.0)
                    self.riverSoluteL.set(0.0)
                    self.riverSolute = np.zeros(self.lpoints, dtype=np.float64)
                    self.riverSoluteToOcean = 0.0
                    # Per-species routing: routed river load by tracer (the total
                    # fields above are Σ over k; the per-species seepage source
                    # `gwSoluteFluxSp` is allocated with the geochem state above).
                    self.riverSoluteSp = np.zeros((self.lpoints, nsp), dtype=np.float64)
                    self.riverSoluteToOceanSp = np.zeros(nsp, dtype=np.float64)
                    self.riverSoluteLost = np.zeros(nsp, dtype=np.float64)  # in-transit
                    # Marine coupling: cumulative delivered inventory per tracer
                    # (m³, integrated over time) + per-node coastal input this step.
                    if getattr(self, "gwMarineCoupling", False):
                        self.marineSolute = np.zeros(nsp, dtype=np.float64)
                        self.marineSoluteInput = np.zeros(self.lpoints, dtype=np.float64)

        return

    def updateGroundwater(self):
        """
        Per-step groundwater / duricrust update.

        1. **Recharge** ``R = f_infil * max(0, rain - evap)`` (m/yr), held at 0
           under standing water (marine ``seaID`` or a ponded continental lake
           ``pitIDs>-1 & lFill>hl``) and under ice, where rain does not infiltrate.
           Stored in ``self.rechargeL`` for output.
        2. **Head solve** (``_solveHead``): the implicit Dupuit-Boussinesq water
           table from that recharge; sets ``self.wtDepth = z - h``.
        3. **Duricrust** (``_updateDuricrust``, only when ``duriOn``): evolve the
           capillary-fringe crust ``duriH`` and the induration ``duriF`` / armor
           multiplier ``duriKarmor`` from the new water-table depth.
        4. **Geochemistry** (``_updateSolute``, only when ``gwGeochemOn`` — Level
           B): per tracer, dissolve → transport (down ``q = -T∇h``) → precipitate
           at the fringe (feeding ``duriH``, replacing the Level-A supply) →
           export; a closed mass budget with per-species crust typing, optional
           solute-source provenance and a dissolved ocean flux.
        5. **River routing** (``_routeRiverSolute``, only when ``gwRiverLoad``):
           route the exported solute down the surface drainage network to the
           coast (per species; optional in-transit loss + marine coupling).
        6. **Stratigraphic archive** (``_recordInduration``, when ``duriOn`` and
           stratigraphy is on): record the induration degree and, with geochem,
           the crust's dominant species / source per layer (exhumation re-arm +
           formation write-down).

        No-op when ``gwOn`` is off. The head, solute-transport and river-routing
        solves are collective (KSP); recharge, duricrust and the archive are
        purely rank-local (per-node).
        """
        if not getattr(self, "gwOn", False):
            return

        t0 = process_time()
        rain = self.rainVal
        evap = getattr(self, "evapVal", None)
        net = np.maximum(0.0, rain if evap is None else (rain - evap))

        # Effective infiltration fraction, optionally modulated (opt-in) by:
        #  - surface lithology (coarse infiltrates more than fine) via the exposed
        #    coarse fraction and `fine_infil_factor` (dual lithology only);
        #  - terrain slope (steeper ⇒ less infiltration, more runoff) via
        #    f/(1 + slope/infil_slope_ref). See DESIGN §3.
        finf = self.gwInfiltration
        if self.gwFineInfilFactor != 1.0 and getattr(self, "stratLith", False):
            fc = self._surfaceComposition()                # exposed coarse fraction
            finf = finf * (fc + (1.0 - fc) * self.gwFineInfilFactor)
        if self.gwInfilSlopeRef > 0.0 and getattr(self, "rcvID", None) is not None:
            finf = finf / (1.0 + self._surfaceSlope() / self.gwInfilSlopeRef)
        R = finf * net

        # No rain-recharge where the surface is not subaerial land:
        #  - standing water (marine `seaID`, or a ponded continental lake) — the
        #    head is pinned to the surface there;
        #  - ice-covered land — precipitation falls as snow/ice and does not
        #    infiltrate the ground (parallels the soil ice-freeze gate). Subglacial
        #    meltwater recharge is added back below (opt-in).
        sub = np.zeros(self.lpoints, dtype=bool)
        sub[self.seaID] = True
        pitIDs = getattr(self, "pitIDs", None)
        lFill = getattr(self, "lFill", None)
        if pitIDs is not None and lFill is not None:
            sub |= (pitIDs > -1) & (lFill > self.hLocal.getArray())
        if getattr(self, "iceOn", False):
            iceHL = getattr(self, "iceHL", None)
            if iceHL is not None:
                sub |= iceHL.getArray() > ICE_COVER_MIN
        R = np.where(sub, 0.0, R)

        # Subglacial-meltwater recharge (opt-in): a fraction of the glacial
        # meltwater (`iceMeltRiverL`, m³/yr) infiltrates the aquifer where it melts
        # out — the one recharge path allowed under ice (the ice gate above zeroes
        # the rain path). Converted to a rate (m/yr) by the cell area; never on the
        # sea. See DESIGN_WATERTABLE_DURICRUST.md §3.
        if self.gwSubglacial > 0.0 and getattr(self, "iceOn", False):
            imr = getattr(self, "iceMeltRiverL", None)
            if imr is not None:
                Rsub = self.gwSubglacial * imr.getArray() / self.larea
                Rsub[self.seaID] = 0.0
                R = R + Rsub

        self.rechargeL.setArray(R)

        # Phase 2: solve the implicit water-table head from this recharge.
        self._solveHead()

        # Phase 3: evolve the capillary-fringe duricrust from the new water table
        # (breakdown/decay; the Level-A formation supply is off when geochem is on).
        if getattr(self, "duriOn", False):
            self._updateDuricrust()
        # Level-B (G2): dissolve → transport → precipitate the solute; the
        # precipitation feeds the crust `duriH` (replacing the Level-A supply).
        if getattr(self, "gwGeochemOn", False):
            self._updateSolute()
            # ext 2: route the exported solute down the rivers to the shoreline
            # (river dissolved load). Uses the flow matrix built earlier this step.
            self._routeRiverSolute()
        # Phase 6: sync the live crust with the per-layer stratigraphic induration
        # archive (exhumation re-arm + formation write-down) after both sources.
        if getattr(self, "duriOn", False):
            self._recordInduration()

        if MPIrank == 0 and self.verbose:
            print(
                "Update Groundwater Table (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )
        return

    def _makeGWKSP(self):
        """
        Cached KSP for the elliptic head solve: fgmres + **hypre BoomerAMG**
        (algebraic multigrid). The head operator ``I + (Δt/S)·L(T)`` is a stiff,
        high-condition 2-D elliptic operator; block-Jacobi/ILU does NOT converge
        on it (it stalls at ``DIVERGED_ITS`` and returns a garbage iterate that
        drifts to the surface) — an elliptic solve needs a multigrid PC. Same
        rationale as the flexure biharmonic wanting a strong solver. Options
        prefix ``gw_`` scopes any override (env ``-gw_pc_type ...``).
        """
        ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
        ksp.setType("fgmres")
        ksp.getPC().setType("hypre")
        ksp.setTolerances(rtol=1.0e-10, max_it=500)
        ksp.setInitialGuessNonzero(True)
        ksp.setOptionsPrefix("gw_")
        ksp.setFromOptions()
        return ksp

    def _solveHead(self):
        r"""
        Implicit (backward-Euler) Dupuit-Boussinesq water-table head solve
        (DESIGN_WATERTABLE_DURICRUST.md §2/§3). One elliptic solve per step:

        .. math::
            (I + \tfrac{\Delta t}{S}\,L(T))\,h = h_{old} + \tfrac{\Delta t}{S}\,R

        with transmissivity :math:`T = K_h\max(h - z_{bed}, b_{min})` (unconfined),
        :math:`L(T)=-\nabla\cdot(T\nabla)` the FV neg-Laplacian (from
        ``jacobiancoeff``, area-normalised — the same operator the marine Picard
        solver uses). Reduces to the steady elliptic limit as :math:`\Delta t/S`
        grows, so it is quasi-steady at century steps.

        - **Unconfined non-linearity** ``T(h)``: ``gwPicardIts`` Picard sweeps
          (lag ``T`` on the previous iterate; operator rebuilt per sweep).
        - **Seepage free boundary** ``h <= z``: Dirichlet ``h = z`` at the base
          seepage set (marine ``seaID`` + ponded lakes + open ``outletIDs``) via
          ``zeroRowsLocal``; after each Picard block, clip ``h = min(h, z)`` and
          add any newly-over-topped cell to the Dirichlet set, up to
          ``gwSeepagePasses`` passes (break on an ``Allreduce``'d new-seepage
          count — collective-safe). ``aquifer_base`` is prescribed (scalar/map)
          OR, with ``aquifer_base: from_soil``, tied to the bedrock elevation
          ``z_bed = lHbed − bedrock_depth`` (Phase 5 soil coupling; ``_gwZbed``).
        - **Baseflow closure** (opt-in ``conserve_baseflow``): the seepage-return
          discharge is accounted into ``self.baseflowL`` (``_baseflowClosure``).

        Scratch: uses ``self.tmp`` (global rhs); leaves it defined. Writes
        ``self.headL/headG`` and ``self.wtDepth = z - h``.
        """
        z = self.hLocal.getArray()
        zbed = self._gwZbed(z)                  # scalar/map, or lHbed − d_bedrock
        bmin = float(self.gwMinSatThick)
        gwdt = self.dt / float(self.gwSpecificYield)          # Δt/S
        R = self.rechargeL.getArray()
        inIDs = self.inIDs
        zeroKp = np.zeros(self.lpoints, dtype=np.float64)
        IntType = petsc4py.PETSc.IntType

        # Base seepage (Dirichlet h = z): standing water + open outlets.
        seep = np.zeros(self.lpoints, dtype=bool)
        seep[self.seaID] = True
        pitIDs = getattr(self, "pitIDs", None)
        lFill = getattr(self, "lFill", None)
        if pitIDs is not None and lFill is not None:
            seep |= (pitIDs > -1) & (lFill > z)
        outletIDs = getattr(self, "outletIDs", None)
        if outletIDs is not None and len(outletIDs) > 0:
            seep[outletIDs] = True

        if self._ksp_gw is None:
            self._ksp_gw = self._makeGWKSP()
        ksp = self._ksp_gw

        self.dm.globalToLocal(self.headG, self.headL)
        hloc = self.headL.getArray().copy()
        hold = hloc.copy()

        npass = 0
        seep_converged = False
        for _ in range(int(self.gwSeepagePasses)):
            npass += 1
            for _ in range(int(self.gwPicardIts)):
                T = self.gwKsat * np.maximum(hloc - zbed, bmin)   # transmissivity (m²/yr)
                coeffs = gwdt * jacobiancoeff(hloc, T, zeroKp)     # (Δt/S)·L(T)
                coeffs[:, 0] += 1.0                                # I + (Δt/S)·L(T)
                M = self._assembleDiffMatCSR(coeffs)
                owned_seep = np.where(seep & (inIDs == 1))[0].astype(IntType)
                M.zeroRowsLocal(owned_seep, diag=1.0)             # collective; h = rhs there

                rhs = hold + gwdt * R
                rhs[seep] = z[seep]                               # pin h = z at seepage
                self.headL.setArray(rhs)
                self.dm.localToGlobal(self.headL, self.tmp)       # rhs (global)
                self.headL.setArray(hloc)
                self.dm.localToGlobal(self.headL, self.headG)     # nonzero guess
                ksp.setOperators(M, M)
                ksp.solve(self.tmp, self.headG)
                M.destroy()
                self.dm.globalToLocal(self.headG, self.headL)
                hloc = self.headL.getArray().copy()

            # Free boundaries: seepage (upper, h <= surface — discover new
            # seepage cells) and the dry-aquifer floor (lower, h >= z_bed — the
            # head cannot drop below the impermeable base).
            over = (hloc > z) & (~seep)
            n_new = MPI.COMM_WORLD.allreduce(
                int((over & (inIDs == 1)).sum()), op=MPI.SUM
            )
            seep |= over
            hloc = np.clip(hloc, zbed, z)
            if n_new == 0:
                seep_converged = True
                break

        if MPIrank == 0 and self.verbose:
            # `reason > 0` = KSP converged; a non-converged elliptic solve (e.g.
            # a too-weak PC) silently corrupts the water table, so surface it.
            reason = int(ksp.getConvergedReason())
            print(
                "[gw] head solve: %d seepage pass(es)%s, last KSP reason %d%s"
                % (
                    npass,
                    "" if seep_converged else " (seepage set still growing)",
                    reason,
                    "" if reason > 0 else " -- WARNING: KSP did not converge",
                ),
                flush=True,
            )

        self.headL.setArray(hloc)
        self.dm.localToGlobal(self.headL, self.headG)
        # Water-table depth, BOUNDED by the aquifer thickness (z − z_bed): a water
        # table cannot sit below its aquifer base. Where the head drains below the
        # base (a dry aquifer — kept solvable by the `min_sat_thickness` floor on
        # T), `z − h` would otherwise report a spurious depth far exceeding the
        # aquifer (e.g. tens of m in a few-m regolith with `aquifer_base:
        # from_soil`), which also zeroes the capillary-fringe favourability Φ and
        # wrongly suppresses the duricrust there. Clamp to [0, z − z_bed].
        self.wtDepth = np.clip(z - hloc, 0.0, np.maximum(z - zbed, 0.0))

        # Phase 5: account the seepage-return (baseflow) discharge (opt-in).
        if getattr(self, "gwConserveBaseflow", False):
            self._baseflowClosure(hold, hloc, seep)
        # Opt-in lake ↔ aquifer volume coupling: the signed across-bed flux.
        if getattr(self, "gwLakeExchange", False):
            self._lakeExchangeFlux(hloc, zbed, bmin, zeroKp)
        return

    def _lakeExchangeFlux(self, hloc, zbed, bmin, zeroKp):
        r"""
        Signed groundwater flux exchanged with the surface at every node (m³/yr),
        stored in ``self.gwLakeFlux`` for the lake volume coupling (§15). It is the
        FV divergence of the lateral groundwater flow, ``∇·(T∇h)·A = −(L·h)·A``
        (``L`` = the area-normalised neg-Laplacian from ``jacobiancoeff``, so
        ``L·h = −∇·(T∇h)``): **positive = the aquifer discharges INTO the surface**
        (a gaining lake), **negative = the surface leaks INTO the aquifer** (a
        losing lake). One extra operator assemble + mat-vec, only when
        ``lake_exchange`` is on. Uses the *un-zeroed* operator (the seepage-Dirichlet
        rows would otherwise null the flux exactly at the lake nodes we need).
        """
        T = self.gwKsat * np.maximum(hloc - zbed, bmin)
        Lop = self._assembleDiffMatCSR(jacobiancoeff(hloc, T, zeroKp))  # pure L
        self.headL.setArray(hloc)
        self.dm.localToGlobal(self.headL, self.headG)
        Lop.mult(self.headG, self.tmp)                     # tmp = L·h (global)
        Lop.destroy()
        self.dm.globalToLocal(self.tmp, self.headL)
        Lh = self.headL.getArray()
        self.gwLakeFlux = -Lh * self.larea                 # >0 aquifer→surface
        # restore headG/headL to the head (tmp/headL were reused as scratch)
        self.headL.setArray(hloc)
        self.dm.localToGlobal(self.headL, self.headG)
        return

    def _gwZbed(self, z):
        r"""
        Aquifer-base elevation ``z_bed`` (m). Prescribed by ``aquifer_base``
        (scalar or per-vertex map, ``z_bed = z − aquifer_base``), OR — with
        ``aquifer_base: from_soil`` (Phase-5 soil coupling) — tied to the bedrock
        elevation ``z_bed = lHbed − bedrock_depth``, the *"permeable regolith over
        impermeable bedrock"* model. ``from_soil`` needs ``soilSPL`` (``lHbed``);
        it falls back to the surface (a self-consistent no-aquifer floor) with a
        one-time warning if soil is off.

        **Depositional basins (`from_soil` + stratigraphy).** Under Option-2.5
        ``lHbed = z − Lsoil`` is only the base of the thin weathering regolith;
        deposited sediment lives in the stratigraphy, not ``Lsoil``. So in a
        filled basin the porous sediment column *is* the aquifer and its base is
        the bottom of the (non-sentinel) stratigraphic pile, deeper than
        ``lHbed``. The base is therefore taken as the **deeper** (lower) of
        ``lHbed − bedrock_depth`` and ``z − Σ sediment thickness`` — so uplands
        (thin pile) keep the regolith base while basins deepen to the fill base.
        Rank-local.
        """
        base = self.gwAquiferBase
        if isinstance(base, str):                       # 'from_soil'
            lHbed = getattr(self, "lHbed", None)
            if lHbed is None:
                if MPIrank == 0 and self.verbose:
                    print(
                        "[gw] aquifer_base: from_soil needs soilSPL (lHbed) — "
                        "falling back to the surface as the aquifer base.",
                        flush=True,
                    )
                return z.copy()
            zbed = lHbed.getArray() - float(self.gwBedrockDepth)
            if self.stratNb > 0 and self.stratH is not None:
                lo = int(getattr(self, "bedrockLay", 0))    # skip bedrock sentinel
                top = self.stratStep + 1
                if top > lo:
                    sed = self.stratH[:, lo:top].sum(axis=1)  # porous fill thickness
                    zbed = np.minimum(zbed, z - sed)          # basin fill = aquifer
            return zbed
        return z - base

    def _syncHalo(self, arr):
        r"""
        Return ``arr`` with its halo/ghost entries overwritten by the **owning
        rank's** value, via a local→global→local round-trip (``localToGlobal``
        INSERT takes the owner's value; ``globalToLocal`` broadcasts it back).

        For output fields computed on OWNED nodes only (e.g. the seepage export
        ``seep_sink·c·A``, or the area-weighted baseflow), the local array's halo
        stays 0 / stale, which renders as a **seam at partition boundaries** in
        the per-rank output. Syncing removes the seam. Serial runs have no halo,
        so this is a no-op there. Rank-local vectors; the round-trip is collective.
        """
        self.tmpL.setArray(np.ascontiguousarray(arr, dtype=np.float64))
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        return self.tmpL.getArray().copy()

    def _subaerialMask(self):
        r"""
        Boolean per-node mask of **subaerial land** — not marine (``seaID``) and
        not under a ponded continental lake (``pitIDs>-1 & lFill>z``). The
        capillary-fringe duricrust forms, and Level-B solute precipitates, ONLY
        here — there is no capillary fringe under standing water. Rank-local.
        """
        z = self.hLocal.getArray()
        sub = np.zeros(self.lpoints, dtype=bool)
        sub[self.seaID] = True
        pitIDs = getattr(self, "pitIDs", None)
        lFill = getattr(self, "lFill", None)
        if pitIDs is not None and lFill is not None:
            sub |= (pitIDs > -1) & (lFill > z)
        return ~sub

    def _baseflowClosure(self, hold, hnew, seep):
        r"""
        Baseflow (seepage-return) accounting (DESIGN_WATERTABLE_DURICRUST.md §3
        step 7), opt-in via ``conserve_baseflow``. Over the quasi-steady step the
        net recharge that does not go into aquifer storage discharges back to the
        surface-water network at the seepage nodes:

        .. math:: Q_{seep} = \sum_{owned}(R\,A) - \sum_{owned} S\,(h-h_{old})\,A/\Delta t

        (m³/yr). The global budget is reduced across ranks (``Allreduce``) and
        distributed over the owned seepage nodes weighted by cell area, so
        ``Σ_owned baseflowL ≈ Σ_owned recharge`` in the steady limit (``ΔS→0``) —
        river discharge stays ``≈ rain − evap``. Stored in ``self.baseflowL``
        (diagnostic/output); re-injection into the surface-flow source (making
        rivers physically baseflow-fed, which redistributes erosion) is the next
        increment. Collective (two ``Allreduce``); no rank-local collective gate.
        """
        A = self.larea
        S = float(self.gwSpecificYield)
        owned = self.inIDs == 1
        R = self.rechargeL.getArray()
        # Net discharge rate = recharge volume − storage-change volume (m³/yr).
        rvol = np.where(owned, R * A, 0.0).sum()
        svol = np.where(owned, S * (hnew - hold) * A, 0.0).sum() / self.dt
        Qtot = MPI.COMM_WORLD.allreduce(float(rvol - svol), op=MPI.SUM)

        bf = np.zeros(self.lpoints, dtype=np.float64)
        # Distribute the return flow over the SUBAERIAL seepage nodes only — the
        # springs / lakes / land outlets that actually feed the rivers. Marine
        # (`seaID`) nodes are seepage BCs too, but baseflow injected there is lost
        # to the sea (it never reaches a river), which would break the water
        # conservation this closure exists to enforce (and paints "baseflow" over
        # the ocean in the output). So exclude the marine seepage nodes.
        marine = np.zeros(self.lpoints, dtype=bool)
        marine[self.seaID] = True
        seep_owned = seep & owned & (~marine)
        wsum = MPI.COMM_WORLD.allreduce(
            float(np.where(seep_owned, A, 0.0).sum()), op=MPI.SUM
        )
        if wsum > 0.0:
            bf[seep_owned] = Qtot * A[seep_owned] / wsum
        # `bf` is set on OWNED nodes only (0 on the halo). Sync the halo from the
        # owning rank so shared/ghost nodes carry the owner's value; otherwise the
        # local `baseflowL` has a 0 ring at partition boundaries — a visible seam
        # in the output. The re-injection is unaffected (localToGlobal INSERT
        # already takes the owner's value); this just fixes the halo/output.
        self.baseflowL.setArray(self._syncHalo(bf))
        return

    def _arrhenius(self):
        r"""
        Optional Arrhenius temperature factor for the weathering supply,
        ``exp(Ea/Rg · (1/T_ref − 1/T))`` (dimensionless, 1 at ``T = T_ref``).

        Returns ``1.0`` (no temperature dependence) when ``weather_Ea <= 0`` OR
        no temperature map is available — so duricrust runs soil-free by default.
        Reuses the soil ``tempMap`` (``self.tempFile``/``tempData``/``tempRef``)
        when present; loaded once and cached. Rank-local.
        """
        if self.duriWeatherEa <= 0.0:
            return 1.0
        if self._gwTempK is None:
            tempFile = getattr(self, "tempFile", None)
            if tempFile is None:
                if MPIrank == 0 and self.verbose:
                    print(
                        "[gw] duricrust weather_Ea > 0 but no soil tempMap — "
                        "Arrhenius term disabled (factor 1).",
                        flush=True,
                    )
                self.duriWeatherEa = 0.0          # don't retry every step
                return 1.0
            data = np.load(tempFile)
            self._gwTempK = data[self.tempData][self.locIDs] + 273.15
        Tref = getattr(self, "tempRef", 15.0) + 273.15
        return np.exp(self.duriWeatherEa * (1.0 / Tref - 1.0 / self._gwTempK) / 8.314)

    def _surfaceSourceClass(self):
        r"""
        Per-node **current surface lithology** class — the dominant provenance
        class of the **top non-empty stratigraphic layer** (from ``stratP``),
        falling back to the static bedrock ``source_class`` where a column has no
        sediment. As erosion exhumes deeper layers (or deposition buries the
        surface under transported sediment of a different provenance), the top
        layer — and hence this label — changes, so a ``surface_class``-driven
        weatherability tracks the rock actually exposed each step.

        Uses the same top-non-empty-layer scan as ``_recordInduration`` /
        ``_surfaceComposition``. Needs provenance (``stratP``); returns the plain
        ``source_class`` (or ``None``) when the per-layer class record is absent.
        Rank-local; no collective.
        """
        base = getattr(self, "source_class", None)
        stratP = getattr(self, "stratP", None)
        if base is None or stratP is None or self.stratNb == 0:
            return base
        top = self.stratStep + 1
        H = self.stratH[:, :top]
        rev = (H > 0)[:, ::-1]
        valid = rev.any(axis=1)                       # columns with any sediment
        label = base.copy()
        if valid.any():
            top_idx = (H.shape[1] - 1 - np.argmax(rev, axis=1))[valid]
            rows = np.arange(H.shape[0])[valid]
            label[rows] = stratP[rows, top_idx, :].argmax(axis=1)
        return label

    def _resolveGeoWeather(self):
        r"""
        Resolve the per-species weatherability to a list of **scalar or
        ``(lpoints,)`` array** entries (Level-B extension 1 — lithology →
        chemistry). Lets lithology control *which* species each region yields:

        - **(a) per-vertex map** — a species' ``weatherability: [file, key]`` is
          loaded and subset to the local partition (``self.locIDs``), exactly
          like ``duriWeatherability`` / ``gwInfiltration``.
        - **(b) table by a standalone lithology map** — a species'
          ``weatherability_by_class: [...]`` is gathered by a per-vertex integer
          ``lithology: [file, key]`` map (independent of provenance).
        - **(c) table by the provenance label** — the same
          ``weatherability_by_class`` gathered by a per-vertex class label:
          ``weatherability_from: source_class`` uses the **static bedrock**
          class, while ``weatherability_from: surface_class`` uses the
          **dynamic current surface** class (the top stratigraphic layer,
          recomputed each step, so weatherability tracks exhumation/burial). Both
          need no new input beyond provenance.
        - otherwise the scalar (the default; byte-identical to the old path).

        The lithology label for (b)/(c) is the ``lithology:`` map when given,
        else the provenance class (static ``source_class`` or dynamic
        ``surface_class``). Resolved lazily (first ``_updateSolute``) so it runs
        after provenance seeds ``source_class``; the dynamic ``surface_class``
        form is re-resolved every step. Rank-local, partition-exact, no
        collective.
        """
        label = None
        lithomap = getattr(self, "_gwLithoMap", None)
        wfrom = getattr(self, "gwWeatherFrom", None)
        if lithomap is not None:                          # (b) standalone map
            d = np.load(lithomap[0] + ".npz")
            label = d[lithomap[1]][self.locIDs].astype(np.int64)
        elif wfrom == "surface_class":                    # (c) dynamic top layer
            label = self._surfaceSourceClass()
            if label is None and MPIrank == 0 and self.verbose:
                print(
                    "[gw] geochem weatherability_from='surface_class' needs "
                    "provenance on — falling back to scalar weatherability.",
                    flush=True,
                )
        elif wfrom == "source_class":                     # (c) static bedrock
            label = getattr(self, "source_class", None)
            if label is None and MPIrank == 0 and self.verbose:
                print(
                    "[gw] geochem weatherability_from='source_class' needs "
                    "provenance on — falling back to scalar weatherability.",
                    flush=True,
                )
        byclass = getattr(self, "gwGeoWeatherByClass", None)
        out = []
        for k in range(int(self.gwNspecies)):
            wab = self.gwGeoWeather[k]
            tbl = byclass[k] if byclass is not None else None
            if isinstance(wab, (list, tuple)):            # (a) per-vertex map
                d = np.load(wab[0] + ".npz")
                out.append(d[wab[1]][self.locIDs].astype(np.float64))
            elif tbl is not None and label is not None:   # (c) table by label
                tblA = np.asarray(tbl, dtype=np.float64)
                out.append(tblA[np.clip(label, 0, len(tblA) - 1)])
            else:                                          # scalar (unchanged)
                out.append(float(wab))
        return out

    def _weatheringSupply(self):
        r"""
        Solute-supply rate ``Ψ`` feeding in-situ fringe precipitation
        (DESIGN_WATERTABLE_DURICRUST.md §3a). Three modes (``duriWeatherMode``);
        all plug in at the same place, so the formation ODE is mode-agnostic.
        Supply-only (no mass debit) — see §3a Level-A caveat. Rank-local.

        - **proxy** (default): ``Ψ = max(0, rain − evap)^p · arrhenius(T)`` — a
          climate/temperature stand-in, no new inputs.
        - **rate** (Level A): Maher-Chamberlain kinetic×thermodynamic law
          ``W = R·C_eq·(1 − exp(−Dw/(R·L)))·arrhenius·weatherability`` driven by
          the recharge ``R`` the head solve already computes; ``L = Lsoil`` when
          ``soilSPL`` is on (capped by the regolith production supply), else the
          prescribed ``path_length``.
        - **prodsoil**: reuse the temperature-scaled ``soilSPL.prodSoil`` directly,
          scaled by water availability (the chemical∝physical congruency). Falls
          back to the proxy when soil is not tracked.
        """
        rain = self.rainVal
        evap = getattr(self, "evapVal", None)
        net = np.maximum(0.0, rain if evap is None else (rain - evap))
        mode = self.duriWeatherMode

        if mode == "rate":
            R = self.rechargeL.getArray()
            if getattr(self, "cptSoil", False) and getattr(self, "Lsoil", None) is not None:
                L = np.maximum(self.Lsoil.getArray(), 1.0e-3)   # regolith residence length
            else:
                L = max(float(self.duriWeatherL), 1.0e-3)
            # weatherability: scalar OR a lazily-resolved per-vertex map.
            wab = self.duriWeatherability
            if isinstance(wab, (list, tuple)):
                if self._duriWeatherArr is None:
                    d = np.load(wab[0] + ".npz")
                    self._duriWeatherArr = d[wab[1]][self.locIDs].astype(np.float64)
                wab = self._duriWeatherArr
            else:
                wab = float(wab)
            # R→0 ⇒ Dw/(R·L)→∞ ⇒ frac→1 ⇒ W→0 (guard the divide).
            with np.errstate(divide="ignore", invalid="ignore"):
                frac = 1.0 - np.exp(-self.duriWeatherDw / (R * L))
            W = np.where(R > 0.0, R * self.duriWeatherCeq * frac, 0.0) * wab
            W = W * self._arrhenius()
            prod = getattr(self, "prodSoil", None)
            if prod is not None:                                 # regolith supply cap
                W = np.minimum(W, prod * rain)
            return W

        if mode == "prodsoil":
            prod = getattr(self, "prodSoil", None)
            if prod is None:
                if MPIrank == 0 and self.verbose:
                    print(
                        "[gw] duricrust weathering mode 'prodsoil' needs soil "
                        "production — falling back to the climate proxy.",
                        flush=True,
                    )
                self.duriWeatherMode = "proxy"
            else:
                return prod * net

        # proxy (default / fallback)
        return net ** float(self.duriSupplyExp) * self._arrhenius()

    def _updateDuricrust(self):
        r"""
        Evolve the generic capillary-fringe duricrust over ``self.dt``
        (DESIGN_WATERTABLE_DURICRUST.md §3, step 5). Per-node, rank-local — no
        collective; ``localToGlobal`` on ``duriH`` at the end for the halo.

        - **Fringe favourability** ``Φ = exp(−((wt − d0)/w)²)`` — a Gaussian band
          on the water-table depth ``wt = z − h`` centred on the mean fringe depth
          ``d0`` (half-width ``w``); →1 at the fringe, →0 far above/below. With the
          opt-in ``discharge_gate`` (``_dischargeWeight``) ``Φ`` is additionally
          restricted to groundwater discharge zones (absolute-accumulation crust).
        - **Formation** ``dduriH/dt = k_form·Φ·Ψ·(1 − duriH/duriH_max)``
          (self-limiting to ``duriH_max``); ``Ψ`` from ``_weatheringSupply``. When
          soil is tracked (``cptSoil``) the formation supply is additionally
          **regolith-limited** — capped by ``_regolithSupplyRate`` (chemical crust
          growth cannot outpace physical regolith production; §8 soil coupling).
        - **Breakdown**: the per-step surface incision (``z_last − z > 0``) strips
          the crust top (``−k_break·incision``), plus a slow disequilibrium decay
          away from the fringe (``−k_decay·(1−Φ)·duriH``). Clipped to
          ``[0, duriH_max]`` (a crust eroded through resets to 0).

        Writes ``duriH`` (``duriHL``/``duriHG``), the induration ``duriF =
        duriH/duriH_max`` and the armor multiplier ``duriKarmor = 1 −
        armor_max·duriF`` (consumed by the Phase-4 erodibility hook).
        """
        z = self.hLocal.getArray()
        wt = self.wtDepth                                        # z − h ≥ 0
        Hmax = float(self.duriMaxThick)

        Phi = np.exp(-(((wt - self.duriFringeDepth) / self.duriFringeWidth) ** 2))
        # No capillary fringe under standing water: zero the favourability off
        # subaerial land so the crust neither forms nor is held under the sea /
        # a ponded lake (a submerged crust then decays via the (1-Phi) term).
        Phi = np.where(self._subaerialMask(), Phi, 0.0)
        Psi = self._weatheringSupply()

        if getattr(self, "gwGeochemOn", False):
            # Level-B conservative geochemistry provides the crust source (the
            # transported, precipitated solute in `_updateSolute`), so the Level-A
            # proxy/rate supply is switched off here to avoid double-counting.
            supply = np.zeros(self.lpoints, dtype=np.float64)
        else:
            # Absolute-accumulation gate (opt-in): restrict the fringe supply to
            # groundwater discharge zones (lateral convergence of q = −T∇h) rather
            # than forming wherever the water table is shallow.
            if getattr(self, "duriDischargeGate", False):
                h = self.headL.getArray()
                T = self.gwKsat * np.maximum(
                    h - self._gwZbed(z), float(self.gwMinSatThick)
                )
                _, divq = self._soluteAdvecCoeffs(h, T)
                Phi = Phi * self._dischargeWeight(divq, self.rechargeL.getArray())
            supply = self.duriFormRate * Phi * Psi
            if getattr(self, "cptSoil", False):                  # regolith-limited
                supply = np.minimum(supply, self._regolithSupplyRate())

        duriH = self.duriHL.getArray().copy()
        duriH += self.dt * supply * (1.0 - duriH / Hmax)
        incision = np.maximum(0.0, self.gwZlast - z)             # surface lowered
        duriH -= self.duriBreakRate * incision
        duriH -= self.dt * self.duriDecayRate * (1.0 - Phi) * duriH
        duriH = np.clip(duriH, 0.0, Hmax)

        self.duriHL.setArray(duriH)
        self.dm.localToGlobal(self.duriHL, self.duriHG)
        # Refresh the halo from the owners (localToGlobal INSERT took the owned
        # values; globalToLocal broadcasts them back) so the crust `duriH`, the
        # induration `duriF` and the armor `duriKarmor` are halo-consistent and the
        # `duricrust`/`induration`/`Karmor` outputs show no partition seam.
        self.dm.globalToLocal(self.duriHG, self.duriHL)
        duriH = self.duriHL.getArray()
        self.duriF = duriH / Hmax
        self.duriKarmor = 1.0 - self.duriArmorMax * self.duriF
        self.gwZlast = z.copy()
        return

    def _regolithSupplyRate(self):
        r"""
        Rate (m/yr) at which weathering-produced regolith can feed the duricrust,
        the cap on chemical crust formation when soil is tracked (§8): the soil
        production rate ``prodSoil · rain`` (chemical∝physical weathering
        congruency). Returns ``+inf`` (no cap) when ``prodSoil`` is unavailable,
        so a soil-off run is unaffected. Rank-local.
        """
        prod = getattr(self, "prodSoil", None)
        if prod is None:
            return np.inf
        return prod * self.rainVal

    def _surfaceSlope(self):
        r"""
        Per-node steepest-descent slope (m/m, ≥ 0) from the flow-direction
        receivers built by ``flowAccumulation`` — ``(z − z_rcv)/dist`` to the
        primary receiver ``rcvID[:,0]``. Flats / sinks / outlets (``dist=0`` or a
        self-receiver) return 0. A cheap proxy for the optional slope-dependent
        infiltration (goSPL stores no explicit slope). Rank-local.
        """
        hl = self.hLocal.getArray()
        rcv0 = self.rcvID[:, 0]
        d = self.distRcv[:, 0]
        valid = (rcv0 >= 0) & (d > 0.0)
        idx = np.where(valid, rcv0, np.arange(self.lpoints))
        slope = np.where(valid, (hl - hl[idx]) / np.where(d > 0.0, d, 1.0), 0.0)
        return np.maximum(slope, 0.0)

    def _recordInduration(self):
        r"""
        Sync the live per-node induration ``duriF`` with the per-layer
        stratigraphic archive ``stratDuri`` (DESIGN_WATERTABLE_DURICRUST.md §9).
        Two directions:

        - **Exhumation (read-up)** — a previously buried, indurated layer now at
          the surface (its overburden eroded through last step) re-arms the live
          crust: ``duriF = max(duriF, stratDuri[top])`` at the **top non-empty
          layer** (found as in ``_surfaceComposition``), refreshing ``duriKarmor``.
          The stacked-duricrust / relief-inversion behaviour of cratonic laterites.
        - **Formation (write-down)** — the live crust has a real thickness
          ``duriH`` that spans a **depth range**, not just the surface layer, so
          the induration is written into **every layer whose top lies within
          ``duriH`` below the surface**: ``stratDuri[k] = max(stratDuri[k], duriF)``
          for all ``k`` with ``depth_above(k) < duriH`` (``depth_above`` = the
          summed thickness of the layers above ``k``). A thick crust therefore
          indurates several thin layers; each is **preserved** when later buried
          (``deposeStrat`` fresh layers start at 0) and re-arms the surface when
          re-exhumed — so an exhumed crust resists incision over its full
          thickness, not one layer's worth.

        With Level-B geochemistry on, the same write-down also stamps the
        crust's **dominant solute species** into ``stratCrustType`` (and, with
        in-model provenance, its **dominant source region** into
        ``stratCrustSource``) across the crust's depth range, so a stratigraphic
        section preserves *what* each crust layer is and *where* its chemistry
        came from — the categorical companions to the ``stratDuri`` degree.
        These are archive-only (no exhumation read-up).

        No-op (surface-only ``duriF``, no archive) when ``stratDuri`` is
        unallocated (``stratNb == 0``). Composition-only — no geometry change.
        Rank-local (per-node); no collective.
        """
        if getattr(self, "stratDuri", None) is None or self.stratNb == 0:
            return
        top = self.stratStep + 1
        H = self.stratH[:, :top]
        rev = (H > 0)[:, ::-1]
        valid = rev.any(axis=1)                      # columns with any sediment
        if not valid.any():
            return

        # Read-up: the exposed (top non-empty) layer re-arms the live crust.
        top_idx = (H.shape[1] - 1 - np.argmax(rev, axis=1))[valid]
        rows = np.arange(H.shape[0])[valid]
        self.duriF[valid] = np.maximum(self.duriF[valid], self.stratDuri[rows, top_idx])
        # The read-up is rank-local and reads `stratDuri`, so sync `duriF` from the
        # owners (owner wins) before deriving the armor, keeping the induration /
        # Karmor outputs free of a partition seam even when a crust is exhumed.
        self.duriF = self._syncHalo(self.duriF)
        self.duriKarmor = 1.0 - self.duriArmorMax * self.duriF

        # Write-down: record duriF into every layer within `duriH` of the surface.
        # depth_above(k) = Σ_{j>k} H[j] (thickness overlying layer k, so the top
        # non-empty layer has depth_above 0 and is included whenever duriH > 0).
        duriH = self.duriHL.getArray()
        depth_above = np.cumsum(H[:, ::-1], axis=1)[:, ::-1] - H
        within = (H > 0) & (depth_above < duriH[:, None])
        self.stratDuri[:, :top] = np.where(
            within,
            np.maximum(self.stratDuri[:, :top], self.duriF[:, None]),
            self.stratDuri[:, :top],
        )

        # Level-B chemistry archive: stamp the crust's DOMINANT solute species
        # (and, with provenance, its dominant source region) into every layer
        # of the crust's depth range, wherever crust has actually precipitated
        # (gwCrustType >= 0). Categorical, so a later step overwrites rather than
        # accumulates — the layer carries the current dominant while it is within
        # the live crust, then freezes when buried below duriH. Archive-only (no
        # read-up); the live gwCrustType/gwCrustSource are recomputed each step.
        if getattr(self, "stratCrustType", None) is not None:
            paint = within & (self.gwCrustType[:, None] >= 0)
            self.stratCrustType[:, :top] = np.where(
                paint, self.gwCrustType[:, None].astype(np.float64),
                self.stratCrustType[:, :top],
            )
        if getattr(self, "stratCrustSource", None) is not None:
            paint = within & (self.gwCrustSource[:, None] >= 0)
            self.stratCrustSource[:, :top] = np.where(
                paint, self.gwCrustSource[:, None].astype(np.float64),
                self.stratCrustSource[:, :top],
            )
        return

    # ------------------------------------------------------------------ #
    #  Level-B geochemistry — solute transport (G1: operator + solver).   #
    #  Steady advection ∇·(q c) of a lumped conservative tracer along the #
    #  groundwater flux q = −T∇h, reusing the head operator's face        #
    #  conductances. NOT wired into updateGroundwater yet (inert); the    #
    #  dissolution source, fringe precipitation and export land in G2/G3. #
    # ------------------------------------------------------------------ #

    def _makeSoluteKSP(self):
        """
        Cached KSP for the solute-transport solve: a **direct LU factorisation**
        (``preonly`` + ``lu``; MUMPS in parallel), ``gw_solute_`` prefix.

        The upwind advection operator is a non-symmetric M-matrix, but its
        diagonal is only weakly dominant: at a flat, poorly-drained *saturated*
        node the lateral outflow, the seepage sink and the recharge floor are all
        small, so the operator is **badly conditioned**. A Krylov + block-Jacobi
        solve then STALLS there for the strong (high-weatherability) tracers —
        it returns a non-converged iterate with a concentration ~10× too large,
        which silently breaks the per-step ``dissolved = precipitated + exported``
        mass balance (the ocean/marine budget then drifts by an order of
        magnitude). A direct solve is exact, so it is both robust to the
        conditioning and conservative to round-off; the operator is the same size
        / sparsity as the flexure biharmonic that already uses this path, and the
        per-species factorisation is reused across the provenance right-hand
        sides. Env-overridable via the ``gw_solute_`` prefix (e.g. request an
        iterative solver for meshes too large to factorise).
        """
        ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
        ksp.setType("preonly")
        pc = ksp.getPC()
        pc.setType("lu")
        if MPIsize > 1:
            pc.setFactorSolverType("mumps")              # parallel direct solve
        ksp.setOptionsPrefix("gw_solute_")
        ksp.setFromOptions()
        return ksp

    def _soluteAdvecCoeffs(self, h, T):
        r"""
        Upwind finite-volume coefficients for the steady solute advection
        ``∇·(q c)`` by the groundwater flux ``q = −T∇h`` (area-normalised, per
        cell). Reuses ``jacobiancoeff`` — its off-diagonals are the
        area-normalised face conductances ``C_ik/A_i`` (correct for **flat and
        global** meshes alike), so the signed face flux ``i→k`` is
        ``f_ik = (C_ik/A_i)·(h_i − h_k)``. First-order upwinding puts the
        outflow (``f>0``) on the diagonal (carries ``c_i``) and the inflow
        (``f<0``) on the neighbour column (carries ``c_k``). Returns
        ``(adv, divq)``: ``adv`` is the ``(lpoints, 1+maxnb)`` array for
        ``_assembleDiffMatCSR`` (col 0 = diagonal), and ``divq = Σ_k f_ik`` is the
        lateral divergence (used for the seepage sink / export).
        """
        zeroKp = np.zeros(self.lpoints, dtype=np.float64)
        lap = jacobiancoeff(h, T, zeroKp)                 # area-norm neg-Laplacian
        ncol = lap.shape[1] - 1
        cond = -lap[:, 1:]                                # C_ik/A_i ≥ 0 (off-diag)
        flux = np.zeros((self.lpoints, ncol), dtype=np.float64)
        for k in range(ncol):
            flux[:, k] = cond[:, k] * (h - h[self.FVmesh_ngbID[:, k]])
        adv = np.zeros((self.lpoints, 1 + ncol), dtype=np.float64)
        adv[:, 0] = np.maximum(flux, 0.0).sum(axis=1)     # outflow → c_i (diagonal)
        adv[:, 1:] = np.minimum(flux, 0.0)                # inflow  → c_k (neighbour)
        # Vertical exchange closes the balance: the lateral divergence
        # div q = Σ_k f_ik equals recharge (source, >0) minus seepage (sink, <0).
        # At a DISCHARGE node (net lateral inflow, div q < 0) the solute leaves the
        # aquifer to the surface — a diagonal sink `−div q` that makes the operator
        # a well-posed, diagonally-dominant M-matrix (pure lateral advection has no
        # sink there and blows up). Recharge nodes carry their source in the RHS.
        divq = flux.sum(axis=1)
        adv[:, 0] += np.maximum(-divq, 0.0)
        return adv, divq

    def _dischargeWeight(self, divq, R):
        r"""
        Absolute-accumulation gate (opt-in ``discharge_gate``): a per-node weight
        ``G ∈ [0,1]`` that restricts crust formation to groundwater **discharge**
        zones. The lateral groundwater flux ``q = −T∇h`` **converges** into a
        discharge cell (``div q < 0``), delivering dissolved load that emerges at
        the surface and precipitates — the *absolute-accumulation* (lateral) style
        of valley/footslope ferricrete, as opposed to the in-situ / relative
        accumulation that indurates wherever the water table is shallow.

        ``G = (−div q)⁺ / ((−div q)⁺ + R)`` — the fraction of the cell's upward
        discharge that was **imported laterally** (vs supplied by local recharge
        ``R``). Both terms are area-normalised rates (m/yr), so ``G`` is
        dimensionless and self-scaling (no new tuning constant): ``G→1`` at a
        strongly convergent valley floor, ``G→0`` at a recharge (divergent) rise.
        Returns all-ones when the gate is off (backwards-compatible). Rank-local.
        """
        if not getattr(self, "duriDischargeGate", False):
            return np.ones(self.lpoints, dtype=np.float64)
        dis = np.maximum(-divq, 0.0)                       # lateral convergence
        return dis / (dis + np.maximum(R, 0.0) + 1.0e-12)

    def _solveSoluteTransport(self, source, dmask, dval):
        r"""
        Solve one steady tracer transport ``M c = source`` with Dirichlet nodes
        ``dmask`` pinned to ``dval`` (``M`` = upwind advection of §``_soluteAdvecCoeffs``
        at the current head). The seepage set is pinned (``c`` leaves the aquifer
        there — the export sink), which anchors the M-matrix. G1 machinery: the
        physical dissolution source / fringe precipitation / baseflow export are
        added in G2/G3. Collective (KSP); returns the local concentration array.
        """
        z = self.hLocal.getArray()
        h = self.headL.getArray()
        T = self.gwKsat * np.maximum(h - self._gwZbed(z), float(self.gwMinSatThick))
        adv, _ = self._soluteAdvecCoeffs(h, T)
        M = self._assembleDiffMatCSR(adv)
        IntType = petsc4py.PETSc.IntType
        owned_d = np.where(dmask & (self.inIDs == 1))[0].astype(IntType)
        M.zeroRowsLocal(owned_d, diag=1.0)                # c = dval on Dirichlet rows

        rhs = np.asarray(source, dtype=np.float64).copy()
        rhs[dmask] = dval[dmask]
        if self._ksp_solute is None:
            self._ksp_solute = self._makeSoluteKSP()
        ksp = self._ksp_solute
        self.soluteL.setArray(rhs)
        self.dm.localToGlobal(self.soluteL, self.tmp)     # rhs (global)
        self.dm.localToGlobal(self.soluteL, self.soluteG)  # nonzero guess = rhs
        ksp.setOperators(M, M)
        ksp.solve(self.tmp, self.soluteG)
        M.destroy()
        self.dm.globalToLocal(self.soluteG, self.soluteL)
        return self.soluteL.getArray().copy()

    def _soluteSolveRHS(self, M, rhs, guess):
        r"""
        Solve ``M c = rhs`` for the (pre-assembled) transport operator ``M`` with
        a warm-start ``guess``, returning the local concentration (clipped ≥ 0).
        A thin wrapper so ``_updateSolute`` can reuse one assembled operator for
        several right-hand sides (the per-source-class provenance solves, G5).
        Collective (KSP).
        """
        if self._ksp_solute is None:
            self._ksp_solute = self._makeSoluteKSP()
        ksp = self._ksp_solute
        self.soluteL.setArray(rhs)
        self.dm.localToGlobal(self.soluteL, self.tmp)          # rhs (global)
        self.soluteL.setArray(guess)
        self.dm.localToGlobal(self.soluteL, self.soluteG)      # nonzero guess
        ksp.setOperators(M, M)
        ksp.solve(self.tmp, self.soluteG)
        self.dm.globalToLocal(self.soluteG, self.soluteL)
        return np.maximum(self.soluteL.getArray().copy(), 0.0)

    def _updateSolute(self):
        r"""
        Level-B per-step solute update (G2), per tracer: **dissolve → transport →
        precipitate → export**, conservatively accounted.

        1. **Dissolution** — a chemical-weathering source (the Level-A
           ``_weatheringSupply`` scaled per tracer by ``weatherability``), on
           subaerial land only, **debiting the conserved source pool**.
        2. **Transport** — the steady ``M c = D`` solve of ``_soluteAdvecCoeffs``
           (upwind advection by ``q = −T∇h`` + the vertical seepage sink), plus a
           **precipitation sink** at the capillary fringe added to the diagonal
           (``k_p·Φ`` where the tracer is super-saturated — a one-step-lagged
           ``c > c_sat`` gate).
        3. **Precipitation** — the sink mass ``k_p·Φ·c`` feeds the duricrust
           ``duriH`` (Level-B thus **replaces** the Level-A proxy supply as the
           crust source; the induration/armoring then follow unchanged).
        4. **Export** — by domain mass balance, the solute that is neither
           precipitated nor left in solution has discharged to the surface network
           (→ ocean; the flux is formalised in G3).

        Budget per tracer (owned nodes): ``dissolved = precipitated + exported +
        Δ(in solution)`` — accumulated in ``gwDissolved``/``gwPrecip``/``gwOceanFlux``
        for the conservation guard. Rank-local accounting (KSP solve collective).
        """
        if not getattr(self, "gwGeochemOn", False):
            return
        z = self.hLocal.getArray()
        h = self.headL.getArray()
        T = self.gwKsat * np.maximum(h - self._gwZbed(z), float(self.gwMinSatThick))
        A = self.larea
        dt = self.dt
        owned = self.inIDs == 1

        adv, divq = self._soluteAdvecCoeffs(h, T)  # advection + seepage sink (G1)
        seep_sink = np.maximum(-divq, 0.0)         # lateral-convergence discharge
        # Vertical discharge to the surface = recharge − lateral divergence, from the
        # quasi-steady balance ∇·q = R − seepage ⇒ seepage = R − ∇·q. Apply it at
        # EVERY node, not just the saturated ones: wherever the groundwater cannot
        # transmit the dissolved load it must leave with the recharge/runoff to the
        # surface network, and there are TWO such near-singular populations —
        #   (i)  a flat, fully SATURATED seepage node (∇·q ≈ 0, no lateral outflow);
        #   (ii) a DRY node (water table at the aquifer base, so `T` is pinned at
        #        `min_sat_thickness` and the face conductances — hence the transport
        #        diagonal — collapse), common over arid interiors.
        # Both are left with a ~zero transport diagonal by the `−div q` sink alone,
        # so the steady concentration (and the ocean/marine flux) spikes. Adding the
        # recharge term floors the diagonal at ~R. At a well-drained node ∇·q ≈ R, so
        # the extra term is ~0 and the solute advects laterally exactly as before.
        R = self.rechargeL.getArray()
        extra = np.maximum(np.maximum(R - divq, 0.0) - seep_sink, 0.0)  # only ADD discharge
        adv[:, 0] += extra                          # keep the operator diagonal in sync
        seep_sink = seep_sink + extra
        self.gwSoluteFlux[:] = 0.0                 # per-node baseflow export (G3)
        self.gwSoluteFluxSp[:] = 0.0               # per-species export (output + routing)
        prov = getattr(self, "provOn", False) and getattr(self, "source_class", None) is not None
        W = self._weatheringSupply()               # base weathering rate (Level A)
        # ext 1: per-species weatherability. Resolved once and cached, EXCEPT the
        # dynamic surface-lithology form, which is re-resolved every step so the
        # weatherability follows the top stratigraphic layer (exhumation/burial).
        if self._gwGeoWeatherArr is None or getattr(self, "gwWeatherFrom", None) == "surface_class":
            self._gwGeoWeatherArr = self._resolveGeoWeather()
        # Subaerial land only — no dissolution AND no fringe precipitation under
        # the sea / a ponded lake (there is no capillary fringe under standing
        # water, so a duricrust cannot form there).
        subaerial = self._subaerialMask()
        # Precipitation favourability = the capillary fringe (needs the
        # duricrust), zeroed off subaerial land.
        if getattr(self, "duriOn", False):
            Phi = np.exp(
                -(((self.wtDepth - self.duriFringeDepth) / self.duriFringeWidth) ** 2)
            )
            Phi = np.where(subaerial, Phi, 0.0)
            # Absolute-accumulation gate: precipitate only where the groundwater
            # discharges (converges), not everywhere the fringe is shallow.
            Phi = Phi * self._dischargeWeight(divq, R)
        else:
            Phi = np.zeros(self.lpoints, dtype=np.float64)

        duriH = self.duriHL.getArray().copy()
        Hmax = float(self.duriMaxThick)
        for k in range(int(self.gwNspecies)):
            # 1. Dissolution — debit the source pool.
            # weatherability: scalar OR a per-vertex array (ext 1) — broadcasts.
            Drate = np.where(subaerial, self._gwGeoWeatherArr[k] * W, 0.0)
            want = Drate * A * dt
            diss = np.minimum(want, self.gwSourcePool[:, k])
            self.gwSourcePool[:, k] -= diss
            Deff = diss / (A * dt)
            # One-time warning when the finite source pool starts to limit
            # dissolution (the tracer will go inert as it empties) — a silent
            # exhaustion otherwise looks like the geochemistry has stopped. Raise
            # `source_pool` for that species to keep weathering rate-limited. The
            # `not _solutePoolWarned` gate is uniform across ranks (the flag is set
            # from a reduced count), so the Allreduce below is reached collectively.
            if not self._solutePoolWarned:
                limited = MPI.COMM_WORLD.allreduce(
                    int(((want > diss + 1.0e-30) & (self.inIDs == 1)).sum()),
                    op=MPI.SUM,
                )
                if limited > 0:
                    self._solutePoolWarned = True
                    if MPIrank == 0:
                        print(
                            "[gw] geochem: source pool exhausting for '%s' at %d "
                            "node(s) — dissolution now supply-limited (raise "
                            "`source_pool` to keep it rate-limited)."
                            % (self.gwGeoName[k], limited),
                            flush=True,
                        )

            # 2. Precipitation sink at the fringe — a LINEAR removal `k_p·Φ`
            # (proportional to the local concentration), **self-limiting** as the
            # crust fills: `·(1 − duriH/duriH_max)` (mirrors the Level-A ODE), so a
            # crust at `max_thickness` REJECTS further solute — precipitation → 0
            # there and the solute is exported to the rivers instead of piling into
            # an already-full crust. Uses the running `duriH` (updated per species),
            # so it also shares the remaining room across tracers. Kept linear so
            # the operator is constant under constant forcing → the solute reaches
            # a true per-step steady state (a hard `c > c_sat` gate would
            # oscillate; `c_sat` is a documented nonlinear refinement, DESIGN §3).
            p = self.gwGeoPrecip[k] * Phi * np.maximum(0.0, 1.0 - duriH / Hmax)

            # 3. Transport: (advection + seepage sink + precip sink) c = Deff.
            coeffs = adv.copy()
            coeffs[:, 0] += p
            M = self._assembleDiffMatCSR(coeffs)
            if prov:
                # G5 solute-source provenance: by linearity of the (fixed)
                # operator, transport the solute dissolved in each source-rock
                # class separately (same M, class-restricted RHS), sum to the
                # total, and attribute the precipitated crust to each source.
                c = np.zeros(self.lpoints, dtype=np.float64)
                for r in range(int(self.provNb)):
                    Deff_r = np.where(self.source_class == r, Deff, 0.0)
                    c_r = self._soluteSolveRHS(M, Deff_r, self.gwSolute[:, k])
                    c += c_r
                    self.gwCrustProv[:, r] += p * c_r * dt * self.gwGeoVsolid[k]
            else:
                c = self._soluteSolveRHS(M, Deff, self.gwSolute[:, k])
            M.destroy()
            c = np.maximum(c, 0.0)                               # guard tiny negatives
            self.gwSolute[:, k] = c

            # 4. Sinks: precipitation feeds the crust; seepage exports to the
            # surface (baseflow). The steady operator is exactly conservative
            # (upwind internal faces cancel), so per step
            # dissolved = precipitated + exported to the solver tolerance — no
            # storage term (a steady solve maintains, not accumulates, the
            # standing concentration).
            precip_mass = p * c * A * dt                     # → crust
            export_mass = seep_sink * c * A * dt             # → surface / ocean
            duriH += precip_mass / A * self.gwGeoVsolid[k]
            self.gwDissolved[k] += float(diss[owned].sum())
            self.gwPrecip[k] += float(precip_mass[owned].sum())
            self.gwOceanFlux[k] += float(export_mass[owned].sum())
            # Per-node baseflow export rate (m³/yr): total (summed over tracers)
            # and the per-species split (output + the river-routing source).
            self.gwSoluteFlux += seep_sink * c * A
            self.gwSoluteFluxSp[:, k] = seep_sink * c * A
            # Per-node crust contributed by this tracer (G4 typing).
            self.gwCrustBySpecies[:, k] += precip_mass / A * self.gwGeoVsolid[k]

        # The seepage-export fields are `seep_sink·c·A`; `seep_sink` (from the
        # local FV advection stencil) is only correct on OWNED nodes, so sync the
        # halo from the owner (local -> global -> local) — otherwise the
        # `soluteflux` / `soluteflux_<name>` outputs show a partition seam (like
        # baseflow). The river routing is unaffected (it uses localToGlobal INSERT,
        # owner wins). `gwSolute` / `riverSolute` come from global solves and are
        # already halo-consistent.
        self.gwSoluteFlux = self._syncHalo(self.gwSoluteFlux)
        for k in range(int(self.gwNspecies)):
            self.gwSoluteFluxSp[:, k] = self._syncHalo(self.gwSoluteFluxSp[:, k])

        duriH = np.clip(duriH, 0.0, Hmax)
        self.duriHL.setArray(duriH)
        self.dm.localToGlobal(self.duriHL, self.duriHG)
        # Refresh the halo from the owners so the crust / induration / armor are
        # halo-consistent (no partition seam in the outputs); see _updateDuricrust.
        self.dm.globalToLocal(self.duriHG, self.duriHL)
        if getattr(self, "duriOn", False):
            duriH = self.duriHL.getArray()
            self.duriF = duriH / Hmax
            self.duriKarmor = 1.0 - self.duriArmorMax * self.duriF

        # G4 typing: the dominant crust-forming tracer per node (−1 = no crust).
        tot = self.gwCrustBySpecies.sum(axis=1)
        self.gwCrustType = np.where(
            tot > 0.0, self.gwCrustBySpecies.argmax(axis=1), -1
        ).astype(np.int32)
        # G5 provenance: the dominant SOURCE-ROCK class of the crust per node.
        if prov:
            ptot = self.gwCrustProv.sum(axis=1)
            self.gwCrustSource = np.where(
                ptot > 0.0, self.gwCrustProv.argmax(axis=1), -1
            ).astype(np.int32)

        if self.verbose:
            tot = MPI.COMM_WORLD.allreduce(
                float(self.gwSoluteFlux[owned].sum()), op=MPI.SUM
            )
            if MPIrank == 0:
                per = ", ".join(
                    "%s=%0.3g" % (self.gwGeoName[k], self.gwOceanFlux[k])
                    for k in range(int(self.gwNspecies))
                )
                print(
                    "[gw] dissolved solute flux to surface: %0.4g m3/yr "
                    "(cumulative per tracer: %s)" % (tot, per),
                    flush=True,
                )
        return

    def _routeRiverSolute(self):
        r"""
        Route the groundwater-exported (seepage/baseflow) solute **down the
        surface drainage network** to the shoreline — the river dissolved load
        (Level-B extension 2, DESIGN_GEOCHEM_EXTENSIONS.md §2), **per species**.

        For each tracer an implicit accumulation solve reuses the **cached
        flow-accumulation matrix** ``fMati`` (the ``_getSedFlux`` pattern) with
        that species' per-node seepage export ``gwSoluteFluxSp[:, k]`` (already an
        extensive m³/yr flux — no area weighting) as the source:

        - **conservative** (``river_decay = 0``): ``(I − Wᵀ) L = s`` — a passive
          tracer; solute reaching a coast/outlet leaves the continent, solute
          into a closed basin is trapped (evaporite), both from the filled-topo
          matrix with no special handling.
        - **in-transit loss** (``river_decay = κ > 0``): a first-order removal
          along the network adds a diagonal sink, ``(I − Wᵀ + κ I) L = s`` (built
          as ``fMati.shift(κ)``); the lost mass ``κ·Σ L`` is a per-species
          diagnostic (in-channel precipitation / uptake).

        Per-species fields ``riverSoluteSp`` / ``riverSoluteToOceanSp`` /
        ``riverSoluteLost`` are summed into the totals ``riverSolute`` /
        ``riverSoluteToOcean``. With ``marine_coupling`` on, the delivered coastal
        flux accumulates into a per-species marine reservoir ``marineSolute``
        (m³, integrated) and the per-node coastal input ``marineSoluteInput``.

        Called at the end of ``updateGroundwater`` — ``flowAccumulation`` (which
        caches ``fMati``) runs earlier in the step, so the matrix is fresh and no
        reordering is needed. No-op unless ``river_load`` is on; the flow matrix
        must exist (guarded for the first call). Collective (KSP per tracer).
        """
        if not getattr(self, "gwRiverLoad", False):
            return
        if getattr(self, "fMati", None) is None:
            return                                    # no flow matrix yet
        owned = self.inIDs == 1
        exitm = np.zeros(self.lpoints, dtype=bool)    # coast / outlet exits
        exitm[self.seaID] = True
        outl = getattr(self, "outletIDs", None)
        if outl is not None:
            exitm[outl] = True

        self.riverSolute[:] = 0.0
        for k in range(int(self.gwNspecies)):
            # RHS = this species' per-node seepage export (m³/yr, extensive).
            self.soluteL.setArray(self.gwSoluteFluxSp[:, k])
            self.dm.localToGlobal(self.soluteL, self.tmp)
            kappa = float(self.gwGeoRiverDecay[k])
            if kappa > 0.0:                           # in-transit first-order loss
                M = self.fMati.copy()
                M.shift(kappa)                        # (I − Wᵀ) + κ I
                self._solve_KSP(False, M, self.tmp, self.riverSoluteG)
                M.destroy()
            else:                                     # conservative accumulation
                self._solve_KSP(False, self.fMati, self.tmp, self.riverSoluteG)
            self.dm.globalToLocal(self.riverSoluteG, self.riverSoluteL)
            Lk = np.maximum(self.riverSoluteL.getArray().copy(), 0.0)
            self.riverSoluteSp[:, k] = Lk
            self.riverSolute += Lk
            self.riverSoluteToOceanSp[k] = MPI.COMM_WORLD.allreduce(
                float(Lk[owned & exitm].sum()), op=MPI.SUM
            )
            # In-transit loss (mass balance: Σs = Σ_terminal L + κ·Σ L).
            self.riverSoluteLost[k] = (
                MPI.COMM_WORLD.allreduce(float(kappa * Lk[owned].sum()), op=MPI.SUM)
                if kappa > 0.0 else 0.0
            )
        self.riverSoluteToOcean = float(self.riverSoluteToOceanSp.sum())

        # Marine coupling: the delivered coastal flux feeds a per-species marine
        # reservoir (integrated over time) + a per-node coastal-input field.
        if getattr(self, "gwMarineCoupling", False):
            self.marineSoluteInput[:] = 0.0
            self.marineSoluteInput[exitm] = self.riverSolute[exitm]
            self.marineSolute += self.riverSoluteToOceanSp * self.dt

        if MPIrank == 0 and self.verbose:
            per = ", ".join(
                "%s=%0.3g" % (self.gwGeoName[k], self.riverSoluteToOceanSp[k])
                for k in range(int(self.gwNspecies))
            )
            print(
                "[gw] river dissolved load to ocean: %0.4g m3/yr (per tracer: %s)"
                % (self.riverSoluteToOcean, per),
                flush=True,
            )
        return
