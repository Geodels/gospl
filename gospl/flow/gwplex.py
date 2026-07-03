import os
import petsc4py
import numpy as np

from time import process_time
from mpi4py import MPI

from gospl.tools.constants import ICE_COVER_MIN

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import jacobiancoeff

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()


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
                # Per-species parameters (parser lists -> numpy, len n_species).
                self.gwGeoWeather = np.asarray(self.gwGeoWeather, dtype=np.float64)
                self.gwGeoCsat = np.asarray(self.gwGeoCsat, dtype=np.float64)
                self.gwGeoPrecip = np.asarray(self.gwGeoPrecip, dtype=np.float64)
                self.gwGeoVsolid = np.asarray(self.gwGeoVsolid, dtype=np.float64)
                # Groundwater solute concentration and the dissolvable source pool
                # (rank-local, per tracer). The pool is seeded as an ample per-area
                # reservoir (`1e6 · cell area`), so dissolution is rate-limited
                # (not exhausted) over normal runs; a `source_pool` YAML refinement
                # can later tie it to the actual weatherable rock mass.
                self.gwSolute = np.zeros((self.lpoints, nsp), dtype=np.float64)
                self.gwSourcePool = 1.0e6 * self.larea[:, None] * np.ones(
                    (self.lpoints, nsp), dtype=np.float64
                )
                # Cumulative mass budget per tracer (m³-equiv, owned nodes) for the
                # conservation guard / diagnostics: dissolved = precipitated +
                # exported (ocean) + currently in solution.
                self.gwDissolved = np.zeros(nsp, dtype=np.float64)
                self.gwPrecip = np.zeros(nsp, dtype=np.float64)
                self.gwOceanFlux = np.zeros(nsp, dtype=np.float64)
                # Per-node baseflow-carried solute export (m³/yr, summed over
                # tracers) — the spatial output field (G3).
                self.gwSoluteFlux = np.zeros(self.lpoints, dtype=np.float64)
                # Per-node cumulative crust precipitated by each tracer, and the
                # dominant crust-forming tracer (G4 typing: -1 = no crust).
                self.gwCrustBySpecies = np.zeros((self.lpoints, nsp), dtype=np.float64)
                self.gwCrustType = np.full(self.lpoints, -1, dtype=np.int32)
                # Scratch Vec pair for the per-species transport solve (G1+).
                self.soluteL = self.hLocal.duplicate()
                self.soluteG = self.hGlobal.duplicate()
                self.soluteL.set(0.0)
                self.soluteG.set(0.0)

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

        No-op when ``gwOn`` is off. The head solve is collective (KSP); the
        recharge and duricrust steps are purely rank-local (per-node).
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
        self.wtDepth = z - hloc

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
        seep_owned = seep & owned
        wsum = MPI.COMM_WORLD.allreduce(
            float(np.where(seep_owned, A, 0.0).sum()), op=MPI.SUM
        )
        if wsum > 0.0:
            bf[seep_owned] = Qtot * A[seep_owned] / wsum
        self.baseflowL.setArray(bf)
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
          ``d0`` (half-width ``w``); →1 at the fringe, →0 far above/below.
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
        Psi = self._weatheringSupply()

        if getattr(self, "gwGeochemOn", False):
            # Level-B conservative geochemistry provides the crust source (the
            # transported, precipitated solute in `_updateSolute`), so the Level-A
            # proxy/rate supply is switched off here to avoid double-counting.
            supply = np.zeros(self.lpoints, dtype=np.float64)
        else:
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
        Cached KSP for the solute-transport solve: ``fgmres`` + block-Jacobi
        (``gw_solute_`` prefix, pivot shift). The upwind advection operator is
        non-symmetric but diagonally dominant (an M-matrix, anchored by the
        seepage Dirichlet sink), so — unlike the stiff elliptic *head* operator
        — a Krylov + ILU solve converges quickly (same class as the orographic
        advection solver). Env-overridable via the ``gw_solute_`` prefix.
        """
        ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
        ksp.setType("fgmres")
        ksp.getPC().setType("bjacobi")
        ksp.setTolerances(rtol=1.0e-10, max_it=500)
        ksp.setInitialGuessNonzero(True)
        ksp.setOptionsPrefix("gw_solute_")
        petsc4py.PETSc.Options()["gw_solute_sub_pc_factor_shift_type"] = "nonzero"
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
        seep_sink = np.maximum(-divq, 0.0)         # discharge-to-surface coefficient
        self.gwSoluteFlux[:] = 0.0                 # per-node baseflow export (G3)
        W = self._weatheringSupply()               # base weathering rate (Level A)
        # Subaerial land only (no dissolution under sea / ponded lake).
        sub = np.zeros(self.lpoints, dtype=bool)
        sub[self.seaID] = True
        pitIDs = getattr(self, "pitIDs", None)
        lFill = getattr(self, "lFill", None)
        if pitIDs is not None and lFill is not None:
            sub |= (pitIDs > -1) & (lFill > z)
        subaerial = ~sub
        # Precipitation favourability = the capillary fringe (needs the duricrust).
        if getattr(self, "duriOn", False):
            Phi = np.exp(
                -(((self.wtDepth - self.duriFringeDepth) / self.duriFringeWidth) ** 2)
            )
        else:
            Phi = np.zeros(self.lpoints, dtype=np.float64)

        duriH = self.duriHL.getArray().copy()
        Hmax = float(self.duriMaxThick)
        for k in range(int(self.gwNspecies)):
            # 1. Dissolution — debit the source pool.
            Drate = np.where(subaerial, self.gwGeoWeather[k] * W, 0.0)
            diss = np.minimum(Drate * A * dt, self.gwSourcePool[:, k])
            self.gwSourcePool[:, k] -= diss
            Deff = diss / (A * dt)

            # 2. Precipitation sink at the fringe — a LINEAR removal `k_p·Φ`
            # (proportional to the local concentration). Kept linear so the
            # operator is constant under constant forcing → the solute reaches a
            # true per-step steady state (a hard `c > c_sat` on/off gate instead
            # oscillates). The `c_sat` saturation threshold is a documented
            # refinement needing a nonlinear/Picard treatment (see DESIGN §3).
            p = self.gwGeoPrecip[k] * Phi

            # 3. Transport: (advection + seepage sink + precip sink) c = Deff.
            coeffs = adv.copy()
            coeffs[:, 0] += p
            M = self._assembleDiffMatCSR(coeffs)
            if self._ksp_solute is None:
                self._ksp_solute = self._makeSoluteKSP()
            ksp = self._ksp_solute
            self.soluteL.setArray(Deff)
            self.dm.localToGlobal(self.soluteL, self.tmp)
            self.soluteL.setArray(self.gwSolute[:, k])          # warm start
            self.dm.localToGlobal(self.soluteL, self.soluteG)
            ksp.setOperators(M, M)
            ksp.solve(self.tmp, self.soluteG)
            M.destroy()
            self.dm.globalToLocal(self.soluteG, self.soluteL)
            c = np.maximum(self.soluteL.getArray().copy(), 0.0)  # guard tiny negatives
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
            # Per-node baseflow export rate (m³/yr), summed over tracers, for output.
            self.gwSoluteFlux += seep_sink * c * A
            # Per-node crust contributed by this tracer (G4 typing).
            self.gwCrustBySpecies[:, k] += precip_mass / A * self.gwGeoVsolid[k]

        duriH = np.clip(duriH, 0.0, Hmax)
        self.duriHL.setArray(duriH)
        self.dm.localToGlobal(self.duriHL, self.duriHG)
        if getattr(self, "duriOn", False):
            self.duriF = duriH / Hmax
            self.duriKarmor = 1.0 - self.duriArmorMax * self.duriF

        # G4 typing: the dominant crust-forming tracer per node (−1 = no crust).
        tot = self.gwCrustBySpecies.sum(axis=1)
        self.gwCrustType = np.where(
            tot > 0.0, self.gwCrustBySpecies.argmax(axis=1), -1
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
