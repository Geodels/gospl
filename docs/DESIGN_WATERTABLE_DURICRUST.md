# DESIGN: water table (groundwater) + generic duricrust

Opt-in near-surface hydrology and chemical armoring for goSPL. A parallel implicit
Dupuit–Boussinesq **water-table head** solve on the DMPlex drives a **generic
duricrust** induration state that forms in the fluctuating **capillary fringe**
around the water table and **armors erodibility** (drives relief inversion). Ships
disabled by default; ice-model / dual-lithology convention — **byte-identical to
the current code when off**.

Companion to `DESIGN_DUAL_LITHOLOGY.md`, `DESIGN_ICE_SHEET.md`,
`DESIGN_PROVENANCE.md`. Honors the invariants in `AGENTS.md` (MPI contract,
KSP/SNES lifecycle, scratch-vector contract, `destroy_DMPlex` registration).

> **Status: IMPLEMENTED (branch `feat/watertable-duricrust`).** All phases −1…7
> of the plan in §13 are done and merged behind the opt-in `groundwater:` block:
> the implicit water-table solve (fgmres + hypre AMG, analytic-Dupuit validated),
> recharge, seepage/baseflow, the duricrust ODE + K-armoring hook, the `stratDuri`
> stratigraphic record, soil coupling (`from_soil`, regolith limiter), restart of
> `head`/`duriH`, and the full user/tech/API docs. Guarded by the
> `test_groundwater_*` / `test_duricrust_*` suite (serial + np=2). Since then the
> **baseflow re-injection** into the surface flow (§3 step 7), the **`from_soil`
> basin base** from the stratigraphy (§15), the opt-in **lake ↔ aquifer volume
> coupling** (§15), and the smaller **recharge refinements** — subglacial-meltwater
> recharge, and lithology-/slope-modulated `f_infil` (§3) — have also landed (all
> opt-in, default off). **Remaining deferred increments:** the multi-layer
> formation **depth range** (§9) and the geochemical **Level-B** solute transport
> (§3a/§15) — the latter also being the prerequisite for duricrust **solute-source
> provenance** (§11). The sections below are the original design narrative,
> annotated with "as built" notes where the implementation refined a choice.

---

## 1. Scope & decisions (locked)

| # | Decision | Choice |
|---|---|---|
| 1 | Water-table model | **Full Dupuit–Boussinesq PDE**, solved **implicitly (backward-Euler) once per step** as an elliptic/quasi-steady head field. NOT transient/explicit — no CFL. |
| 2 | Recharge | **`R = f·max(0, rain − evap)`** — reuses the existing per-node `rainVal` / `evapVal` forcing. `f` = infiltration fraction (YAML). |
| 3 | Duricrust representation | **Per-node induration state** (`duriH` thickness + `duriF` degree ∈ [0,1]) that **armors K** through the single `_surfaceLithoK` hook (§5). |
| 4 | Duricrust chemistry | **Generic** — one chemistry-agnostic mechanism with tunable parameters; ferricrete/laterite/calcrete/silcrete emulated by parameter choice, not code branches. |
| 5 | Formation rule | **Water-table fringe** — induration accumulates where the surface sits within a depth band of the water table (the fluctuating capillary fringe), scaled by a temperature/water weathering-supply proxy. |
| 6 | Opt-in | `groundwater:` YAML block. Absent ⇒ every path gated out ⇒ bit-identical to current `dev`. |
| 7 | Target regime | Δt ≈ 10²–10³ yr, Δx ≈ 500 m – km, runs of 10⁶–10⁷ yr, annual-mean rain + temperature. |
| 8 | Soil dependence | Targets the **Option-2.5 soil model** (`DESIGN_SOIL_REGOLITH.md`: regolith-only `Lsoil`, deposits in stratigraphy, subaerial lake/sea gate). **Ships both**: runs **soil-independent** on any config, AND couples to `soilSPL` when tracked (regolith limiter + optional `aquifer_base = from_soil`). See §8. |
| 9 | Stratigraphic record | Armoring is **recorded per layer** (`stratDuri` induration) when `stratNb>0`, so buried crusts are preserved, advect/compact with the pile, and **re-armor on exhumation** (multi-cycle relief inversion). See §9. |
| 10 | Compatibility | **Composes with dual-lithology and provenance**, individually or all three at once: multiplicative erodibility/diffusivity hooks with dual, passive-label reflection with provenance; `stratDuri` is intensive so the conservation guards stay green. See §10. |

---

## 2. Background — why an implicit elliptic solve (and why 500 m is fine)

The unconfined water-table equation is

```
S ∂h/∂t = ∇·(T(h) ∇h) + R − Q_seep       T(h) = K_h·(h − z_bed)   (transmissivity)
```

with `S` specific yield, `K_h` hydraulic conductivity, `R` recharge, `z_bed` the
aquifer base, and a seepage sink `Q_seep` where `h` would exceed the surface `z`.

**Timescale argument.** The water-table diffusivity is `D = T/S = K_h·b/S`
(`b` = saturated thickness). The re-equilibration time over a flow path of length
`L` is `τ_gw ≈ L²/D`. For effective regolith/aquifer values (`S≈0.1`, `b≈20 m`):

| Flow system | L | `K_h≈0.1 m/d` | `K_h≈1 m/d` | `K_h≈10 m/d` |
|---|---|---|---|---|
| local / hillslope | 1 km | ~140 yr | ~14 yr | ~1 yr |
| intermediate | 10 km | ~1.4×10⁴ yr | ~1400 yr | ~140 yr |
| regional | 100 km | ~1.4×10⁶ yr | ~1.4×10⁵ yr | ~1.4×10⁴ yr |

At Δt ≈ 10²–10³ yr, local and most intermediate systems have `τ_gw ≪ Δt` — the
table is in **quasi-equilibrium** with recharge and topography every step. Regional
systems lag, but topography itself only evolves over 10⁴–10⁶ yr, so the table is
quasi-static relative to the landscape too.

**One formulation covers both regimes.** Backward-Euler in time:

```
(S/Δt)·h − ∇·(T ∇h) = (S/Δt)·h_old + R
```

- Δt ≫ τ_gw ⇒ the `S/Δt` term vanishes ⇒ steady elliptic `∇·(T∇h) = −R` (equilibrium limit).
- Δt ≲ τ_gw ⇒ a proper transient lag.
- **Unconditionally stable** — no CFL, no substepping. Discretized as `(a·I + L)h = rhs`
  with `a = S/Δt`, which is exactly the cached FV operator form goSPL already assembles
  (`hillslope._buildDiffMat` / `_assembleDiffMatCSR`, `_makeDiffusionKSP`).

**500 m resolution.** Stability is independent of Δx (implicit). Finer Δx resolves
*shorter* flow paths (smaller `τ_gw`), pushing further into the quasi-equilibrium
regime — never toward stiffness. 500 m sharpens the valley-shallow / ridge-deep
water-table structure that governs duricrust, at the usual per-node cost (one extra
diffusion solve of the mesh's size per step; well within goSPL's demonstrated
5.9–23.7 M-node HPC envelope). This is the **opposite** of the SIA ice case removed
in `DESIGN_ICE_SHEET.md`: that failed on the `H≥0` obstacle margin, not the diffusion;
the groundwater seepage boundary is far milder (§4).

**What we deliberately do NOT model** (validity, not feasibility): sub-500 m local
flow cells (sub-grid), confined / multi-layer hydrostratigraphy (single unconfined
layer), and calibrated absolute heads. Parameters are **effective/lumped** (one `K_h`
per lithology, effective `b`, `S`); the robust, usable output is the **water-table
depth-below-surface pattern**, which is what drives induration.

---

## 3. The algorithm (per goSPL step)

New mixin `_GWMesh` in `gospl/flow/gwplex.py` (mirrors `flow/iceplex.py`). Runs in
`runProcesses` **after `flowAccumulation` (`model.py:389`) and before erosion
(`model.py:392`)** — flow gives the drainage graph (seepage nodes) and discharge;
groundwater sets the armoring state that erosion then reads.

`updateGroundwater()`:

1. **Recharge.** `R = f_infil · max(0, rainVal − evapVal)` (m/yr), per node, from the
   existing forcing arrays. `R = 0` where the surface is **not subaerial land**: `seaID` and
   ponded continental lakes (head pinned; see step 3) **and ice-covered land** (`iceHL >
   ICE_COVER_MIN` — rain falls as snow/ice and does not infiltrate the ground; parallels the
   soil ice-freeze gate. **As built:** an opt-in `subglacial_recharge` fraction lets the ice
   model's `iceMeltRiverL` infiltrate under ice — the one recharge path allowed there). `f_infil`
   may be a **scalar or a per-vertex map** `[file, key]` (loaded in `_GWMesh`) and is optionally
   modulated (opt-in, both default off) by **surface lithology** (`fine_infil_factor` on the exposed
   coarse fraction, dual lithology) and by **slope** (`f/(1+slope/infil_slope_ref)`, the
   steepest-descent gradient) — DESIGN §3 recharge refinements.
2. **Seepage set.** Nodes where the table is pinned to the surface: rivers/lakes
   (drainage-connected, from the flow graph) + coast/sea (`seaID`) + open boundary
   outlets (`outletIDs`). Dirichlet `h = z` there (partition-invariant — derived from
   the already-invariant drainage arrays, `AGENTS.md` mechanism-#2 work).
3. **Head solve** (`_solveHead`): build `(a·I + L(T))h = (a·I)·h_old + R·Δt` with
   `a = S/Δt`, `L(T)` the FV neg-Laplacian weighted by transmissivity `T = K_h·(h−z_bed)`.
   - **Unconfined nonlinearity** in `T(h)`: 2–3 **Picard** iterations lagging `T` on the
     previous iterate (exactly the `_diffuseImplicitPicard` pattern). Break on an
     `Allreduce`d head-change norm.
   - **Seepage free boundary** (`h ≤ z`): after each Picard solve, clip `h = min(h, z)`;
     any newly-clipped node is added to the Dirichlet set and the solve repeated. A
     mild, well-posed fixed point (2–4 outer passes typical). The "re-solve?" decision
     is **reduced across ranks** (`allreduce(any_new_seep, MPI.LOR)`) before it gates
     the collective re-solve (`AGENTS.md` #1 deadlock rule).
   - Cached `gw_`-prefixed KSP (`_makeGWKSP`: **fgmres + hypre BoomerAMG**); operator
     rebuilt per step (T varies), KSP object reused. AMG is **required**: block-Jacobi/ILU
     stalls (`DIVERGED_ITS`) on this stiff 2-D elliptic operator and returns a garbage
     iterate that drifts the water table up to the surface (same reason the flexure
     biharmonic wants a strong solver).
4. **Water-table depth** `wtL = z − h` (≥ 0 by the clip), the field the duricrust reads.
5. **Duricrust update** (`_updateDuricrust`, per-node, rank-local ODE over Δt):
   - **Fringe favourability** `Φ = exp(−((wt − d0)/w)²)` — a Gaussian band centred on the
     mean capillary-fringe depth `d0` with half-width `w` (both YAML). Φ→1 when the
     surface sits at the fringe, →0 far above/below the table.
   - **Weathering supply** `Ψ` — the solute-supply rate feeding precipitation. Two modes,
     selected by the optional `weathering:` sub-block (§3a); both plug in at the *same* place:
     - **proxy (default)** `Ψ = clip(rain − evap, 0, ·)^p · arrhenius(T_annual)` — a
       temperature/water proxy (generic; emulates the humid-Fe vs arid-carbonate/silica
       contrast via `p` and the Arrhenius scaling). No new inputs.
     - **explicit rate (Level A, opt-in)** `Ψ = W` — a physically-based chemical-weathering
       rate driven by the groundwater recharge `R`, temperature, and (when `soilSPL` is on)
       regolith thickness. See §3a for the law and requirements.
   - **Formation** `dduriH/dt = k_form · Φ · Ψ · (1 − duriH/duriH_max)` (self-limiting to
     a maximum crust thickness). Induration degree `duriF = duriH / duriH_max ∈ [0,1]`.
   - **Breakdown / exhumation**: when the surface incises into or below the crust
     (`cumED < 0` stripping the top, or the crust emerges above the fringe into the
     dissolution/erosion zone), `duriH` is reduced by `k_break·|incision|` and by a slow
     `k_decay·(1−Φ)` disequilibrium term. A crust fully eroded through resets to 0.
6. **Armor K.** `duriF` feeds the erodibility hook (§5). No elevation change here —
   duricrust modifies *rate*, not geometry, so flow/routing are unchanged this step
   (like dual-lithology deposition being composition-only).
7. **Baseflow closure (opt-in).** Seepage discharge `Q_seep = Σ_owned (R − ΔS)`
   accounted at the seepage nodes so total river discharge stays `≈ rain − evap`
   over the quasi-steady step (`Allreduce`d budget). **As built:** the discharge
   is stored in `self.baseflowL` (conserved diagnostic, distributed over owned
   seepage nodes by cell area, `Σ baseflow ≈ Σ recharge` at steady state, written
   as the `baseflow` output) **and re-injected into the surface-flow source** —
   `applyForces` builds `bL = rain·A − recharge·A + baseflow`, so the infiltrated
   recharge leaves surface runoff and returns at the seepage nodes (net-neutral
   globally, rivers physically baseflow-fed, mirroring ice `iceMeltRiverL`).
   Applied at the single per-step `bL` reset so the two per-step
   `flowAccumulation` calls do not double-count the (non-idempotent) subtraction.
8. **Sync.** `localToGlobal` on `head`, `duriH` before the next collective (erosion).

**Water sources — rainfall and lakes (complementary roles).** The table is fed by *both*, but
they enter differently:
- **Rainfall** is the distributed **source**: `R = f·(rain − evap)` drives the head up on subaerial
  land (step 1).
- **Lakes / rivers / sea** are **fixed-head boundaries** (`h = z`, step 2); `R = 0` there (the head
  is pinned, not recharged). This makes the exchange **two-way and directional** — a lake ringed by
  lower head **leaks into** the aquifer (recharge), one ringed by higher head **receives** groundwater
  discharge (gaining lake). The head solve gets the flux sign automatically.

**Known limitation — lake ↔ aquifer *volume* coupling.** Lakes are fixed-head BCs (standard for
regional groundwater), so the water table responds to lakes correctly and rainfall recharges it
correctly — **but the lake's own volumetric budget** (goSPL's fill / evaporation / spill in
`flowplex._potentialLakeEvap` / `_distributeDownstream`) **is not yet reconciled with the
across-bed groundwater flux.** If the aquifer leaks into or drains a lake, the lake's fill/spill
budget does not see that exchange. Acceptable and conventional for a first version; full coupling
(lake leakage debits the lake, groundwater discharge fills it — the lake analogue of the river
baseflow closure, step 7) is an open decision (§15).

### 3a. Weathering-supply coupling (`Ψ`) — proxy vs explicit rate

The solute supply `Ψ` in the formation term (step 5) has two modes; both plug in at the same
point, so the rest of the pipeline is unchanged. Selected by an optional `weathering:` sub-block
inside `duricrust:` (absent ⇒ proxy).

**Proxy (default).** `Ψ = clip(rain − evap, 0, ·)^p · arrhenius(T_annual)`. No new inputs, no
conservation claim — a climate/temperature stand-in for solute availability.

**Level A — explicit chemical-weathering rate (opt-in, `weathering: mode: rate`).** Replace the
proxy with a physically-based rate `W`, still precipitating in-situ at the fringe (no solute
transport, no mass debit — supply-only, like the proxy but grounded in physics). The key enabling
input is the **recharge `R` the water-table solve already computes** — the water flux through the
weathering zone is the first-order control that the proxy lacks. Supported law (generic, one
tunable form):

```
W = R · C_eq · (1 − exp(−Dw / (R·L)))          # Maher–Chamberlain: kinetic × thermodynamic
    · arrhenius(T_annual)                       # optional extra Arrhenius (weather_Ea)
    · weatherability(lithology)                 # per-node/-lithology multiplier (from dual comp.)
```

with the two limits recovered by parameter choice: `Dw ≫ R·L` ⇒ kinetic/supply-limited
`W ≈ C_eq·Dw/L`; `Dw ≪ R·L` ⇒ thermodynamic/transport-limited `W ≈ R·C_eq`. `L` is the flow-path /
regolith length: **when `soilSPL` is on**, `L = Lsoil` (regolith residence controls the thermostat)
and `W` is additionally capped by the regolith production supply; **when off**, `L` is a prescribed
scale and `W` falls back to the recharge/kinetic form on bare rock. **Cheapest variant**
(`weathering: mode: prodsoil`): reuse `soilSPL.prodSoil` (already a temperature-scaled production
rate) directly as `W`, scaled by water availability — the "chemical ∝ physical weathering"
congruency, essentially free when soil is tracked.

```yaml
    duricrust:
        # ... formation/armor keys ...
        weathering:
            mode: proxy          # proxy (default) | rate | prodsoil
            C_eq: 1.0            # thermodynamic solute ceiling (rate mode)
            Dw: 1.0              # kinetic length scale (rate mode)
            path_length: 20.0    # L when soilSPL off (m); uses Lsoil when on
            weather_Ea: 0.0      # optional Arrhenius activation energy (0 ⇒ off)
            weatherability: 1.0  # scalar | map | per-lithology multiplier
```

**Level A is supply-only and NOT mass-conservative** — it does not remove dissolved solid from a
source pool. Making the crust mass a closed geochemical budget (dissolve → transport solute along
`q = −T∇h` → precipitate at the fringe → export via baseflow, with conservation guards) is **Level
B**, a separate geochemical solute-transport module scoped in §15 — the ingredients exist
(groundwater flux, FV advection kernels, per-class strata bookkeeping) but it is a major, distinct
feature with its own design doc.

---

## 4. State variables (all gated on `gwOn`)

**PETSc land** (persistent, halo-synced, in `destroy_DMPlex`, written/restored for restart):

| Vec | Meaning | Units |
|---|---|---|
| `self.headL` / `self.headG` | water-table head `h` (elevation of saturated surface) | m |
| `self.duriHL` / `self.duriHG` | duricrust thickness `duriH` | m |
| `self.rechargeL` | net recharge `R` this step (diagnostic/output) | m/yr |
| `self.baseflowL` | seepage return to rivers (opt-in) | m³/yr |

**Numpy land** (rank-local, no halo): `self.wtDepth` (`z − h`, m), `self.duriF`
(`duriH/duriH_max`, 0..1), `self.duriKarmor` (the K multiplier ≤ 1), `self.gwSeepIDs`
(Dirichlet seepage node indices). Cached operator/solver: `self._gwMat`, `self._ksp_gw`
(both in `destroy_DMPlex`).

**Stratigraphy (only when `gwOn and stratNb>0`):** `self.stratDuri[node, layer]` — per-layer
induration degree ∈ [0,1] (§9). Allocated with the other strata arrays in `_STRAMesh`, advected
with `stratHf`, written/restored in the stratal HDF5; **not** thickness-rescaled by compaction
(it is an intensive property, like `phiS`).

**Restart:** `head`, `duriH` are model memory (the crust integrates over My) and MUST
survive restart — written to the output HDF5 and restored in `outmesh.py` like `cumED`.

---

## 5. Erodibility armoring — the single hook

All three eroders funnel `K` through `stratplex._surfaceLithoK()` (`sed/stratplex.py:544`),
consumed at `SPL.py:62`, `nlSPL.py:252`, `soilSPL.py:277`. Duricrust armoring is a
**multiplicative factor injected there**, so it modulates SPL / nlSPL / soilSPL with
**no branching in the eroders**:

```python
def _surfaceArmoringK(self):
    if not self.gwOn:
        return 1.0                       # off ⇒ no-op ⇒ byte-identical
    # duriF ∈ [0,1]; armor_max ∈ [0,1) is the max erodibility reduction
    return 1.0 - self.armor_max * self.duriF
```

and in `_surfaceLithoK` the return becomes `(fc + ff·fine_k_factor) · _surfaceArmoringK()`.
`duriF=0` ⇒ factor 1.0 ⇒ unchanged. A fully indurated crust (`duriF=1`) cuts K by
`armor_max` (e.g. 0.9 ⇒ 10× more resistant), producing relief inversion. Optionally the
same `duriF` scales hillslope diffusivity `Cd` via `_surfaceLithoD` (crusts resist creep
too) — off by default, one YAML flag.

**No change to the `Eb`/`EbLocal` sign convention or the rcvID snapshot contract** —
armoring only rescales an existing coefficient.

---

## 6. YAML opt-in

```yaml
groundwater:
    Ksat: 1.0            # hydraulic conductivity K_h (m/yr), scalar | map | per-lithology
    specific_yield: 0.1  # S
    aquifer_base: 50.0   # z_bed depth below surface (m): scalar | map | `from_soil`
    bedrock_depth: 0.0   # d_bedrock: permeable weathered/fractured-rock depth below lHbed
                         #   (m), only used when aquifer_base: from_soil
    min_sat_thickness: 1.0    # b_min floor so transmissivity T > 0 near the base (m)
    infiltration: 0.3    # f_infil: fraction of (rain − evap) that recharges
    conserve_baseflow: True   # return seepage to rivers (re-injected into bL)
    lake_exchange: False # opt-in lake ↔ aquifer volume coupling (§15)
    subglacial_recharge: 0.0  # fraction of glacial meltwater infiltrating under ice
    fine_infil_factor: 1.0    # f_infil multiplier for the fine end-member (dual litho)
    infil_slope_ref: 0.0      # slope reference for f/(1+slope/ref) (0 = off)
    picard_its: 3
    seepage_passes: 4

    duricrust:
        form_rate: 1.0e-4     # k_form (m/yr at Φ=Ψ=1)
        max_thickness: 5.0    # duriH_max (m)
        fringe_depth: 3.0     # d0 — centre of the capillary fringe below surface (m)
        fringe_width: 2.0     # w — Gaussian half-width (m)
        supply_exp: 1.0       # p on (rain − evap)
        weather_Ea: 0.0       # Arrhenius activation energy (0 ⇒ no T dependence)
        armor_max: 0.9        # max fractional K reduction (0..1)
        armor_diffusion: False  # also armor hillslope Cd
        break_rate: 1.0       # k_break per unit incision
        decay_rate: 1.0e-6    # k_decay disequilibrium (1/yr)
```

Parsed by `_readGroundwater` → `_extraGroundwater` in `tools/inputparser.py`, following
the `_readIce`/`_extraIce` continuation pattern (`AGENTS.md` "`_extra*` are mandatory
continuations"). Missing block ⇒ `self.gwOn = False` and every parameter inert.
`Ksat`/`aquifer_base`/`infiltration` accept scalar, per-vertex map (npz), or a
`climate`-style time series (paleo-hydrology), resolved in `applyForces`.

**`aquifer_base: from_soil`** ties the aquifer floor to the bedrock elevation `lHbed`:
`z_bed = lHbed − bedrock_depth` (see §8 for the physics and the soil↔groundwater↔duricrust
feedback). It **requires `soilSPL`** (`lHbed` exists only when soil is tracked) — the parser
raises if `from_soil` is set without a `soil:` block. On soil-free / bare-bedrock cells
(`Lsoil ≈ 0`) it falls back to the prescribed depth, and the `min_sat_thickness` floor keeps the
transmissivity `T = K_h·max(h − z_bed, b_min)` positive everywhere.

---

## 7. New outputs

Written to the mesh HDF5 + `gospl.xdmf` (following `HOW_TO_ADD_OUTPUT.md`, registered in
`destroy_DMPlex` and the XDMF writer). Emitted only when `gwOn`:

| Field | XDMF name | Units | Meaning |
|---|---|---|---|
| water-table head | `wtable` | m | saturated-surface elevation `h` |
| water-table depth | `wtdepth` | m | `z − h`, depth below surface (the duricrust driver) |
| duricrust thickness | `duricrust` | m | `duriH` |
| induration degree | `induration` | – | `duriF ∈ [0,1]` |
| armoring factor | `Karmor` | – | effective K multiplier `1 − armor_max·duriF` (≤1) |
| net recharge | `recharge` | m/yr | `R` (diagnostic; verifies the climate → hydrology coupling) |
| seepage / baseflow | `baseflow` | m³/yr | returned to rivers (only if `conserve_baseflow`) |

`EDrate`/`cumED` already reflect the armored erosion (armoring rescales K in-place), so no
new erosion-rate field is needed. Post-processing: `wtdepth` + `duricrust` map directly in
the existing `gospl-grid` NetCDF export and `gospl-section` stratigraphic tools; a crust
recorded in the strata (see §9) is viewable as a per-layer `induration` property in `gospl-strata-volume`.

---

## 8. Soil-production dependence — and soil-free cells (the key question)

**This design assumes the Option-2.5 soil model** (`DESIGN_SOIL_REGOLITH.md`): `Lsoil` is the
*weathering-produced regolith* (subaerial, on bedrock), deposited sediment lives in the
**stratigraphy** (not routed into soil), a few *lumped profile scalars* (weathering-front depth,
weathering degree) ride the regolith, and subaerial exposure is gated by the lake/sea mask. The
duricrust couples to that model — and, crucially, **degrades gracefully** so it still runs on any
config.

**Duricrust does NOT require soil production (`soilSPL`) to be on.** The crust is a property of
the *near-surface material*, not of the tracked regolith specifically. The water-table + duricrust
mixin stays independent of `cptSoil`. Host medium / supply reference and where the fringe depth is
measured, by configuration:

| Run config | Host / weathering-supply reference | Fringe from |
|---|---|---|
| **regolith ON** (`cptSoil`) | the weathering regolith `Lsoil` over bedrock `lHbed` (Option-2.5 mantle) — the crust indurates the regolith; formation supply is regolith-limited (`min(k_form·Φ·Ψ, regolith rate)`) | surface `z` |
| **regolith OFF / bare bedrock** (`Lsoil ≈ 0`) | the surface material in-situ — bedrock, or the deposited-sediment top from the **stratigraphy**. No regolith reference, so supply falls back to the climate/temperature proxy `Ψ` (or the Level-A rate, §3a) | surface `z` |

So **where there is no regolith, the duricrust still forms** — it uses the `Ψ` proxy / Level-A rate
on the exposed surface (silcrete in bedrock, calcrete in sediment — generically via `Ψ`+`Φ`); it is
physically reasonable that duricrusts form both in regolith and by in-situ replacement of bedrock.

**Subaerial gate (Option-2.5).** Formation — like soil production — operates **only where the
surface is subaerial**: **not marine (`seaID`) and not ponded** (a continental lake, `pitIDs≥0`
with `lFill>hl`; `DESIGN_SOIL_REGOLITH.md` §3). A pit/lake deposit is **subaqueous soft sediment
until it fills to its spillover**; once emergent, regolith and the duricrust fringe begin operating
on it. The water table added here **generalizes** the gate to the seepage condition `h ≥ z`
(wetlands, near-surface table), but the first-order lake/sea mask exists already.

The water-table solve itself is independent of soil — it uses `z_bed` from the `aquifer_base`
parameter, with the `lHbed` tie an opt-in refinement described next.

**Shipped decision — both, always.** The feature ships **soil-independent** so it runs on any
mesh/config (bare bedrock, sediment-only, or no `soilSPL` at all), **and** it automatically
**couples to `soilSPL`** wherever soil is tracked. The same formation code path handles both:

```python
supply = k_form * Phi * Psi
if self.cptSoil:                      # soil ON  → regolith-supply-limited
    supply = np.minimum(supply, regolith_supply_rate(self.Lsoil))
duriH += supply * (1.0 - duriH/duriH_max) * dt
```

- `cptSoil` **on**: formation is limited by available regolith, and `aquifer_base` **may** be
  tied to `lHbed` (opt-in `aquifer_base: from_soil`) so the aquifer sits in the weathered
  regolith rather than at a prescribed depth.
- `cptSoil` **off**, or any cell with `Lsoil ≈ 0`: the `min(...)` limiter is skipped, supply
  falls back to the climate/temperature proxy `Ψ`, and `aquifer_base` uses the prescribed value.

There is **no configuration in which duricrust is unavailable** — the soil coupling is an
automatic refinement, never a prerequisite. Guarded by `test_duricrust_soilfree` (forms with
`cptSoil=False`) and a soil-on run.

**Aquifer floor tied to bedrock (`aquifer_base: from_soil`).** When regolith is tracked, the
aquifer base can follow the bedrock elevation: `z_bed = lHbed − bedrock_depth` (with `bedrock_depth`
the permeable weathered/fractured-rock zone below the regolith, and a `min_sat_thickness` floor so
`T = K_h·max(h − z_bed, b_min) > 0`). This is the physically apt *"permeable regolith over
impermeable bedrock"* model of cratonic/laterite terrains. **Under Option-2.5 this is now
well-defined**: `lHbed` is the base of the *weathering mantle* (not a fill-inflated surface), so the
depocenter pathology that made `lHbed` unusable under the old lumped-soil model is resolved (that
was `DESIGN_SOIL_REGOLITH.md` §2 item 4). **In depositional basins the aquifer base comes from the
stratigraphy** (the porous fill *is* the aquifer), not from `lHbed`; bare-rock / soil-off cells fall
back to the prescribed `z − aquifer_base`. It closes a genuine (slow, explicit, stable) feedback
loop:

```
soil production → lowers lHbed → deepens aquifer → shifts water-table / fringe depth
   → changes duricrust formation → armors K → changes erosion & soil exposure → soil production …
```

Deep weathering literally deepens the aquifer (right for laterite profiles). Every step in the
loop is slow (10⁴–10⁶ yr) and the coupling is **explicit/sequential within a step** (soil →
groundwater → duricrust → erosion, each reading the previous step's state), so there is no stiff
intra-step feedback — the same stability argument as goSPL's other explicit couplings.

**Stratigraphic coupling** (when `stratNb>0`) — how the armoring becomes part of the rock
record and evolves with erosion/deposition — is described in §9.

---

## 9. Stratigraphic integration — the erodibility record

Direct answers to the three questions: **is the armoring recorded in the stratigraphy?** yes,
when `stratNb>0`. **Does it vary with time?** yes. **Does it evolve with erosion/deposition?**
yes — burial preserves it and exhumation re-activates it. Details below.

### Two erodibility concepts, kept distinct
- **`stratK`** (existing): the **depositional** erodibility multiplier of each layer — fixed when
  sediment is deposited, exposed via the `erodeStrat` forward-fill (`stratplex.py:804`).
- **`stratDuri`** (new): the **diagenetic induration** degree ∈ [0,1] of each layer — an in-situ,
  post-depositional overprint from crust formation. **Not** a depositional property.

They **multiply** at the exposed surface. The effective erodibility of the current top layer is

```
K_eff = K · surfLithoK(stratK) · (1 − armor_max·stratDuri_top)
```

so a layer can be intrinsically soft (high `stratK`) yet crusted (high `stratDuri`) — the crust
dominates. Keeping the two separate avoids mutating a deposition-time property post-hoc, lets
breakdown/dissolution *reverse* the overprint without corrupting the depositional-K record, and
keeps the armoring hook (§5) a clean multiplier. (Baking armoring straight into `stratK` was
considered and rejected for exactly these reasons.)

### How `stratDuri` evolves — time, deposition, erosion
Per step, gated on `gwOn and stratNb>0`:

- **Formation (varies with time).** After `_updateDuricrust`, the induration is written down into
  the near-surface stratigraphy: `stratDuri[node, top]` is raised toward the live `duriF`
  (`_recordInduration`). On a stable, non-eroding surface the crust thickens over 10⁴–10⁶ yr → the
  near-surface layer's `stratDuri` **grows with time**. **As built:** the write-down targets the
  **top non-empty layer** of each column (found as in `_surfaceComposition`); distributing it over a
  multi-layer crust/fringe *depth range* is a possible refinement, unnecessary for the exhumation
  behaviour below (which reads the top non-empty layer).
- **Deposition (burial → preservation).** `deposeStrat` adds a new top layer with `stratDuri = 0`
  (fresh, uncemented sediment). The previously indurated layer keeps its `stratDuri` and is now
  **buried and preserved** — a relict crust locked into the record.
- **Erosion (exhumation → re-armor).** `erodeStrat` strips layers top-down; a crust eroded through
  is removed with its layer. When erosion exposes a **previously buried** indurated layer, its
  preserved `stratDuri` becomes the new surface armoring — the live `duriF` is re-seeded to
  `max(decayed duriF, stratDuri[new top])`, so the relict crust **re-armors** and resists further
  incision. This is the multi-cycle, stacked-duricrust / relief-inversion behaviour of cratonic
  (e.g. Australian laterite) landscapes.
- **Advection / compaction.** `stratDuri` is an **intensive** per-layer property (like `phiS`/
  `phiF`), so it **advects with the pile** (piggybacks the `stratHf` advection in `stratalRecord`)
  and is **compaction-neutral**: burial reduces thickness/porosity but not the induration degree,
  so `getCompaction` does **not** rescale it (contrast `stratP`, a thickness partition, which is
  rescaled).

### Surface (active) vs archive
The live per-node `duriF` is the **active** armoring applied this step; `stratDuri` is the
**archived** per-layer value. They sync at the top layer — formation writes `duriF` *down* into the
top layers; exhumation reads the exposed layer *up* into `duriF`. With `stratNb == 0` there is **no
archive**: `duriF` is surface-only, aggradation resets it to 0 (fresh material buries it with no
memory), and a re-incised surface re-forms the crust from scratch. That is the accepted
reduced-fidelity mode when stratigraphy is off.

### No geometry change
Recording induration is **composition-only** — exactly like dual-lithology deposition and
provenance, it never edits layer thickness or elevation, so flow / routing / mass balance are
untouched the step the crust forms or is recorded. The landscape response (relief inversion)
emerges through the *normal* erosion pathway, because armored cells simply erode more slowly.

### Output
`stratDuri` is written to the stratal HDF5 (like `stratHf`/`stratP`) and restored on restart, and
is exposed as a new per-layer `induration` field in `gospl-strata-volume` (alongside
`porosity`/lithology/provenance) — so a `gospl-section` cross-section shows buried and exhumed
crusts directly.

---

## 10. Compatibility with dual-lithology & provenance

**Verdict: fully compatible with each, and with all three enabled together.** The design was
built to compose — armoring is a multiplicative factor at the *same* erodibility hook the other
two already use, and `stratDuri` is one more *independent, intensive* per-layer field alongside
the mass/composition fields they add. Nothing about duricrust perturbs the conservation invariants
they are guarded by. (Confirmed against `stratplex.py`: `phiF`/`phiS` are per-layer intensive
`(lpoints, stratNb)` arrays that already ride erode/deposit/advect/compact as passengers —
`stratDuri` slots in identically.)

### Shared stratigraphic fields are orthogonal by construction
Each feature adds per-layer strata state; they do not alias:

| Feature | Per-layer field(s) | Kind | Compaction | Advection |
|---|---|---|---|---|
| dual-lithology | `stratHf`, `phiF` | fine **mass** + porosity | per-fraction (`_depthPorosityDual`) | `stratHf` advected (2nd `strataonesed`) |
| provenance | `stratP[·,·,class]` | class **mass** partition (Σ = `stratH`) | **rescaled** to keep Σ = `stratH` | per-class, renormalised |
| **duricrust** | `stratDuri` | **intensive** induration ∈ [0,1] | **neutral** (not rescaled) | like `phiF` (intensive, *not* renormalised) |

Because `stratDuri` is **intensive and carries no volume** (like porosity — unlike the masses
`stratHf`/`stratP`), adding it **cannot break** `test_dual_fine_conservation`,
`test_provenance_conservation`, or `test_mass_conservation`. `erodeStrat`/`deposeStrat` thread it
as a passenger: a removed layer loses its induration; a *partially* eroded layer keeps its degree
(intensive); a newly deposited layer starts at 0 — exactly how `phiF`/`phiS` already ride those
routines. This is the safest kind of field to add (a conserved one would be riskier).

### Erodibility hook composes multiplicatively (duricrust × dual-lithology)
The single `_surfaceLithoK` (`stratplex.py:544`) already carries the dual-lithology term; the
duricrust factor multiplies onto it:

```
surfLithoK = (fc + ff·fine_k_factor)  ·  (1 − armor_max·duriF)
             └─ lithology (coarse/fine) ─┘   └─ diagenetic armor ─┘
```

Orthogonal controls — lithology sets the *intrinsic* rock strength, the crust sets the *diagenetic*
overprint. A fine (mud) layer cemented into a ferricrete hardpan correctly reads
`fine_k_factor · (1 − armor_max·duriF)` (intrinsically weak, but armored) — the physically right
"hardpan over mud" behaviour. Same multiplicative composition for the optional diffusivity armor
via `_surfaceLithoD` (`fc + ff·fine_diff_factor` × the crust factor).

### Provenance rides the armored erosion automatically (duricrust × provenance)
Provenance is a **passive label** (`AGENTS.md`: "no K/D/porosity/sorting feedback; it rides the
total-sediment routing"). Duricrust changes *erosion rate* through K, so an armored source region
erodes less and contributes proportionally less to downstream `stratP` — the provenance signal
reflects the armoring **with no extra code**. What v1 deliberately does **not** do is give the
*duricrust itself* a clastic provenance class: the crust is precipitated from groundwater solutes
(a diagenetic, not detrital, product), so it has no clastic source. Attributing the Fe/Si/CaCO₃
solute source would need a solute-transport tracer — a possible future extension, out of scope here.

### All three at once
The three per-layer fields advect, compact, erode and deposit independently; the erodibility hook
composes; provenance stays passive. Memory is modest (`stratDuri` is one `(lpoints, stratNb)` array,
same size as `stratHf`). The only new bookkeeping is adding `stratDuri` to the layer-carrying loops
in `erodeStrat`/`deposeStrat`/`stratalRecord`/`getCompaction` — the same passenger pattern already
used, and the conservation guards stay green because `stratDuri` carries no mass. A combined
`test_duricrust_dual_provenance` (all three on) pins that the conservation tests still hold and the
armoring composes.

### Optional couplings (enhancements, all opt-in, default off)
- **Composition-dependent formation** — let the surface fine fraction `ff` (dual) modulate `k_form`
  (clay-rich hosts cement differently). One extra multiplier inside `Ψ`; default 1.0 (no-op).
- **Crust preserves fines** — already automatic: a cemented cap armors whatever it caps, including a
  fine-rich layer, through the shared K hook. No code.
- **Solute provenance** — attribute the duricrust's chemical source (future; needs solute routing).

---

## 11. Parallelisation strategy

The whole capability is built on machinery `AGENTS.md` already certifies as
partition-safe; the design adds no new collective-gating hazards.

- **Head solve = standard elliptic FV solve on the DMPlex.** Same partition, halo
  exchange, `lgmap`, and owned-rows-only assembly (`self.glIDs`) as the hillslope/marine
  diffusion operators — one of the three *safe* assembly patterns in `AGENTS.md` §"#2
  partition-dependence" (additive FV-Laplacian, `ADD_VALUES`). Partition-exact head to
  KSP tolerance; the cached `gw_` KSP uses **fgmres + hypre BoomerAMG** (algebraic
  multigrid — the stiff elliptic operator needs it; block-Jacobi/ILU does not converge),
  env-overridable via the `gw_` prefix.
- **Seepage Dirichlet set is partition-invariant.** It is derived from the drainage
  network (rivers/lakes) + `seaID` + `outletIDs`, all built from the partition-invariant
  `locIDs`-keyed drainage arrays (the mechanism-#2 fixes). Applied by `zeroRows`
  (Dirichlet), exactly like the clamped FEM-flexure BC.
- **Picard + seepage outer loops are collective-safe.** Every iteration runs on all ranks;
  the convergence/continue tests break on an `Allreduce`d scalar (head-change norm; new-
  seepage `LOR` flag) — identical on every rank, so no rank-local `.any()` gates a
  collective (the #1 deadlock rule). Cap the Picard (`picard_its`) and seepage
  (`seepage_passes`) counts as a backstop, like `_cascade_max_it`.
- **Operator rebuild.** `T(h)` changes every step, so the operator is rebuilt each step via
  the single-pass `_assembleDiffMatCSR` (no rank-local rebuild gate to reduce — it is
  unconditional). The KSP object is cached and reused; only PCSetUp repeats.
- **Duricrust update is embarrassingly parallel** — a per-node rank-local ODE on
  `wtDepth`, `rain`, `evap`, `cumED` (all already local), then a single `localToGlobal`
  to refresh `duriH` halos before erosion reads it. No reductions needed except the
  optional `conserve_baseflow` budget (`Allreduce` sum over owned nodes, like
  `evapLoss` / ice `melt_conserve`).
- **Cost & scaling.** One extra cached diffusion solve (× 2–3 Picard × 2–4 seepage passes,
  but each pass warm-starts from the last so later passes are cheap) per step — the same
  order as the existing hillslope solve, and it inherits goSPL's strong-scaling behaviour.
  No new serial bottleneck (unlike the pit-graph or global flexure).
- **Restart / determinism.** `head`/`duriH` in `destroy_DMPlex` + the output/restore path;
  the solve is deterministic given the (invariant) seepage set, so np=1 vs np=N agree to
  the usual KSP-noise floor guarded by `test_parallel_correctness`.

---

## 12. AGENTS.md contracts to honor

- **Opt-in, byte-identical when off** — `gwOn=False` gates every path; guard with a bitwise
  `test_groundwater_opt_in` (a non-`groundwater` run reproduces current output exactly).
- **`destroy_DMPlex`** — register `headL/headG`, `duriHL/duriHG`, `rechargeL`, `baseflowL`,
  `_gwMat`, `_ksp_gw` (the `unstructuredmesh.py:1096-1206` list).
- **MPI #1 (deadlock)** — every Picard/seepage continue-flag is `Allreduce`d before gating a
  collective re-solve.
- **MPI #2 (assembly)** — head operator sets **owned rows only** (`self.glIDs`), additive
  FV-Laplacian with `ADD_VALUES`.
- **KSP lifecycle** — `_ksp_gw` is CACHED (hot path, every step), created lazily, never
  `destroy()`d mid-run; in the destroy list.
- **Scratch vectors** — the head/duricrust update may use `tmp`/`tmpL`/`tmp1` but MUST leave
  them defined; document which in each method docstring. `head`/`duriH` are persistent
  state, never scratch.
- **`_extra*` continuation** — `_readGroundwater`→`_extraGroundwater` sets all attributes
  other modules read; never delete/rename without the call chain.
- **Mixin init order** — insert `_GWMesh.__init__` **after `_FAMesh` (`model.py:254`)** (needs
  the flow graph + diffusion machinery) and before the eroders, so `duriF`/the armoring hook
  exist when `_surfaceLithoK` is first called.
- **Constants** — the fringe/supply/armor literals live in the YAML block; any hard sentinel
  (e.g. a minimum saturated thickness floor to keep `T>0`) goes in `tools/constants.py`.

---

## 13. Phased implementation plan (branch per phase, PR into `dev`)

**All phases below are DONE** on `feat/watertable-duricrust` (7 feature commits;
full `tests/` 139 passed, serial + np=2). The status/notes are inline per row.

| Phase | Deliverable | Guard test |
|---|---|---|
| −1 | **DONE (soil PR #482).** Option-2.5 consistency fixes — subaerial lake/sea gate + submarine coherence + ice freeze-inert (`DESIGN_SOIL_REGOLITH.md` §5). | `test_soil_subaerial_gate` |
| 0 | **DONE.** `_readGroundwater` parser + `gwOn`/`duriOn` flags + `_GWMesh` state alloc (head/duriH/recharge/baseflow Vecs, per-node state, cached `_gwMat`/`_ksp_gw`) + `destroy_DMPlex`; init after `_FAMesh`. | `test_groundwater_opt_in` (off ⇒ inert; on ⇒ state + head seeded) |
| 1 | **DONE.** Recharge `R = f·max(0, rain−evap)` from existing forcing, zeroed under **water AND ice** (§3); `f_infil` scalar **or** per-vertex map; `recharge` output. | `test_groundwater_recharge` (humid⇒f·(P−E); arid⇒0; under-water/ice⇒0; per-vertex f honoured) |
| 2 | **DONE + analytically validated.** Implicit head solve `_solveHead`: `(I + (Δt/S)·L(T))h = h_old + (Δt/S)·R` via `jacobiancoeff`/`_assembleDiffMatCSR` (like the marine Picard), Picard on `T=Kh·max(h−z_bed, b_min)`, **two free boundaries** — seepage `h≤z` (Dirichlet `zeroRowsLocal` at `seaID`∪ponded-lake∪`outletIDs`, + clip-and-discover, `Allreduce`'d new-seepage break) and dry-aquifer floor `h≥z_bed`; cached **fgmres + hypre BoomerAMG** `gw_` KSP (block-Jacobi/ILU does **not** converge on this stiff 2-D elliptic operator — it stalls at `DIVERGED_ITS` and the garbage iterate drifts to the surface; AMG is required, cf. the flexure biharmonic); `aquifer_base` prescribed — scalar **or per-vertex map** (`from_soil` = Phase 5); `wtable`/`wtdepth` outputs. | `test_watertable_solve` (bounded `z_bed≤h≤z`, finite, non-trivial); `test_watertable_steady` (repeated solves contract onto the quasi-steady Dupuit fixed point); `test_watertable_parallel` (np=1-vs-2 agree within the partition-drift floor). **`benchmarks/test_dupuit.py`: full analytic Dupuit-parabola benchmark** — west-draining ramp, flat base via an `aquifer_base` map, `bc='wwwf'`; steady head matches `s²=s_d²+(R/K)(2Wx−x²)` to **RMSE 0.00 %, R²=1.0**. |
| 3 | **DONE.** Duricrust ODE (`_updateDuricrust`, rank-local): fringe favourability `Φ = exp(−((wt−d0)/w)²)`, weathering supply `Ψ` (`_weatheringSupply`: **proxy** default `max(0,P−E)^p·arrhenius`; opt-in Level-A **rate** Maher–Chamberlain `W=R·C_eq·(1−exp(−Dw/(R·L)))·…`; **prodsoil** reuse of `soilSPL.prodSoil` — both fall back to proxy when soil is off), self-limiting formation `+k_form·Φ·Ψ·(1−duriH/duriH_max)`, breakdown (per-step incision `z_last−z` strips the top; slow `k_decay·(1−Φ)` decay); optional Arrhenius (`_arrhenius`, reuses the soil tempMap, off when `weather_Ea=0`); writes `duriH`, `duriF=duriH/duriH_max`, `duriKarmor=1−armor_max·duriF`; `duricrust`/`induration` outputs (HDF5+XDMF). Soil-independent by default. | `test_duricrust_forms_at_fringe` (forms at `wt≈d0`, not away); `test_duricrust_soilfree` (forms with no `soil:` block); `test_duricrust_weathering_rate` (Level-A rate ↗ with `R`, 0 at `R=0`). |
| 4 | **DONE.** Armor hook `_surfaceArmoringK` (`sed/stratplex.py`) composed multiplicatively into `_surfaceLithoK` → reaches **all three eroders** (SPL/nlSPL/soilSPL, no branching) as `1 − armor_max·duriF`; scalar `1.0` no-op when off (byte-identical). Optional creep armoring in `_surfaceLithoD` behind `armor_diffusion`. Rate-only (no geometry change → routing untouched). | `test_duricrust_armors_K` (indurated cell K reduced by `armor_max`, erodes ≪ bare; scalar-1.0 no-op when off). |
| 5 | **DONE.** Baseflow conservation (opt-in `conserve_baseflow`, `_baseflowClosure`): seepage-return discharge `Q_seep = Σ(R·A) − ΔS/Δt` (`Allreduce`'d), distributed over owned seepage nodes by cell area into `self.baseflowL`, so `Σ baseflow ≈ Σ recharge` in the steady limit (`baseflow` output); re-injection into the surface-flow source is the next increment. Soil coupling: **regolith supply limiter** (`_regolithSupplyRate = prodSoil·rain` caps formation when `cptSoil`) and **`aquifer_base: from_soil`** (`_gwZbed`: `z_bed = lHbed − bedrock_depth`, `lHbed` now init'd in `soilSPL.__init__`). | `test_groundwater_baseflow_conserves` (Σ baseflow ≈ Σ recharge to 2 %); `test_groundwater_from_soil` (`z_bed = lHbed − d_bedrock`, head bounded); `test_duricrust_regolith_limited` (formation capped by `prodSoil·rain`). |
| 6 | **DONE.** Stratigraphic induration record `stratDuri` (`(lpoints, stratNb)`, §9), allocated when `gwOn and stratNb>0`. `_recordInduration` (after `_updateDuricrust`): **write-down** (record live `duriF` into the top non-empty layer) + **read-up** (an exhumed buried crust re-arms `duriF`/`duriKarmor`). Burial preserves it (`deposeStrat` fresh layers = 0); `erodeStrat` zeroes emptied layers (**no forward-fill** — 0 is a valid uncemented value); advected as an INTENSIVE field (extra `strataonesed`, like `phiS`), compaction-neutral; written/restored in the stratal HDF5; `induration` per-layer field in `gospl-strata-volume` (`stratamesh`). | `test_duricrust_strata_exhumation` (buried crust re-armors on re-exposure; full suite incl. dual-lithology/provenance/restart green — `getattr` guards for bare-`STRAMesh` unit tests). |
| 7 | **DONE.** Restart of `head`/`duriH` (model memory — written to the per-step HDF5, restored in `readData` like `cumED`/`soilH`; `wtDepth`/`duriF`/`duriKarmor` rebuilt); new `Karmor` output. **Documentation (§14):** `groundwater:` block + all keys in `surfproc.rst` (the convention home for process blocks — `climate:`/`ice:` live there too, not `inputfile.rst`); `outputs.rst` fields; `running.rst` (`--field induration`); new `tech_guide/groundwater.rst` (Dupuit–Boussinesq, implicit solve + why, seepage/Picard, fringe formation, `stratDuri`, compatibility) wired into `tech_guide/index.rst`; new `api_ref/gw_ref.rst` (autosummary + automethod) wired into `api_ref/index.rst`; `_readGroundwater` on `in_ref.rst`, `_surfaceArmoringK` on `stra_ref.rst`; `gwplex` module/class docstrings de-staled. `petsc4py`/`gospl._fortran` already mocked → autodoc imports cleanly. | `test_groundwater_restart` (head/duriH survive restart); RST underline/label lint clean; full `tests/` 139 passed. |

---

## 14. Documentation updates (ship with the feature)

Documentation is part of the deliverable (Phase 7), per the goSPL convention that every
user-facing feature updates the **input-file reference**, the **technical guide**, and the
**API reference**. Concretely:

### Input-file reference (`docs/user_guide/`)
- **`inputfile.rst` — a new `groundwater:` block section** documenting every key: `Ksat`,
  `specific_yield`, `aquifer_base` (scalar | map | `from_soil`), `infiltration`,
  `conserve_baseflow`, `picard_its`, `seepage_passes`, and the nested `duricrust:` sub-block
  (`form_rate`, `max_thickness`, `fringe_depth`, `fringe_width`, `supply_exp`, `weather_Ea`,
  `armor_max`, `armor_diffusion`, `break_rate`, `decay_rate`). State units, defaults, the opt-in
  semantics (absent ⇒ off ⇒ byte-identical), and the scalar / per-vertex-map / time-series forms.
  Mirror the existing `climate:` / `ice:` block layout.
- **`optfile1.rst` / `optfile2.rst`** — cross-reference from the forcing / erodibility option pages
  where those keys are indexed by convention.
- **`surfproc.rst`** — a "Groundwater & duricrust" subsection in the surface-processes narrative:
  the physical model (fringe formation → armoring → relief inversion) at user level, when to enable
  it, and the soil-independent / soil-coupled behaviour (§8).
- **`outputs.rst`** — document the new fields `wtable`, `wtdepth`, `duricrust`, `induration`,
  `Karmor`, `recharge`, `baseflow` (§7) and the per-layer `induration` in `gospl-strata-volume` (§9).
- **`running.rst`** — note `gospl-strata-volume --field induration`; no new CLI command (in-model).

### Technical guide (`docs/tech_guide/`)
- **A new `tech_guide/groundwater.rst` page** at the depth of `tech_guide/provenance.rst` /
  `strat.rst` / `ice.rst`: the Dupuit–Boussinesq formulation, the implicit backward-Euler
  discretisation and *why* (the τ_gw-vs-Δt argument, §2), the seepage free boundary + Picard
  treatment, the fringe-formation law, the stratigraphic induration record (§9), and the
  dual-lithology / provenance compatibility (§10). Wire it into `tech_guide/index.rst`
  (toctree + card). Cross-link this design doc, as the other tech pages link their `DESIGN_*.md`.

### API reference (`docs/api_ref/`)
- **A new autodoc page `api_ref/gw_ref.rst`** for `_GWMesh` (`gospl/flow/gwplex.py`) — `autoclass`
  plus **hand-maintained** `.. autosummary::` **and** `.. automethod::` lists (both, per the
  AGENTS.md autodoc invariants) for `updateGroundwater`, `_solveHead`, `_updateDuricrust`,
  `_surfaceArmoringK`, etc. Wire into `api_ref/index.rst` (card + toctree).
- **Extend existing API pages** where new public/private methods land: `_surfaceArmoringK` (and the
  `stratDuri` handling) on the `stratplex` page; `_readGroundwater`/`_extraGroundwater` on the
  `inputparser` page (like `_extraStrata`/`_extraProvenance`).
- **Keep the docs build green**: the new module must import under Sphinx autodoc — it reuses only
  already-mocked compiled/MPI deps (`petsc4py`, `gospl._fortran`), so no new `autodoc_mock_imports`
  entry is needed; do **not** mock the docs-installed scientific deps (`numpy`/`scipy`/…), per the
  AGENTS.md "Docs / Read the Docs" invariants.

### AGENTS.md + design docs
- Add a **Milestones** row (date, branch, one-paragraph change log) and update any AGENTS.md
  contract section whose invariant changes (the new cached `_ksp_gw`, the `destroy_DMPlex`
  additions, the new `stratDuri` strata field, the mixin init-order insertion).

---

## 15. Open decisions (defaults chosen, revisit on validation)

- **Soil-model dependency — DECIDED.** This design targets the **Option-2.5** soil model
  (`DESIGN_SOIL_REGOLITH.md`): `Lsoil` = weathering-only regolith, deposited sediment in the
  stratigraphy, subaerial lake/sea gate, lumped profile scalars, `soil_transition` recast as a
  smooth max-weathering-depth. The duricrust still degrades gracefully when soil/stratigraphy are
  off (§8). Implementation order: land the Option-2.5 *consistency fixes* (subaerial gate +
  submarine coherence) first — they de-risk the duricrust coupling and stand alone.
- **Fringe favourability shape** — Gaussian band vs a top-hat `[d0−w, d0+w]`. **As built:**
  Gaussian `exp(−((wt−d0)/w)²)` (smooth gradients, better for the per-node ODE); revisit only if a
  sharp fringe is wanted.
- **Aquifer base `z_bed`** — **prescribed `z − aquifer_base` (default)** OR **`lHbed − bedrock_depth`
  when `aquifer_base: from_soil`** (requires `soilSPL`; §8). The `from_soil` form is the physically
  apt "permeable regolith over impermeable bedrock" model for cratonic/laterite terrains and makes
  deep weathering deepen the aquifer; `bedrock_depth` adds the fractured-rock zone and the
  `min_sat_thickness` floor keeps `T>0`. Under Option-2.5 `lHbed` is the base of the weathering
  mantle (no depocenter inflation), and **in depositional basins the base comes from the
  stratigraphy** (the porous fill is the aquifer), not `lHbed`. Bare-rock / soil-off cells fall
  back to the prescribed depth, which stays the default so the feature runs standalone. **As built
  (DONE):** `_gwZbed` implements all three — scalar/map, `from_soil = lHbed − bedrock_depth`, and
  the basin case as the *deeper* of `lHbed − bedrock_depth` and `z − Σ(non-sentinel sediment)`.
- **Recharge refinements (DONE, opt-in, default off)** — `R = f_infil·max(0, rain − evap)` with
  three optional modulations (§3 step 1): **subglacial-meltwater recharge** (`subglacial_recharge`
  fraction of the ice model's `iceMeltRiverL` infiltrates under ice — the one recharge path the ice
  gate allows), **lithology** (`fine_infil_factor` scales `f_infil` by the exposed coarse fraction,
  dual lithology), and **slope** (`f/(1+slope/infil_slope_ref)` from a steepest-descent proxy). Each
  defaults to a no-op, so the baseline `f_infil` (scalar or per-vertex map) is unchanged.
- **Weathering supply `Ψ`** — three tiers, escalating cost (see §3a for the rate law + YAML):
  - **Proxy (default, shipped)** — climate/temperature stand-in; no new inputs, non-conservative.
  - **Level A (opt-in, DONE)** — explicit chemical-weathering *rate* `W(R, T, Lsoil,
    lithology)` driven by the groundwater recharge `R` (± the `prodsoil` congruency shortcut when
    `soilSPL` is on). Supply-only, still non-conservative. Shipped (`_weatheringSupply`,
    `weathering: mode: rate|prodsoil`); cheap because all inputs already exist.
  - **Level B (future, separate module — the main remaining item)** — conservative geochemistry:
    debit dissolved solid,
    transport solute along `q = −T∇h` (advection-reaction on the DMPlex, reusing the FV advection
    kernels), precipitate at the fringe, export via baseflow, with `Σ dissolved − precipitated −
    exported ≈ 0` guards. Comparable in scope to dual-lithology; needs its own design doc. The
    enabling pieces (groundwater flux, FV advection, per-class strata bookkeeping) already exist.
- **Lake ↔ aquifer volume coupling** — lakes/rivers/sea are **fixed-head** boundaries (`h = z`) in
  the head solve. **As built (opt-in `lake_exchange`, default off):** the lake's volumetric budget
  is now also debited/credited by the across-bed groundwater flux. `_lakeExchangeFlux` computes the
  signed per-node flux `∇·(T∇h)·A = −(L·h)·A` (>0 = aquifer discharges into the lake, <0 = lake
  leaks into the aquifer) after the head solve; `_distributeDownstream` aggregates it per lake
  (partition-invariant `group_by` + `Allreduce`) and feeds it into the fill budget `inV` at step 0 —
  gains added, leakage debited and clamped at the available water (mirrors the evaporation debit).
  Default off ⇒ the conventional fixed-head behaviour, byte-identical. (Lake **level** feedback onto
  the head boundary within the same step remains one-step-lagged, like all the explicit couplings.)
- **Armor of diffusion** — **As built (DONE):** off by default (SPL K only); `armor_diffusion`
  also scales the hillslope diffusivity `Cd` by `1 − armor_max·duriF` (`_surfaceLithoD`).
- **Do we need transient `head` at all, or steady each step?** **As built:** carried as state with
  backward-Euler (robust in both `τ_gw` regimes; analytic-Dupuit validated). If validation ever
  shows the equilibrium limit everywhere, a pure steady solve is a trivial simplification — still
  open, low priority.
- **Multi-layer formation depth range** (§9) — `_recordInduration` currently writes the crust into
  the **top non-empty layer** only. Distributing it over the crust/fringe *depth range* is the
  remaining formation refinement (unnecessary for the exhumation behaviour, which reads the top
  layer). Still open, low priority.
- **Solute-source provenance** (§11) — blocked on **Level B**: attributing the crust's chemical
  source needs the solute-transport tracer. Deferred with Level B.
