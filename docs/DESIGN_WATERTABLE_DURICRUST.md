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
| 8 | Soil dependence | **Ships both**: runs **soil-independent** on any config, AND automatically **couples to `soilSPL`** when soil is tracked (regolith supply limiter + optional `aquifer_base = lHbed`). See §8. |
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
   existing forcing arrays. `seaID` and closed-lake interiors get `R = 0` (their head
   is pinned; see step 3). `f_infil` may be a scalar, per-lithology, or a map.
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
   - Cached `gw_`-prefixed KSP (`_makeDiffusionKSP`); operator rebuilt per step (T varies),
     KSP object reused, PC not reused.
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
7. **Baseflow closure (opt-in).** Seepage discharge `Q_seep = Σ_owned (R − ΔS)` returned
   to the river network at seepage nodes so total river discharge stays `≈ rain − evap`
   over the quasi-steady step (mirrors ice `melt_conserve`; `Allreduce`d budget).
8. **Sync.** `localToGlobal` on `head`, `duriH` before the next collective (erosion).

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
    conserve_baseflow: True   # return seepage to rivers (river-discharge neutral)
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

**Duricrust does NOT require soil production (`soilSPL`) to be on.** The crust is a
property of the *near-surface material*, not of the tracked soil layer specifically. The
design deliberately keeps the water-table + duricrust mixin independent of `cptSoil`.

What the host medium and supply reference are, by run configuration:

| Run config | Host / weathering-supply reference | Fringe depth measured from |
|---|---|---|
| **soil production ON** | regolith thickness `Lsoil` over bedrock `lHbed` — the crust cements the regolith; formation supply is additionally *limited by available regolith* (`min(k_form·Φ·Ψ, regolith supply rate)`) | surface `z` (= `lHbed + Lsoil`) |
| **soil production OFF** (or bare-bedrock cells where `Lsoil ≈ 0`) | the **surface material itself** (bedrock or deposited sediment). No `Lsoil` reference exists, so supply falls back purely to the **weathering proxy `Ψ`** (climate + temperature), i.e. in-situ replacement/cementation of the top of bedrock | surface `z` (= `hLocal`) |

So **where there is no soil, the duricrust still forms** — it simply uses the climate/
temperature weathering proxy `Ψ` for solute supply instead of a regolith-thickness supply
limiter, and indurates the exposed surface (silcrete forming in bedrock, calcrete in
sediment, etc. — all captured generically by `Ψ` + the fringe favourability `Φ`). This is
physically reasonable: duricrusts form both in regolith and by in-situ replacement of
bedrock.

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

**Aquifer floor tied to bedrock (`aquifer_base: from_soil`).** When soil is tracked, the aquifer
base can follow the bedrock elevation: `z_bed = lHbed − bedrock_depth` (with `bedrock_depth` the
permeable weathered/fractured-rock zone below the regolith, and a `min_sat_thickness` floor so
`T = K_h·max(h − z_bed, b_min) > 0`). This is the physically apt *"permeable regolith over
impermeable bedrock"* model of cratonic/laterite terrains. Bare-rock / soil-off cells fall back to
the prescribed `z − aquifer_base`. It closes a genuine (slow, explicit, stable) feedback loop:

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
  the stratigraphic layers lying within the crust/fringe depth: `stratDuri[node, top layers]` is
  raised toward the live `duriF`. On a stable, non-eroding surface the crust thickens over
  10⁴–10⁶ yr → those near-surface layers' `stratDuri` **grows with time**.
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
  KSP tolerance; the cached `gw_` KSP uses `fgmres`/`bjacobi` or CG (env-overridable),
  the same class as the flow-accumulation solver.
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

| Phase | Deliverable | Guard test |
|---|---|---|
| 0 | `_readGroundwater`/`_extraGroundwater` parser + `gwOn` flag + state alloc + `destroy_DMPlex` | `test_groundwater_opt_in` (bitwise off) |
| 1 | Recharge `R = f·(rain−evap)` from existing forcing; `recharge` output | `test_groundwater_recharge` (arid⇒0, humid⇒f·(P−E)) |
| 2 | Implicit head solve (`_solveHead`): Picard `T(h)` + seepage clip, cached `gw_` KSP; `wtable`/`wtdepth` outputs | `test_watertable_solve` (analytic Dupuit hillslope; np=1-vs-2) |
| 3 | Duricrust ODE (`_updateDuricrust`): fringe Φ, supply Ψ (proxy default; opt-in Level-A explicit-rate `W` via `weathering:`, §3a), formation + breakdown; `duricrust`/`induration` outputs | `test_duricrust_forms_at_fringe`; `test_duricrust_soilfree` (forms with `cptSoil=False`); `test_duricrust_weathering_rate` (Level-A rate responds to `R`) |
| 4 | Armor hook `_surfaceArmoringK` into `_surfaceLithoK`; relief-inversion behaviour | `test_duricrust_armors_K` (indurated cell erodes ≪ bare) |
| 5 | Baseflow conservation (opt-in); soil coupling (regolith supply limiter + `aquifer_base=lHbed`) | `test_groundwater_baseflow_conserves` (Σ baseflow ≈ Σ recharge) |
| 6 | stratigraphic induration record (`stratDuri`, §9): write-on-formation, burial via `deposeStrat` (new layer = 0), exhumation re-armor via `erodeStrat`, advect like `stratHf`, compaction-neutral; `induration` field in `gospl-strata-volume` | `test_duricrust_strata_exhumation` (buried crust re-armors on re-exposure) |
| 7 | **Documentation (§14)**: `groundwater:` block in `inputfile.rst`, `surfproc.rst` + `outputs.rst` updates; new `tech_guide/groundwater.rst`; new `api_ref/gw_ref.rst` (+ `stratplex`/`inputparser` API pages); AGENTS.md milestone | docs build green (autodoc imports cleanly) |

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

- **Fringe favourability shape** — Gaussian band vs a top-hat `[d0−w, d0+w]`. Gaussian
  chosen for smooth gradients (better for the KSP-free per-node ODE); revisit if a sharp
  fringe is wanted.
- **Aquifer base `z_bed`** — **prescribed `z − aquifer_base` (default)** OR **`lHbed − bedrock_depth`
  when `aquifer_base: from_soil`** (requires `soilSPL`; §8). The `from_soil` form is the physically
  apt "permeable regolith over impermeable bedrock" model for cratonic/laterite terrains and makes
  deep weathering deepen the aquifer; `bedrock_depth` adds the fractured-rock zone and the
  `min_sat_thickness` floor keeps `T>0`. Bare-rock / soil-off cells fall back to the prescribed
  depth. Prescribed stays the default so the feature runs standalone.
- **Weathering supply `Ψ`** — three tiers, escalating cost (see §3a for the rate law + YAML):
  - **Proxy (default, shipped)** — climate/temperature stand-in; no new inputs, non-conservative.
  - **Level A (opt-in, this design)** — explicit chemical-weathering *rate* `W(R, T, Lsoil,
    lithology)` driven by the groundwater recharge `R` (± the `prodsoil` congruency shortcut when
    `soilSPL` is on). Supply-only, still non-conservative. Recommended pairing; cheap because all
    inputs already exist.
  - **Level B (future, separate module)** — conservative geochemistry: debit dissolved solid,
    transport solute along `q = −T∇h` (advection-reaction on the DMPlex, reusing the FV advection
    kernels), precipitate at the fringe, export via baseflow, with `Σ dissolved − precipitated −
    exported ≈ 0` guards. Comparable in scope to dual-lithology; needs its own design doc. The
    enabling pieces (groundwater flux, FV advection, per-class strata bookkeeping) already exist.
- **Armor of diffusion** — off by default (SPL K only); enable `armor_diffusion` if crusts
  should also resist hillslope creep.
- **Do we need transient `head` at all, or steady each step?** Carried as state with
  backward-Euler (robust in both `τ_gw` regimes). If validation shows the equilibrium limit
  everywhere, a pure steady solve (`a=0`) is a trivial simplification.
