# DESIGN: soil / near-surface regolith representation

A review of goSPL's current soil model, the inconsistencies it carries, and the
options for improving it **before** the water-table + duricrust capability
(`DESIGN_WATERTABLE_DURICRUST.md`) is finalized — because the duricrust weathering
supply, the `aquifer_base` coupling, and where induration forms all depend on how
"soil" / near-surface regolith is represented.

Status: **decision note** (no implementation yet). Companion to the `DESIGN_*.md`
series. Honors the `AGENTS.md` invariants.

---

## 1. Current model, as built (corrected)

Soil lives in the `_soilSPL` mixin (`gospl/eroder/soilSPL.py`), active when
`soil:` is configured (`self.cptSoil`). This section is the **verified** behaviour
(an earlier reading of it in the duricrust discussion was wrong on two points —
both corrected here).

### State
| Field | Meaning |
|---|---|
| `self.Lsoil` / `self.Gsoil` | soil thickness (local / global PETSc Vec) |
| `self.lHbed` / `self.gHbed` | bedrock elevation `= h − Lsoil` (recomputed in `erodepSPLsoil`) |
| `self.prodSoil` | production-rate coefficient (rainfall-scaled if `tempFile`, per Norton et al. 2013) |
| `self.Hs` | soil-production e-folding depth |
| `self.h_star` | soil-shielding depth in the *erosion* split |
| `self.Ksoil` | soil erodibility coefficient |
| `self.soil_transition` | `−ln(Sperc)·Hs` — used as a **hard maximum** soil thickness |
| `BEDROCK_EXPOSED` (`1e-1` m) | soil below this ⇒ treated as 0 (bedrock exposed) |

### Production (subaerial only)
Heimsath-type exponential, in `_form_residual_soil`:
```python
hSoil += dt · prodSoil · exp(−hSoil / Hs)          # × rainVal if tempFile set
```
Declines with thickness; **capped** at `soil_transition` on write-back.

### What becomes soil — essentially *all* deposition
`updateSoilThickness()` (`soilSPL.py:414`) does `Lsoil += Δh` (the just-applied
elevation change), floored at 0, capped at `soil_transition`. It is called from
**three** sites, so soil absorbs every depositional source:

| Source | → `Lsoil`? | Path |
|---|---|---|
| Subaerial weathering | yes | production term in `_form_residual_soil` |
| Fluvial transport-limited deposition | yes | `(h − hOldArray)` term in the soil solve |
| **Lake / pit infill** | yes | `sedplex._updateSinks:1049` → `updateSoilThickness` |
| **Marine deposition** | yes, then zeroed | `seaplex:829` → `updateSoilThickness`; then `_solveSoil` re-applies `nHsoil[seaID]=0` next step |
| Soil creep (diffusion) | yes | `soilSPL:686` (`diffuseSoil`) |

So today **`Lsoil` is a lumped "soft erodible cover"** that receives weathering
*and* all transported sediment, clamped at `soil_transition`.

### Erodibility split
```python
res += Kbr    · exp(−hSoil/h_star) · Sⁿ    # bedrock-dominated when soil thin
res += K_soil · (1 − exp(−hSoil/h_star)) · Sⁿ  # soil-dominated when soil thick
```
`Kbr = K · surfK · … `, where `surfK = _surfaceK() · _surfaceLithoK()`.
`_surfaceK()` (`stratplex.py:483`) returns the **top stratigraphic layer's `stratK`
multiplier** (1.0 when `stratNb==0`). So deposited sediment's erodibility already
enters `Kbr` **via the stratigraphy** — but only when `stratNb>0`. Both `Kbr` and
`K_soil` are zeroed at `seaID` (fluvial erosion is off underwater; marine erosion/
transport is handled by `seaplex`).

---

## 2. Problems / inconsistencies

1. **Conflation.** `Lsoil` mixes *pedogenic soil* (in-situ weathering product) with
   *transported sediment* (fluvial/lake/marine). A 50 m marine clinoform or a deep
   lake fill is not "soil"; the `soil_transition` clip hides this by saturating the
   contribution, but the concept is muddled.
2. **Submarine add-then-zero.** `seaplex` adds marine deposition to `Lsoil`, then the
   next fluvial `_solveSoil` zeroes `nHsoil[seaID]`. The submarine soil value is
   therefore transient and inconsistent (largely inconsequential today because
   `Kbr`/`K_soil` are zeroed at `seaID`, but it corrupts any downstream reader of
   soil thickness — e.g. a duricrust/regolith model).
3. **Hard cap `soil_transition`.** Defined as `−ln(Sperc)·Hs` (`Sperc = soil.bedrockConv`,
   default `1e-4` ⇒ ≈ `9.2·Hs`; `100 m` when `Sperc=0`) — the depth at which production
   has decayed to fraction `Sperc` of its surface rate. It is used **only** as a max-soil
   clip in two write-backs (`_solveSoil:363`, `updateSoilThickness:424`) — *not* in the
   erosion split (`h_star`) or bedrock-exposure (`BEDROCK_EXPOSED`). It is **not required**
   numerically (post-SNES clip, no effect on convergence); it is only **load-bearing today
   as a band-aid** — because all deposition is dumped into `Lsoil`, without it a depocenter
   would register tens of metres of "soil." Under Option 2 (soil = weathering-only) that
   justification disappears. The *concept* — a maximum weathering / mobile-regolith depth —
   is sound and is really a crude version of the Option-2.5 weathering front; it should be
   **recast as a smooth maximum-weathering-depth** on the regolith, not a discontinuous clip.
4. **`lHbed` in depocenters.** `lHbed = h − Lsoil` with `Lsoil` capped means a thick
   fill pushes `lHbed` *up* with the surface — the "bedrock" rises through a porous
   sediment pile. This is the specific reason `aquifer_base = lHbed` (duricrust design
   §8) is wrong in depositional basins.
5. **Deposited-sediment softness depends on stratigraphy.** Only `stratNb>0` gives
   deposits a soft `stratK`; with stratigraphy off, lake/marine deposits erode at raw
   bedrock `K`.
6. **Soil vs stratigraphy overlap.** Both reservoirs track "deposited soft material,"
   with no clear division of responsibility.
7. **Hillslope coupling unaudited.** Hillslope creep *is* soil transport, yet
   `diffuseSoil` (soil creep, `soilSPL`) and `getHillslope` (elevation diffusion,
   `hillslope.py`) coexist without a verified joint soil-conservation contract.

---

## 3. The three framing questions

- **No soil in the ocean?** Correct that **no *pedogenic* soil** forms underwater
  (subaerial weathering is off), so zeroing *production* is right. But the seafloor is
  *soft sediment*, not hard bedrock — its erodibility should come from the **marine
  stratigraphic record (`stratK`)**, applied consistently, not from an add-then-zero
  `Lsoil`.
- **Soil in pits / lakes?** The fill is **soft unconsolidated sediment** (stratigraphy,
  soft `stratK`), not pedogenic soil, and while under water it is **subaqueous** (no
  subaerial soil production). **Pedogenic soil forms on the fill only once it is
  subaerial/emergent.** So: deposit first (subaqueous soft sediment), then soil forms on
  top when exposed.
- **Hillslope?** This is where soil coupling should be **strengthened**: creep is soil
  transport, so diffusivity should be soil-dependent, the flux soil-conserving, and
  creep should shut off / go bedrock-limited where soil is stripped.
- **Soil under ice?** No subaerial *production* under ice (frozen, insulated, no biota /
  rain infiltration), but — unlike underwater — the existing regolith is **preserved
  (frozen inert)**, not zeroed: cold-based ice can keep a buried regolith for Myr, so
  glaciation freezes the column rather than removing it. Glacial erosion / till are handled
  by the ice model (abrasion → `Eb`; till → stratigraphy). See "Under ice" below.

### The subaerial gate — when a sink is land (goSPL already provides it)
goSPL represents a sink as a depression with **`lFill > hl`** (bed `hl` below the spill/
water level `lFill` — the pit-fill code calls `lFill` "each depression's final water level";
`pitParams` stores depth `diffh = lFill − hl`). Cells with `pitIDs >= 0` and `lFill > hl`
are the **ponded (subaqueous) lake floor**; the lake-evap budget already uses exactly this
set. Below sea level (`lFill ≤ sealevel`) it is marine (`seaID`), not a continental lake.

A sink **becomes subaerial** when the accommodation `lFill − hl` is consumed or removed:
(1) **sediment fills it to the spillover** — the main path: `_moveDownstream` fills the pit
until `inV > pitVol`, then `hl → lFill`, `pitVol → 0`, and it overspills → former lake floor
is land at spill level; (2) the lake **dries** (evaporation-limited); or (3) the depression is
**removed** (rim/spill eroded, tectonics). So pit/lake deposition is **subaqueous until
filled-to-spill**, then subaerial.

**The gate for soil/duricrust is therefore already computable from existing fields**:
subaerial ⇔ **not `seaID` and not (`pitIDs >= 0` with `lFill > hl`)**. Today soil production
is gated **only** on `seaID`, so continental lakes are wrongly treated as subaerial (and
their fill dumped into `Lsoil`) — that is the gap Option 2 closes. (The water table in
`DESIGN_WATERTABLE_DURICRUST.md` generalizes this to `h ≥ z` seepage — wetlands, near-surface
table — but the first-order lake mask exists now.)

### Under ice — frozen inert (DECIDED)
Ice-covered land is neither subaerial nor subaqueous — it is a third state. Pedogenic soil
production must be suppressed there (no subaerial weathering under ice), but the pre-existing
regolith is **preserved, not zeroed** — the "**freeze inert**" rule (chosen over the
subaqueous *zero* because cold-based ice can preserve a buried regolith for Myr). Glacial
erosion and till are the ice model's job (`_glacialAbrasion` → `Eb`; `glacialTill` →
`deposeStrat`/stratigraphy), so "what erodes/deposits under ice" is already handled — the only
gap was pedogenic soil growing under ice.

Implementation (shipped in step 1): `soilSPL._iceFrozenMask` = `iceOn` and
`iceHL > ICE_COVER_MIN` (`1e-2 m`, `constants.py`), restricted to LAND (subaqueous wins on
overlap — an ice shelf over sea/lake stays a soil-free subaqueous cell). Both `Lsoil`
write-backs (`_solveSoil`, `updateSoilThickness`) hold the ice column at its prior value
(no production, no deposition-into-soil increment). Deglaciation (mask clears) resumes normal
evolution from the preserved column. Duricrust formation (a subaerial weathering process) will
reuse the same exposed-land gate — no induration under ice.

**Refinement (open):** the freeze is a blanket rule; warm-based (fast, erosive) ice actually
strips regolith while cold-based (slow) ice preserves it. The diagnostic ice model carries
basal velocity (`iceUbL`), so a future refinement could strip under fast ice and freeze under
slow ice. Deferred — the blanket freeze is the conservative first choice.

---

## 4. Options

### Option 1 — Keep the lumped "soft cover", make it consistent (minimal)
Accept `Lsoil` = "soft erodible cover" (weathering + all deposition) but fix the
inconsistencies: decide the submarine case **once** and apply it consistently (drop the
add-then-zero); optionally replace the `soil_transition` clip with a smooth saturation;
document `Lsoil` explicitly as "soft cover, not pedogenic soil."
- **Pro:** smallest change; no new state; low risk.
- **Con:** conflation remains; `aquifer_base=lHbed` depocenter problem remains; duricrust
  weathering coupling stays crude.

### Option 2 — Separate weathered regolith (soil) from deposited sediment (recommended)
Two reservoirs with distinct roles:
- **Soil / regolith `Lsoil`** = *only* weathering-produced, subaerial, on bedrock;
  depleted by erosion; transported by hillslope creep. The true pedogenic mantle.
  Deposition **no longer routes into `Lsoil`** (remove the `updateSoilThickness` calls
  from `sedplex`/`seaplex`).
- **Deposited sediment** = the **stratigraphy** (already exists; already feeds the top
  layer's `stratK` into `Kbr` via `_surfaceK()`). Fluvial/lake/marine softness comes
  from `stratK`, not from being relabelled "soil".
- **Unified surface erodibility** = a single function of what is on top: `K_soil` where a
  weathering mantle exists, else top-layer `stratK` where sediment was deposited, else
  bedrock `K`. (Mostly a re-interpretation of the existing `Kbr`/`K_soil`/`stratK`
  machinery, not new kernels.)
- **Submarine:** no pedogenic soil (production off), seafloor soft via marine `stratK` —
  consistent, no add-then-zero.
- **Pit / lake:** fill is sediment; soil grows once emergent.
- **`lHbed`** becomes the base of the *weathered mantle* on bedrock; in depocenters the
  relevant substrate/aquifer base comes from **stratigraphy**, not a fill-inflated `lHbed`.
- **Pro:** physically clean; reuses stratigraphy; fixes the `aquifer_base` depocenter case;
  gives the duricrust a well-defined regolith to form in; strengthens soil↔hillslope.
- **Con:** rework the soil bookkeeping (stop deposition→soil) + the erodibility
  unification + tests; deposited-sediment softness needs `stratNb>0` (else deposits erode
  as bedrock, as today). Not byte-identical — changes results in soil runs (guard with
  before/after comparison + updated regression baselines).

### Option 2.5 — Option 2 plus lumped profile scalars (recommended ceiling)
Option 2, augmented with a **few per-node lumped scalars** that carry *some* vertical
information without a discretized column:
- **weathering-front depth** (≈ `z − lHbed`, already implicit),
- a bulk **weathering-degree / chemical-depletion index** ∈ [0,1],
- the **duricrust horizon depth + thickness** (from the fringe model),
- **`soil_transition` recast** as a *smooth* maximum-weathering-depth on this regolith
  (§2, item 3), replacing the hard clip.

Evolved by simple rate laws (residence time, water flux, temperature) — **no reactive-
transport column**. This is the honest fidelity ceiling at 500 m–km resolution: it says
*where* the front/crust sit and *how weathered* the mantle is, which is all a coarse cell
can support, at trivial cost.

**Feature dependencies (what needs soil vs stratigraphy):**

| Capability | Needs soil (`cptSoil`)? | Needs stratigraphy (`stratNb>0`)? |
|---|---|---|
| Water table (`head`, `wtdepth`) | No (prescribed base); soil optional (`from_soil`) | No |
| Live duricrust induration + K-armoring (this-cycle relief inversion) | No (soil-independent); soil optional | No |
| Weathering-front depth, weathering-degree index (the 2.5 scalars) | **Yes** | No |
| Duricrust burial → exhumation memory (`stratDuri`, multi-cycle) | No | **Yes** |
| Deposited-sediment soft erodibility (via `stratK`) | No | **Yes** |

So the Option-2.5 profile scalars require **soil production**, **not** stratigraphy;
stratigraphy only adds the *burial/exhumation archive* and *deposit softness* (graceful
degradation — see `DESIGN_WATERTABLE_DURICRUST.md` §9).

### Option 3 — Full weathering-profile / regolith-column model (ambitious)
A vertically structured near-surface — fresh rock → saprock → saprolite → mobile soil —
with an explicit weathering front, chemical depletion, and the **duricrust as an
indurated horizon within the profile**. The water table and duricrust live natively here.

**Scale verdict — not justified at goSPL's resolution.** *Temporally* it is well-matched:
weathering fronts descend ~1–100 m/Myr and profiles mature over 10⁵–10⁷ yr, exactly the
run length, and at 100s-yr steps the front moves sub-mm–cm/step (slow, stable). But
*spatially* the resolution defeats it: a regolith column earns its keep by resolving the
vertical horizons **and their lateral variation** at scales where that variation is real
(fracture/catchment/microtopography, metres–tens of metres). At **500 m–km cells** each
node is a huge areal average, so a zoned 1-D profile per cell is **false precision**, its
mineral-kinetic parameters are **unconstrainable** at continental/My scale, and a reactive
column per node × millions of nodes × thousands of steps is a **major cost + stiff coupling**
to a code with no subsurface discretization.
- **Pro:** most physical; duricrust + water table + weathering mutually consistent (the
  Braun-style regolith LEM); the "correct home" for a duricrust.
- **Con:** justified only at **much finer resolution** (hillslope-scale, tens of metres) or
  as a dedicated regional weathering study — **not** a continental/My goSPL run; a major
  research module that subsumes both the soil rework and the duricrust design.

---

## 5. Recommendation

**Target Option 2.5** (Option 2 + lumped profile scalars) — the honest fidelity ceiling at
goSPL's 500 m–km / My scales. **Option 3 is explicitly out of scope** at these scales (see
its scale verdict: temporally apt, spatially unjustified). Staged so the disruptive part is
optional and guarded:

1. **Consistency fixes first (safe, ≈ Option 1) — DONE.** Subaerial gate + submarine
   coherence + ice freeze-inert. `soilSPL._subaqueousMask` (marine + ponded lake → soil 0)
   and `_iceFrozenMask` (ice-covered land → soil preserved, production off). §3.
2. **Separation (opt-in behind a flag, default = current) — DONE.** `soil: mode: regolith`
   (`self.regolithSoil`, default `lumped` → byte-identical). In regolith mode `Lsoil` is the
   **weathering-produced regolith only**; deposited sediment stays in the stratigraphy and
   carries its own soft erodibility there. Mechanism:
   - **Deposition is not routed into `Lsoil`.** `updateSoilThickness(deposition=True/False)`:
     the lake/pit (`sedplex`) and marine (`seaplex`) callers use `deposition=True` → the
     increment is **skipped in regolith mode** (added in lumped); **the subaqueous/ice gates
     still run in both modes** (so a cell newly ponded by this step's deposition is re-zeroed
     consistently). Soil **creep** (`diffuseSoil`) uses `deposition=False` → always applied
     (creep transports the regolith itself, both modes).
   - **Fluvial transport-limited deposition growth** is removed at the `_solveSoil` write-back
     (`nHsoil −= max(0, Δh)`) — **post-solve, so the SNES residual and its smoothness are
     untouched**; erosion still strips soil.
   - **Fresh deposits erode like soil (option C).** A freshly deposited layer is given a soft
     `stratK = Ksoil/K` (`stratplex.deposeStrat`, regolith mode) so the SPL bedrock term
     `Kbr·stratK = Ksoil` — reusing the already-defined `soilK` (regolith mode ⇒ `cptSoil`,
     so `Ksoil` exists; **no new parameter**). Lumped mode keeps `stratK = 1.0` (there the
     deposit becomes soil and gets `Ksoil` via `updateSoilThickness`). Unified erodibility:
     bare bedrock `K`; weathering regolith `Ksoil`; fresh deposit `Ksoil`.

   **Stratigraphy is required — and it is triggered by `time: strat:`, NOT a `strata:` block.**
   `stratNb > 0` ⇔ a stratal time step is set (`inputparser`: `stratNb = (tEnd−tStart)/strat + 1`);
   the `strata:` block is only for *initial* layers / dual-lithology / bedrock sentinel. Regolith
   mode needs `stratNb > 0` so the excluded deposits are recorded with their soft `stratK`;
   without it `_surfaceK`=1.0 and a fresh deposit erodes at raw bedrock `K`. Handled with a
   **rank-0 warning** (`soilSPL.__init__`) pointing at `time: strat:`, not a hard failure
   (erosional / low-deposition runs can still use regolith mode). Guard: `test_soil_mode_regolith`
   on the **soil+stratigraphy** fixture `minimal_soil_strata.yml` (stratigraphy via `time: strat: 10`)
   — one step from an identical state gives **identical elevations** (bookkeeping-only), the
   subaerial gate holds, and fresh deposits carry `stratK = Ksoil/K` in regolith / `1.0` in lumped.
   Docs: `soil: mode:` added to `user_guide/inputfile.rst`.
3. **Lumped profile scalars (Option 2.5):** weathering-front depth, weathering-degree index,
   duricrust horizon; recast `soil_transition` as a **smooth** max-weathering-depth (§2, item 3).
4. **Hillslope coupling audit:** verify/enforce joint soil conservation between
   `diffuseSoil` and `getHillslope`.

Option 3 stays on the horizon only if a much-higher-resolution or dedicated regional
weathering study ever motivates it.

---

## 6. Interaction with the water-table + duricrust design

Option 2 tightens `DESIGN_WATERTABLE_DURICRUST.md` in three places:

- **Weathering supply `Ψ` / Level-A rate (§3a):** draws on the *weathering regolith*
  reservoir and `prodSoil` cleanly, without contamination by marine/lake sediment mislabelled
  as soil.
- **`aquifer_base = lHbed` (§8):** well-defined as the base of the weathered mantle;
  in depositional basins the aquifer base comes from **stratigraphy**, resolving the
  documented depocenter caveat.
- **Duricrust host (§8):** the crust indurates the *regolith* (or in-situ rock), a
  physically coherent medium — not a lump of "soil" that happens to include marine mud.

Under the current (unmodified) model, the duricrust design still works but must keep the
caveats already written in §8/§15 (soil is a lumped soft cover; `aquifer_base=lHbed`
degrades in depocenters).

---

## 7. AGENTS.md / migration concerns

- **Byte-identical default.** Any change ships behind a `soil: mode:` switch (or equivalent)
  whose default reproduces current output; the separation mode is the opt-in. Guard with a
  before/after regression on a soil run.
- **`destroy_DMPlex`.** No new persistent Vecs in Options 1–2 (they re-interpret existing
  state); Option 3 would add regolith-horizon state to the destroy list.
- **MPI contract.** `updateSoilThickness` and the erodibility assembly are per-node/owned-rows
  operations; removing/relocating the `updateSoilThickness` calls must not move a collective
  under a rank-local guard.
- **Regression baselines.** Option 2 changes soil-run results → update the soil regression
  fixtures and document the physical rationale (as dual-lithology/ice did).

---

## 8. Open questions

- **Emergent-fill weathering:** how fast should soil (re)establish on a newly subaerial
  lake bed / marine terrace? Straight `prodSoil` from zero, or a head-start?
- **Deposited-sediment softness without stratigraphy:** if `stratNb==0`, deposits erode as
  bedrock. Acceptable, or should a minimal "fresh-deposit soft cover" exist independent of
  stratigraphy?
- **Smooth vs hard cap:** replace the `soil_transition` clip with a smooth saturation, or
  keep it (post-solve, so harmless to convergence)?
- **Hillslope/soil conservation:** is soil currently conserved across `diffuseSoil` +
  `getHillslope`, or is there double-counting to fix?
- **Scope commitment:** Option 2 now, or hold for Option 3 (full regolith profile) if the
  duricrust is going to motivate that investment anyway?
