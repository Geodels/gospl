# DESIGN: Dual-lithology (coarse/fine) sediment option

Status: **design / not yet implemented**. Target branch: `dev`.
Last updated: 2026-06-12.

This document is the implementation plan for adding an **opt-in dual-lithology**
(coarse + fine) sediment capability to goSPL, active only when stratigraphy is
turned on. When the option is **off (default), behaviour is bitwise-identical to
the current single-fraction code** — every dual-lithology path is gated behind
the flag.

Read this alongside `AGENTS.md` (invariants), `docs/HOW_TO_ADD_FORCING.md`, and
`docs/HOW_TO_ADD_OUTPUT.md`. It supersedes the "single value per layer" caveats
in `docs/tech_guide/strat.rst` once shipped.

---

## 1. Scope & decisions (locked)

| Decision | Choice |
|---|---|
| Number of lithologies | **Two**: coarse (sand) + fine (silt/clay) |
| Composition source | **Per stratigraphic layer** — each layer carries its own coarse/fine split; eroded material inherits the composition of the layer(s) it came from |
| Transport | **Separate** — fines travel farther (distal basin / depocenter); coarse stays proximal |
| Lake/depression behaviour | Coarse builds **inlet deltas**, fines settle in the **depocenter**; fine-enriched overspill may bypass downstream while coarse is trapped |
| Fortran kernel | **Call `strataonesed` twice** (once for the total/coarse pile, once for the fine pile) — no Fortran change. NOT `stratathreesed`: it treats its extra fields as 0–1 fractions and renormalises (`functions.F90:1968`), which conflicts with the bulk-thickness representation used here. |

The dormant `stratathreesed` 3-fraction kernel exists but is never called and
its fraction semantics don't match this design's bulk-thickness `stratHf`;
Phase 5 reuses `strataonesed` (already bound, `stratplex.py:13`) a second time
for the fine pile instead. This design wires up the
2-fraction subset of it.

---

## 2. State variables

Keep `self.stratH` = **total** layer thickness (so all existing compaction /
elevation sums in `getCompaction`/`_depthPorosity` keep working) and add parallel
fraction arrays, all shape `(lpoints, stratNb)`:

| New array | Holds | Mirror of |
|---|---|---|
| `self.stratHf` | fine-fraction layer thickness (coarse = `stratH − stratHf`) | `stratH` |
| `self.phiF` | fine porosity per layer (`phiS` is now coarse porosity) | `phiS` |

Erosion-rate outputs become a pair: `self.thCoarse` **and** `self.thFine`
(both m/yr, solid, uncompacted-equiv). Flag: `self.stratLith` (bool, default
`False`).

npstrata input gains optional `strataHf` + `phiF`; absent → all-coarse, so old
input files still load unchanged.

---

## 3. YAML opt-in

```yaml
strata:
  dual: True                          # default False -> current single-fraction path
  coarse: {phi0: 0.49, z0: 3700.}     # coarse porosity-depth (today's phi0s/z0s = coarse)
  fine:   {phi0: 0.63, z0: 1960., k_factor: 1.0}   # fines: porosity-depth + erodibility ratio
  bedrock_coarse_frac: 0.5            # composition of the layer-0 bedrock sentinel
  # transport / deposition contrast (the "fines travel farther" knobs):
  fine_efficiency: <factor>           # marine: per-fraction Gmar (fine deposits less readily)
  pitInletBias: {coarse: 0.5, fine: 0.0}   # lacustrine: coarse delta vs fine depocenter
  fine_diff_factor: 1.0               # fine diffusivity relative to coarse (1.0 = no contrast)
```

**Diffusivity contrast is a multiplier, not a replacement.** goSPL already has
the base diffusion coefficients `hillslopeKa` (→ `self.Cda`, subaerial linear),
`hillslopeKm` (→ `self.Cdm`, subaqueous linear), and `nonlinKm` (→ `self.nlK`,
non-linear) — all grain-size-agnostic today. `fine_diff_factor` layers the
lithology contrast on top of whichever base coefficient a cell uses:

> `Cd_eff = Cd_base * (fc + ff * fine_diff_factor)`,  with `Cd_base ∈ {Cda, Cdm, nlK}`.

So the existing land/sea/linear/non-linear structure is untouched;
`fine_diff_factor == 1.0` (default) or `dual: False` ⇒ `Cd_eff == Cd_base`
(byte-identical). This mirrors `fine.k_factor` for erodibility — one
`_surfaceComposition()` feeds both blends. (Earlier drafts used absolute
`Dc`/`Df`; dropped because they would duplicate and risk contradicting
`hillslopeKa`/`hillslopeKm`/`nonlinKm`.)

Parsed by a new `_extraStrata` method in `tools/inputparser.py`, following the
**mandatory `_extra*` continuation convention** (add it to the `_extra*` table in
`AGENTS.md`). Sets `self.stratLith`, the two porosity curves, `fine_k_factor`,
`fine_efficiency`, per-fraction inlet bias, and `fine_diff_factor`. All keys
default such that `dual: False` reproduces current behaviour.

---

## 4. The two shared hooks (the core simplification)

Lithology enters the physics through **one new helper plus the existing erode/
depose path** — *not* through per-flavour math:

- **`_surfaceComposition()`** — new helper in `stratplex.py`, modeled exactly on
  the existing `_surfaceK()` (`stratplex.py:117`). Returns the exposed surface
  layer's `(fc, ff)` per node. Returns all-coarse when `stratLith` is off or the
  column is empty.

This single helper feeds **all three** consumers:
1. **Erodibility blend** — effective `K = fc·K_coarse + ff·K_fine`, inherited by
   all three SPL flavours because they already scale `self.K` through
   `_surfaceK`.
2. **Eroded-mass split** — `erodeStrat` splits the eroded solid by the consumed
   layers' composition (Section 6).
3. **Diffusivity blend** — hillslope `Cd_eff = Cd_base*(fc + ff*fine_diff_factor)` (Section 7).

---

## 5. Erosion across the SPL flavours

The SPL solvers produce a **total** erosion thickness rate (`Eb`, m/yr); they do
**not** compute the coarse/fine split themselves. Lithology enters only via the
two shared hooks above.

| Flavour | Change for dual lithology |
|---|---|
| **SPL, detachment-limited** (`fDepa==0`, `SPL._getEroDepRate`) | None in the solve. K-blend + `erodeStrat` split. Deposition entirely in `sedplex` (Phase 3). Trivial. |
| **nlSPL, detachment-limited** (`nlSPL._solveNL`) | Same — the nonlinear `n` acts on total stream power. Trivial. |
| **SPL/nlSPL transport-limited** (`SPL._coupledEDSystem` `fDepa!=0`; `nlSPL._solveNL_ed`) | The one complication — deposition is baked into the erosion solve via a single `G`. **Resolution (A, default):** run SPL erosion-only in dual mode and defer per-fraction deposition (`Gc`/`Gf`) to the `sedplex` routing of Phase 3. **(B, opt-in later):** run the coupled solve twice with `Gc` then `Gf` (doubles the AD-HOC fieldsplit / cached `_snes_ed` cost). |
| **soilSPL** (`soilSPL._solveSoil`) | Most lithology-intrinsic: the soil layer carries its own composition, and soil production weathers bedrock into **fine-enriched regolith** (`prodSoil`). Soil erosion splits by soil composition, exposed-bedrock erosion by `bedrock_coarse_frac`. |

---

## 6. Conservation: what happens when transport demand exceeds the top layer (CRITICAL)

`D_eff` (and `K`) set the **magnitude** of total transport, never the
**composition** of what is removed. Composition is always reconciled against the
strata record by `erodeStrat`, integrating **layer by layer downward** (the
existing `cumThick` consume loop, `stratplex.py:195-219`):

1. The solve yields a **total** erosion depth at each node.
2. `erodeStrat` consumes downward: takes all the fine in the top layer, then the
   remaining demanded volume from the next layer at *that* layer's composition,
   and so on. **A fraction can never be over-drawn** — each layer caps its own
   contribution.
3. The **bedrock sentinel** (`BEDROCK_SENTINEL`, `stratplex.py:99`, composition
   `bedrock_coarse_frac`) is the ultimate backstop, so total volume is always
   realizable — never a true "ran out," just erosion into bedrock.
4. Deposited composition ≡ eroded composition → **per-fraction mass conserved
   exactly**; no negative fractions.

**Residual = a one-timestep explicit lag**, not a conservation error: `D_eff`
used the *surface* composition while the removed mass is *layer-integrated*. It is
bounded (diffusion is CFL/`tsStep`-limited), self-corrects next step (surface
composition then reflects the freshly-exposed layer), and composition is frozen
during the implicit elevation solve (operator-split; avoids a doubly-nonlinear
SNES/TS). A guard + a "diffuse-through-a-compositional-boundary" regression test
will flag if the lag ever matters; `D_eff` substepping is the escape hatch, not
built speculatively.

Mirror the existing negative-thickness clamp (`stratplex.py:226-228`) for
`stratHf`/`phiF`.

---

## 7. Diffusion (hillslope + subaqueous)

Today hillslope is **single-rate** for both grain sizes (`hillslope.py:64`).
goSPL already has the base diffusivities `Cda`/`Cdm` (linear land/marine, YAML
`hillslopeKa`/`hillslopeKm`) and `nlK` (non-linear, `nonlinKm`). Add the
grain-size contrast as a **dimensionless multiplier** `fine_diff_factor` (NOT a
new absolute coefficient — see §3) and evaluate via **composition-weighted
effective diffusivity in a single solve**:

- `Cd_eff = Cd_base * (fc + ff*fine_diff_factor)`, with `fc/ff` from
  `_surfaceComposition()` and `Cd_base ∈ {Cda, Cdm, nlK}` (the cell's existing
  land/marine/non-linear coefficient). The diffused flux carries the **donor
  cell's composition** into `erodeStrat`/`deposeStrat`.
- **Drop-in** to the existing solvers — the FV scheme already supports node-
  varying diffusivity (the nonlinear hillslope varies D by slope via
  `Dlimit`/`dexp`); composition is just another multiplicative factor. **No new
  cached solver** (so no `destroy_DMPlex` churn) and **no two-field nonlinear
  coupling**.
- Neutral by construction: `fine_diff_factor == 1.0` or `dual: False` ⇒
  `Cd_eff == Cd_base`, byte-identical to today. Mirrors `fine.k_factor`.
- A genuine two-field separate-D solve is more accurate but much heavier; offer
  as a future opt-in only.

**Subaqueous diffusion is where the contrast does the real work.** The marine +
lake diffusion path (cached `_ts_marine`, `_diffuseImplicit` in `hillslope.py`)
with `fine_diff_factor > 1` (on the `Cdm`/`nlK` base) is *the mechanism* by which
fines spread to the depocenter while coarse stays proximal — i.e. it implements
"fines travel farther" consistently with the lacustrine inlet-bias and marine
`Gmar`:

| Domain | "Fines travel farther" realized by |
|---|---|
| Subaerial slopes | `Cd_eff` weighting on `Cda` (`fine_diff_factor > 1`) |
| Lakes/depressions | per-fraction inlet bias + subaqueous `Cd_eff` |
| Marine | per-fraction `Gmar` + subaqueous `Cd_eff` (`Cdm`/`nlK` base) |

---

## 8. Lake / depression deposition (`sedplex`)

Continental closed-depression deposition is a third depositional environment,
separate from river transport and the marine (`seaplex`) clinoforms. The existing
machinery already has the two end-members dual lithology needs:

- **`_bottomUpDelta`** (`sedplex.py:230`) — deepest-first fill → the **fine /
  depocenter** path.
- **`_diffuseLargePit`** (`sedplex.py:385`) + inlet bias (`nl_pit_inlet_bias`,
  `inputparser.py:494`) → the **coarse / inlet-delta** path.
- **`_moveDownstream`** (`sedplex.py:84`) — pit fill + overspill cascade.

Dual-lithology deposition is therefore **not new geometry** — it splits the per-
pit deposited volume `depo → depoC, depoF` and routes each through the existing
path, with `pitInletBias` per fraction (coarse high, fine ~0). New work:

1. **Per-fraction pit budget** — carry `vSedC`/`vSedF` through `_moveDownstream`;
   both fractions draw down the shared `self.pitVol` accommodation while placing
   coarse-at-inlet and fine-in-depocenter, conserving total = `pitVol`.
2. **Per-fraction deposit shape** — coarse → `_diffuseLargePit`, fine →
   `_bottomUpDelta`.
3. **Fine-enriched overspill** — the excess `eV` overspilling a filled lake
   (`sedplex.py:118-126`) carries composition; coarse is largely trapped, so
   overspill is fine-enriched (and may bypass downstream).
4. **Fresh porosity** — lacustrine coarse → `phi0c`, fine → `phi0f` (feeds the
   per-fraction compaction).
5. **Per-fraction mass conservation across the cascade** — the spillover cascade
   (`step>0`) and the lake-evap budget (`flowplex.py:364-375`, applied at
   `step==0` only) must conserve each fraction independently. This is the
   lacustrine analog of the documented marine terminal-sink leak (`AGENTS.md` >
   Fixed > Marine sediment leak).

---

## 9. Phased implementation plan

Land **Phases 0–2 + 4–5 first as a co-transport increment** (composition tracked,
fines compact differently, routed *with* coarse), get it green, **then** do the
separate-transport Phase 3 as a focused second PR with a full regression run. This
de-risks the routing surgery.

| Phase | Scope | Files | Risk |
|---|---|---|---|
| **0** | Flag + `_extraStrata` parser + backward-compat guard. *Acceptance: dual-off run bitwise-identical.* | `inputparser.py`, `AGENTS.md` | low |
| **1** | Allocate `stratHf`/`phiF` (both branches of `readStratLayers`); write/restore in I/O + XDMF entries | `stratplex.py:45-103`, `outmesh.py` | low |
| **2** | `_surfaceComposition()` helper; K-blend hook; eroded-mass split in `erodeStrat` | `stratplex.py` | medium |
| **4** | Per-fraction deposition (`deposeStrat`) + per-fraction compaction (`_depthPorosity`/`getCompaction`) | `stratplex.py:146-345` | medium |
| **5** | Stratal advection: second `strataonesed` call for the fine pile (`stratHf`, `phiF`) | `stratplex.py` (`stratalRecord`) | medium |
| **3a** | River routing: split `vSed → vSedC/vSedF` in `updateSedLoad`/`_getSedFlux` | `sedplex.py:46-83` | **high** |
| **3b** | Lacustrine/pit deposition (Section 8) | `sedplex.py` | **high** |
| **3c** | Marine routing: per-fraction `Gmar` in `_distOcean`/`_depMarineSystem` | `seaplex.py` | **high** |
| **6** | Tests (Section 10) | `tests/`, `benchmarks/` | — |

Phase 3 touches the high-risk modules (`flowplex`, `sedplex`, `seaplex`,
`unstructuredmesh`) — it is ~70% of the effort/risk. Phases 0–2, 4–5 are mostly
local to `stratplex`/`inputparser`/`outmesh`.

---

## 10. AGENTS.md contracts to honor

- **`destroy_DMPlex` registration** (`unstructuredmesh.py:738-799`): any new
  persistent Vec (`vSedF`, etc.) must be added to the destroy list or it leaks.
- **`i`-suffix flow arrays**: SPL/sed kernels must use `rcvIDi`/`fMati`/etc.
  (pre-fill snapshot), not the live `rcvID`/`fMat`.
- **`Eb`/`EbLocal` thickness-rate sign convention** (positive deposition): both
  fractions follow it; mirror the `sedplex._getSedFlux` negate-for-upstream quirk
  per fraction.
- **Scratch Vec contract**: document any `self.tmp`/`tmpL`/`tmp1`/`h`/`hl`/`dh`
  reuse in new method docstrings; leave them in a defined state on exit.
- **KSP/SNES/TS lifecycle**: prefer reusing existing cached solvers with node-
  varying coefficients over adding new ones; any genuinely new cached solver →
  add to `destroy_DMPlex`.
- **Forcing/DataFrame named access**, **constants in `constants.py`**, **version
  bump lockstep** — as per the standing rules.

---

## 11. Tests (Phase 6)

1. **Dual-off reproduces single-fraction** bitwise (regression guard for the flag).
2. **Per-fraction mass conservation** on the `minimal` global-sphere fixture
   (closed domain; analog of `test_mass_conservation`) — assert eroded == deposited
   per fraction to the documented floor tolerance.
3. **Fines deposit more distally than coarse** (river → lake/marine partitioning).
4. **Diffuse-through-a-compositional-boundary** fixture — bounds the explicit
   `D_eff` lag (Section 6).
5. `model.destroy()` in `try/finally`; benchmarks marked `@slow @benchmark` with
   `pytest.importorskip` for scipy/matplotlib.

---

## 12. Open decisions

- Whether to ship a dedicated `stratatwosed` Fortran kernel (one pass instead of
  two `strataonesed` calls — interpolates the redundant elevation only once).
  Implemented Phase 5 with two `strataonesed` calls (no f2py change); a fused
  kernel is a later micro-optimization.
- Whether transport-limited SPL gets per-fraction coupled solves (option B,
  Section 5) or stays deferred-to-sedplex (option A). Start with A.
- Marine `fine_efficiency` vs lacustrine `pitInletBias`: kept as separate tunables
  (different energy regimes) but documented as a pair.
