# Design note: carbonate platforms and reefs in goSPL (draft for discussion)

Branch `feat/carbonates` (from `feat/gravity-currents`). Status: design, no code.

**Decisions (maintainer, 2026-11-24):** target both regional platforms / reef terraces
and global carbonate burial; carbonate is a true lithology (third, with its own
properties), plus hemipelagic rain; simple response tables (no fuzzy rule system); no
cellular automaton. Sections 2.4 and 4 below are kept as the record of the evaluation.
Sources read (2026-11-24, in `~/Documents/carb`):

- **pyBadlands** (Salles, Ding & Brocard 2018, PLoS ONE): reef growth from fuzzy logic on
  three controls (depth, wave energy, sedimentation rate), one or more coral assemblages,
  defuzzified by centroid; mixed siliciclastic-carbonate example on the GBR (14 kyr).
- **Pastier et al. 2019** (G-cubed): 2-D profile model of Quaternary reef terrace
  sequences. Effective growth `G = gamma(depth) * zeta(distance to open water) * G*`,
  with `gamma = (1 + cos(pi h / zmax))/2` and a linear cut above the wave depth `zmin`,
  `zeta = (1 + tanh((x_ow - x)/delta))/2`; wave erosion after Anderson et al. (1999),
  `dEr/ds = K Er exp(-h/z0)`, dissipated shoreward, residual carves a notch and cliff;
  eroded material fills lagoons to `z0`, the rest is laid on the fore reef at a rest angle.
  Main finding: vertical land motion and foundation slope control the architecture, the
  growth rate only matters secondarily, and erosion controls the barriers.
- **CarboKitten.jl** (Hidding et al. 2026, GMD): Bosscher & Schlager (1992) production per
  factory, `P_f(w) = g_m tanh(I0 exp(-k w) / I_k)` (euphotic / oligophotic / aphotic,
  optional pelagic as the integral over the water column); Burgess (2013) cellular
  automaton (5x5 neighbourhood, viability 4-10, activation 6-10) to choose the active
  factory per cell; transport of an active layer `q_f = -C_f d_f grad(eta) + C_f v_f(w)`
  fed by production and disintegration, then lithification with a half-life. They note
  that this transport behaves as a diffusion of the topography.

## 1. What goSPL already has

Most of what a carbonate system needs downstream of production is already in goSPL:

| Need | goSPL piece |
|---|---|
| Water depth, sea-level curve, subsidence, flexure | `sealevel` forcing, tectonics, flexure |
| Shelf redistribution of loose sediment | marine routing `_distOcean` + donor marine diffusion |
| Slope and basin redistribution, slope failure | gravity currents (`gravity:`, this branch) |
| Grain-size classes | dual lithology (coarse = carbonate sand, fine = carbonate mud) |
| Identity of the sediment through erosion, transport, strata | provenance classes (`stratP`, machine-exact) |
| Harder reef framework | per-layer erodibility `stratK` |
| Subaerial dissolution of exposed carbonate | Level-B geochemistry (carbonate solute, weatherability per class) |
| Compaction | per-lithology porosity curves |

The missing parts are: (1) in-situ production, (2) the split between framework (stays in
place) and detritus (enters the routing), (3) the controls that modulate production, and
for terraces (4) coastal wave erosion.

## 2. Proposed approach

### 2.1 Production (phase C1)

Per factory `f` (user list, default euphotic / oligophotic / aphotic, optional pelagic),
on each marine node with water depth `w = sealevel - h > 0`:

    P_f = g_f * L_f(w) * prod_j R_fj(x_j) * M_f(node, t)

- `L_f(w) = tanh(I0 exp(-k_f w) / I_k,f)`: Bosscher & Schlager light response (the
  CarboKitten defaults: `g` 500/400/100 m/Myr, `k` 0.8/0.1/0.005 m^-1, `I_k` 60, `I0` 400).
  Pastier's cosine `gamma` can be offered as an alternative shape.
- `R_fj`: response curves to other controls, each a piecewise-linear table of
  (value, factor in [0, 1]) given in the YAML. This is the "fuzzy-lite" form of the
  Badlands approach: a product of membership functions, no rule base or defuzzification,
  which is easier to calibrate and to explain. Candidate controls available in goSPL each
  step: siliciclastic deposition rate and its fine fraction (turbidity, burial),
  river discharge reaching the coast (fresh water, nutrients), seafloor slope, and an
  exposure / open-water index (2.3).
- `M_f`: optional habitat map in time (a forcing like `rainMesh`, e.g. from SST or
  latitude on global runs), following `docs/HOW_TO_ADD_FORCING.md`.

Accommodation: production per step is integrated per node as the ODE
`dh/dt = P(sealevel - h)` over `dt` (sub-steps, local, no communication), so a node
cannot grow above `sealevel - w_min` and keep-up, catch-up and give-up emerge from the
law. goSPL steps are 1-10 kyr, reef rates 1-10 mm/yr, so 1-100 m per step: the cap and
the sub-stepping matter. All of this is pointwise, so it is partition-exact by
construction and costs nothing next to the solves.

### 2.2 Framework versus detritus (phase C1)

Each factory has a framework fraction `phi_f` (reef boundstone: high for euphotic
corals, ~0 for pelagic and most oligophotic producers):

- `phi_f * P_f * dt` stays in place: raises `h`, is deposited in the top stratigraphic
  layer with the factory's provenance class and a reduced erodibility `stratK`.
- `(1 - phi_f) * P_f * dt` is loose sediment, split into coarse/fine by a per-factory
  sand fraction, added as a distributed source to the marine routing (`_distOcean`
  accepts a per-node source) with its provenance class, then moved by the marine
  diffusion and the gravity currents like river sediment.

This is the main design choice, and the reason to do it in goSPL rather than port
CarboKitten: carbonate detritus and siliciclastics share one transport, so mixed systems
(reefs smothered by deltas, carbonate aprons, allodapic turbidites) come out of the same
budget. The transport law of CarboKitten (an active layer that behaves diffusively) is
already covered by the marine diffusion; no new transport solver is needed.

Order in the step: production before the marine deposition, so the detritus joins the
same routing pass; the siliciclastic controls (`R` on deposition rate, turbidity) then
use the previous step's values (one-step lag, acceptable at kyr steps, avoids a second
marine pass).

### 2.3 Exposure / open-water index (phase C2)

Pastier's `zeta` (distance to open water) and Badlands' wave energy both favour the
platform margin over the lagoon. Options, cheapest first:
- bathymetric openness: `h - smooth(h)` with the existing marine smoother, positive on
  rims and highs;
- fetch-weighted depth upwind of a prevailing wave direction (a sweep, needs care in
  parallel);
- distance to deep water by a parallel eikonal or diffusion solve.
Start with the first; it needs no new solver.

### 2.4 Cellular automaton (deferred, optional)

The Burgess CA gives ecological patchiness at the scale of tens to hundreds of metres.
On goSPL meshes (hundreds of metres to tens of kilometres per cell) the patches would be
set by the mesh, not by ecology, and a CA needs per-iteration halo exchanges and a
partition-invariant random seed (by `locIDs` and step). Recommendation: leave it out of
the first phases; offer later, if needed, for fine regional meshes, possibly as
spatially correlated noise on `M_f` rather than a CA.

### 2.5 Coastal wave erosion (phase C3, needed for terraces)

Pastier shows that erosion controls the number of barriers and the terrace widths, so a
terrace study needs it. A depth-only form on the mesh: erosion rate
`E = E* K exp(-d/z0)` for `0 < d < z_wave`, plus a notch rate at the shoreline cells
(the dissipation along a profile does not translate directly to 2-D; a shoreward
attenuation could use the exposure index). Eroded material enters the marine routing.
This is useful beyond carbonates (rocky coasts, wave-cut platforms), so it would be its
own opt-in block.

### 2.6 Exposure on land (later)

Exposed carbonate: karst dissolution through the existing geochemistry (carbonate
weatherability per provenance class), and the reef `stratK` for mechanical erosion.

## 2.7 Lithology architecture (decided: a true carbonate lithology)

Dual lithology stores one extra pile per layer (`stratHf`, fine; coarse is the
remainder `stratH - stratHf`) and threads it through every process: about 160 lines in
`stratplex`, `sedplex`, `iceplex`, `hillslope`, `seaplex`, `outmesh` and the section
tools. A carbonate pile written the same way would duplicate every one of those paths.

Proposal: **generalise the fine pile to a list of lithologies** before adding carbonate.

- Coarse clastic stays the implicit remainder. The extra lithologies are an array axis:
  `stratHl[node, layer, l]`, porosity `phiL[node, layer, l]`, routed sub-fluxes
  `vSedL[l]`, eroded `thL[l]`, deposit fractions `depoFracL[l]` (the provenance tracers
  already work this way for their classes, and conserve exactly).
- Each lithology carries its properties: porosity curve (`phi0`, `z0`), erodibility and
  diffusivity factors (today `fine_k_factor`, `fine_diff_factor`), its grain-size
  behaviour (`coarse` or `fine`: pit retention order, marine fine-fraction rule, gravity
  run-out), and optionally `eroded_as` (reef framework eroded becomes carbonate sand).
- `strata: dual: True` becomes the list `[fine]` and must stay **bitwise identical**
  (the existing guards `test_dual_all_coarse_matches_single_fraction`,
  `test_dual_fine_conservation` and the dual examples, compared before/after).
- The carbonate block adds `[carb_sand (coarse-like), carb_mud (fine-like)]` and
  optionally `reef` (framework: low erodibility, eroded as `carb_sand`). Hemipelagic
  clay is the existing `fine`.

Phase C0 is this refactor alone, with no new physics, so it can be checked bitwise.

## 2.8 Pelagic and hemipelagic rain (phase C2)

Settled, not routed: a thickness rate added on marine cells each step.
- Pelagic carbonate (`carb_mud`): export production from a map in time (like
  `rainMesh`; on global runs from latitude or a productivity / SST field, the PNAS 2025
  carbonate burial regimes), times a dissolution taper with depth between the lysocline
  and the carbonate compensation depth (both may vary in time). Above the CCD it is
  deposited, below it dissolved (budgeted as dissolved, not lost).
- Hemipelagic clay (`fine`): a uniform rate or a map, optionally decreasing with distance
  from the coast.
Optionally the rain can be handed to the marine diffusion instead of deposited in
place (resuspension); first version deposits in place.

## 2.9 Revised phases

**Status (2026-11-24): C0 to C7 done (C0, C1, C2, C2b, C3, C4, C5, C6, C7, C7b)** (C3: `exposure` response control, screened-Poisson open-water index over a length in metres, not the marine smoother whose strength is in cells) (C2: `pelagic:` block, `CarbMesh.pelagicRain`; C2b: `clay_decay` and `resuspend`) (C1: `carbonate:` block, `gospl/sed/carbplex.py`, benchmarks in `benchmarks/test_carbonate.py`; the keep-up depth `w*` turned out to be the give-up threshold, an unstable equilibrium, not an attractor). C0 on `feat/carbonates` (stage A: table, piles, erosion,
deposition, compaction, blends, advection, outputs; stage B: the routed composition
at every deposit). Checks: dual lithology bitwise identical on 13 variants; a fine or
coarse fraction split into two identical lithologies reproduces the dual run to
round-off on all of them, with each lithology conserved; np=1 vs np=3. The check
harness lives in `.c0/` in the working copy (not in git: `baseline.py`, `split.py`).
Found on the way, not part of C0: the groundwater recharge and seepage gates test
`lFill > h` without `FILL_LEVEL_TOL`, so a lake cell at its own surface flips on
round-off (a base-model fix for `dev`).

| Phase | Content | Check |
|---|---|---|
| C0 | Lithology list refactor (fine pile generalised) | dual bitwise identical, full suite, `-m mpi` |
| C1 | Carbonate lithologies, benthic production (factories x response tables, per-node ODE with sea-level cap), framework vs detritus into the marine routing, outputs | keep-up depth, drowning, `dt` independence, conservation, np=1 vs np=3 |
| C2 | Pelagic carbonate and hemipelagic rain, lysocline / CCD | budget per step, global example |
| C3 | Open-water index (marine smoother) | regional platform example |
| C2b | Hemipelagic clay decreasing away from the coast (`decay` length L: `exp(-d/L)` from a screened-Poisson solve `(L^2 lap - 1) c = 0`, c = 1 at the coast, one cached KSP per step, parallel; reuses the clinoform coast distance if it is already computed); optional hand-over of the settled rain to the marine diffusion (resuspension), each part riding the lithology lockstep of `_diffuseOcean` | 1-D strip against exp(-x/L); volume and per-lithology conservation through the diffusion; np=1 vs np=3 |
| C4 | Coastal wave erosion (own opt-in block) | Pastier-type terrace example |
| later | Karst of exposed carbonate through the geochemistry | |

## 2.10 Requested extensions (2026-11-24, proposals)

**C5, temperature (done: a top-level `sst:` block, so the pelagic rain can use it without the carbonate block).** A `sst:` forcing, read in the
carbonate module (local slices, like the pelagic maps, so no new dispatcher in
`unstructuredmesh`): events `{start, uniform | map: [file, key]}` (global: SST
maps through time), or `curve: file` (two columns, time and SST in degrees C,
linearly interpolated; regional), optionally both (map shifted by the curve).
A new response control `temperature` per factory, e.g. reef builders
`x: [16, 18, 28, 32], f: [0, 1, 1, 0]`. It is also available to the pelagic
rain (`carbonate_temperature` table on the rain rate). Pointwise, cheap.

**C6, link to the geochemistry (done: `solute` response, `carbonate_supply:`; yield in the `solid_volume` convention)** (needs `groundwater:` with Level-B
geochemistry and `marine_coupling`). Two couplings, both optional:
- regional, coastal plume: the river dissolved load delivered at the coast
  (`marineSoluteInput`, per species) spread offshore with the screened-Poisson
  operator of C2b (`(I - L^2 lap) c = input`), giving a plume concentration;
  a response control `solute: <species>` lets it raise production (alkalinity,
  Ca) or cut it (nutrients, terrigenous input) near river mouths;
- global, supply-limited burial: the marine reservoir of a species
  (`marineSolute`, fed by weathering through the rivers) caps the carbonate
  buried in a step (benthic + pelagic), with a `yield` (m3 of carbonate per
  unit of solute); production is scaled down uniformly when it would exceed
  the reservoir, and burial is removed from it. This closes the weathering to
  burial loop of the carbonate system; optionally the CCD responds to the
  reservoir (deeper when supply exceeds burial). Units and the stoichiometric
  yield must be settled with the maintainer before coding.

**C7, wave energy and shadows (done: `waves:` block, `gospl/sed/waveplex.py`).** A wave energy field `E` (0 to 1) carried
from the open ocean along the wave directions, absorbed by land and dissipated
in shallow water, solved like the orographic rain shadow (same upwind
advection-relaxation operators, `getfacevelocity` + `advecupwind`, cached):
`(v . grad + kappa(h)) E = 0`, `E = 1` on the inflow edges, `kappa` large on
land (shadow behind islands and headlands) and growing as the water shallows
(friction, Pastier's shoreward dissipation). A wave climate is a list of
directions with weights, `E = sum_k w_k E_k`. A response control `wave`
(0 to 1) concentrates reefs on exposed margins and keeps lagoons and lee sides
sheltered; C4 (coastal wave erosion) then uses the same `E` as its erosional
power. Flat meshes first (the orographic operator is flat-only); on a global
mesh the directions would be a field (e.g. trade winds by latitude), later.

Proposed order: C5 (small, independent), C7 (also prepares C4), C6 (after the
units are agreed), then C4.

**C7 on global meshes (added 2026-11-24).** Bearings become a tangent field
(local east/north on the sphere), with `from: zonal` belts and `from_map`; a
fetch growth term `(d.grad + 1/L_f + kappa) E = 1/L_f` replaces the inflow
edges a sphere does not have (behind an island the waves grow again over the
fetch).

**C7b, refraction (done; option 1).** `refraction: True` in the
`waves:` block (off by default): in shallow water the local direction turns
towards the shore by Snell's law, `sin(theta)/c(h)` constant with
`c = sqrt(g h)` and theta measured from the normal to the depth contours of a
smoothed depth gradient; deep water keeps the offshore direction. The upwind
coefficients are then recomputed each step (one Fortran call per direction).
Exact test: plane beach with straight contours, direction against Snell's
angle at every depth. No caustics (smoothed by upwinding), no diffraction
(a small lateral diffusion could be added). A spectral model (directional
bins) and ray tracing were considered and rejected (cost, parallelism).

Agreed order: C7 global, C6, C7b refraction, C4.

## 3. Outputs, tests, benchmarks

- Outputs: production rate per factory, cumulative carbonate thickness, surface
  carbonate fraction; stratigraphy carries the factory classes (provenance), so the
  section tools can draw carbonate facies.
- Off-path byte identity; conservation: produced = framework + routed deposit + exported.
- Benchmarks with exact answers: (a) keep-up depth under constant subsidence `sigma`,
  `P(w*) = sigma` (Bosscher & Schlager steady state); (b) drowning when the rate of
  relative sea-level rise exceeds `g_m L(0)`; (c) the ODE integration independent of
  `dt`.
- Example: a Pastier-type uplifting / subsiding margin under a Quaternary sea-level curve
  (qualitative comparison with their terrace counts), and a mixed shelf with a delta.

## 4. Questions for the maintainer

1. Target: reef terraces and platform margins at regional scale, or global carbonate
   burial (pelagic production, the PNAS 2025 work), or both? This decides whether phase
   C3 (wave erosion) or a pelagic factory with a habitat map comes first.
2. Carbonate identity through provenance classes (passive label, already exact) plus
   dual lithology for grain size, or a true third lithology with its own porosity and
   transport coefficients (larger refactor of `stratplex`)?
3. Fuzzy-lite response tables, or the full Badlands fuzzy rule system?
4. Cellular automaton: out, or wanted for fine meshes?
