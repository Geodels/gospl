# Karst dissolution (design note, branch `feat/gravity-carbonates-karst`)

Opt-in chemical denudation of exposed carbonate: runoff dissolves the
carbonate lithologies of the stratigraphy, more where the flow concentrates,
and the dissolved carbonate can join the river dissolved load of the
groundwater geochemistry. Built on the carbonates (lithology piles,
`sst:`, the carbonate/geochemistry link).

## 1. Physics

Dissolution of limestone is limited by the water that reaches equilibrium
with it (White 1984; Kaufmann & Braun 2001). The denudation rate is the
runoff times the equilibrium carbonate content per unit volume of rock:

    D = alpha * R * F(A) * f_T          (m/yr of carbonate rock, bulk)

- `R = max(rain - evap, 0)` (m/yr), the local runoff;
- `alpha` (`rate`), the dissolution coefficient: volume of rock dissolved per
  volume of water, `Ceq / rho_rock` (~2e-5 to 1e-4: 50-250 mg/l of CaCO3);
  `alpha = 5e-5` and 1 m/yr of runoff give 50 m/Myr, the range of measured
  karst denudation (Ford & Williams 2007: 10-100 mm/kyr);
- `F(A) = max(Q / (R * A0), 1)^p` (`focus`, `focus_area`), the focusing of
  the dissolution by the flow: `Q` the discharge through the cell, so
  `Q / R` its contributing area for uniform rain, `A0` a reference area
  (default the mean cell area). `p = 0` is uniform lowering (pure runoff
  control); `p > 0` concentrates the dissolution along the drainage
  (dolines, karst valleys);
- `f_T`, an optional response table on temperature (`temperature`, read from
  the `sst:` field, the only temperature map of the model).

Only subaerial land dissolves: not the sea, not lakes (`_underLake`), not
ice-covered cells.

## 2. Lithology selectivity (`STRAMesh.dissolveStrat`)

The demand `d = D * dt` is taken from the carbonate lithology piles (default
every lithology from `carbonate:` / `strata: lithologies:` named `reef`,
`carb_sand`, `carb_mud`; or the `lithologies` list), layer by layer from the
top: in each layer the carbonate bulk `C` gives `min(C, d)`, split among the
carbonate piles in proportion. The insoluble part stays in place as residue,
and dissolution stops below a residue cover of `cover` metres (default
0.5 m): a clastic cover protects the limestone under it. Elevation, `cumED`
and `EbLocal` drop by the bulk removed. The provenance piles of each layer
scale with its bulk.

Conservation: bulk removed per lithology = solid removed / (1 - phi) per
layer; the solid volume is the dissolved volume (`karstDissolved`, per
lithology, owned nodes), reported in the run summary.

## 3. Geochemistry link

With the groundwater geochemistry and `river_load` on and `species:` set, the
dissolved solid volume times `solute_per_m3` (units of solute per m3 of rock,
the inverse of the duricrust `solid_volume` convention) is added to the
routing source of that species in `_routeRiverSolute`: it reaches the coast
(and the marine reservoir with `marine_coupling`, so the carbonate supply of
`carbonate_supply`). Without it the dissolved volume leaves the model and is
only reported.

## 4. Step position

After the first `flowAccumulation` of the step (discharge needed) and before
the groundwater step (so its river routing carries the same step's karst
load), before the erosion.

## 5. Phases

- K1: runoff dissolution, lithology-selective, temperature factor; output
  `karstDiss` (cumulative, float64, restart).
- K2: focusing (`focus`, `focus_area`).
- K3: geochemistry link (`species`, `solute_per_m3`).
- K4: tests (uniform lowering exact for `p = 0`, residue cover, lithology
  budget, off byte-identity, partition invariance) and an example (atoll
  with the long sea-level stack: emergent platforms dissolve during
  lowstands).
