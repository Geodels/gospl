# DESIGN: Level-B geochemistry extensions

Two opt-in extensions of the Level-B solute/duricrust geochemistry
(`DESIGN_WATERTABLE_GEOCHEM.md`, G0–G6+G5b), designed as separate increments on
top of the shipped feature:

1. **Spatial per-species weatherability** — let *lithology* control which solute
   species each region yields (mafic → Fe/silica, carbonate platform →
   carbonate), instead of one global weatherability per species.
2. **River dissolved-load coupling** — route the groundwater-exported solute
   *down the surface drainage network* to the shoreline, the natural continuation
   of the baseflow seepage the model already computes.

Both preserve the standing invariants: **opt-in**, **byte-identical when off**,
partition-exact, MPI-safe, no new hard dependency.

---

## Extension 1 — spatial per-species weatherability (lithology → chemistry)

> **STATUS: DONE** — all three forms built: **(a)** direct per-vertex maps,
> **(b)** a standalone `lithology: [file, key]` map + per-(class, species) table,
> and **(c)** the same table gathered by the provenance class — static
> `source_class` **or** dynamic `surface_class` (top stratigraphic layer,
> re-resolved each step; §1.8). Parser keeps `gwGeoWeather` raw + adds
> `gwGeoWeatherByClass` / `gwWeatherFrom` / `_gwLithoMap`;
> `gwplex._resolveGeoWeather` resolves to scalar-or-`(lpoints,)` per species
> (label = the `lithology:` map if given, else the provenance class); the
> use-site is one line. Guards `test_geochem_spatial_weatherability` (c) +
> `test_geochem_lithology_map` (b) + `test_geochem_surface_lithology` (dynamic);
> full `tests/` 154 passed. **No open refinements.**

### 1.1 Motivation & the current gap
Today the solute species are defined by **global scalars** — one
`weatherability`, `c_sat`, `precip_rate`, `solid_volume` per species, applied
everywhere (`inputparser.py:2162`, resolved to a length-`n_species` numpy array
in `gwplex.py:154`, used as the scalar `self.gwGeoWeather[k]` in the dissolution
step `Drate = np.where(subaerial, self.gwGeoWeather[k] * W, 0.0)`,
`gwplex.py:961`). So every subaerial cell dissolves **all** species with the
same mix; only the dissolution *rate* `W` varies in space (climate + Arrhenius
temperature + soil), and it scales all species together.

`source_class` (provenance) *is* per-vertex and can come from a map or the
initial stratigraphy (`stratplex._initProvenance`, `stratplex.py:487–503`), but
it is only a provenance **tag** (it sets `crust_source`); it does **not** change
which species a region contributes. There is therefore **no lithology→chemistry
link**: a mafic terrane and a limestone platform yield identical solute mixes.

**Goal:** make `weatherability` (the primary knob) a **per-vertex, per-species**
field, so a lithology map or the initial-stratigraphy composition decides which
species each region weathers out.

### 1.2 Precedent to mirror
goSPL already resolves scalar-or-`[file, key]` per-vertex maps in exactly this
pattern; the extension reuses it verbatim:

- **`duriWeatherability`** (the Level-A single weatherability): scalar **or**
  `[file, key]`, lazily resolved to `self._duriWeatherArr =
  np.load(wab[0] + ".npz")[wab[1]][self.locIDs]` (`gwplex.py:590–596`).
- **`gwInfiltration`**: `_gwInfilMap = [file, key]` resolved in `__init__`
  (`gwplex.py:140–145`) — subset to the local partition via `self.locIDs`.
- **`source_class`**: `uniform` scalar or `[file, key]` map
  (`stratplex.py:487–490`), already per-vertex and partition-exact.

### 1.3 Input design — three ways to specify the field
Ordered from most explicit to most convenient; **(a) is the base mechanism**,
(b)/(c) are thin conveniences layered on it.

**(a) Direct per-species map** — each species' `weatherability` may be a
`[file, key]` per-vertex field instead of a scalar:
```yaml
geochem:
  species:
    - {name: carbonate, weatherability: [litho, carb_wab], c_sat: 1.0, ...}
    - {name: silica,    weatherability: [litho, sil_wab],  c_sat: 2.0, ...}
```
Most explicit; the user supplies one weatherability field per species.

**(b) Lithology-class map + a per-(class, species) table** — one per-vertex
integer lithology map plus a small matrix, compact when many species share a few
rock types:
```yaml
geochem:
  lithology: [litho, rock_class]          # per-vertex int, 0..n_litho-1
  species:
    - {name: carbonate, weatherability_by_class: [1.0, 0.0, 0.2]}  # per rock class
    - {name: silica,    weatherability_by_class: [0.1, 0.8, 0.6]}
```
Resolve per node: `wab_k[node] = table[k][litho_class[node]]` (a gather).

**(c) Driven by the initial stratigraphy / provenance** — reuse the existing
per-vertex `source_class` (already map-driven and seeded into `stratP`) as the
lithology label, so **no new input** is needed when provenance is already on:
```yaml
geochem:
  weatherability_from: source_class       # reuse provenance regions as lithology
  species:
    - {name: carbonate, weatherability_by_class: [1.0, 0.0]}
    - {name: silica,    weatherability_by_class: [0.0, 1.0]}
```
`wab_k[node] = table[k][source_class[node]]`. This is (b) with the lithology map
supplied by provenance — the tightest tie to "where the rock is," and it makes
provenance's `crust_source` and the species mix mutually consistent.

**Recommendation:** **(a)** is the primitive (a per-species per-vertex array),
**(c)** the ergonomic default for provenance runs (no extra file), and **(b)**
the same table gathered by a standalone `lithology:` map for runs without
provenance. All three are implemented; the label for (b)/(c) is the `lithology:`
map when given, else `source_class`.

### 1.4 Parsing (`inputparser._readGroundwater`)
`gwGeoWeather` stays a per-species **list**, but each entry becomes **scalar or
`[file, key]`** (don't coerce to `float` unconditionally — keep lists):
```python
self.gwGeoWeather = [s.get("weatherability", 1.0) for s in species] or [1.0]
```
Add, when present: `self.gwGeoWeatherByClass = [s.get("weatherability_by_class")
...]`, `self._gwLithoMap = geo.get("lithology")`, `self.gwWeatherFrom =
geo.get("weatherability_from")`. All default `None` in the `except` branch so an
absent block is byte-identical.

### 1.5 Resolution (`gwplex.__init__`, mirroring `_duriWeatherArr`)
Build `self._gwGeoWeatherArr` — a length-`n_species` list where each entry is
**either a float (scalar path, unchanged) or an `(lpoints,)` array**:
- entry is `[file, key]` → `np.load(file+".npz")[key][self.locIDs]` (option a);
- `weatherability_by_class` + a lithology label → gather
  `table[k][label[self.locIDs]]` where `label` is the lithology map (b) or
  `source_class` (c). `source_class` is allocated by `_initProvenance`, so
  resolve **after** it (guard: provenance must be on for option c).
- otherwise `float(...)` (the current behaviour).

### 1.6 Use site (one line, `gwplex.py:961`)
```python
wab_k = self._gwGeoWeatherArr[k]                 # scalar OR (lpoints,) array
Drate = np.where(subaerial, wab_k * W, 0.0)      # numpy broadcasts either
```
Everything downstream (`diss`, `Deff`, transport, precipitation, typing,
provenance) is already per-node and needs **no change**. A zero-weatherability
region simply yields none of that species (e.g. no carbonate from mafic rock) —
the correct behaviour.

### 1.7 Invariants & edge cases
- **Byte-identical when off:** absent maps → the scalar path is untouched.
- **Partition-exact:** every field is subset by `self.locIDs`, like every other
  goSPL map; no cross-partition coupling, no new collective.
- **MPI:** resolution is rank-local in `__init__`; nothing added to the solve.
- **Source pool:** `gwSourcePool` (seeded `1e6·area`) is unchanged — spatial
  weatherability changes the *rate*, not the reservoir. (Optionally scale the
  pool by lithology later; out of scope here.)
- **`c_sat` / `precip_rate` maps:** the same mechanism could map these too, but
  weatherability is the lithology→chemistry primary; keep the others scalar in
  v1 to limit surface area.

### 1.8 Dynamic surface lithology — DONE
The prescribed forms (a)/(b) and `source_class` (c-static) fix the weatherability
to the **present-day bedrock**; but as erosion exhumes deeper stratigraphy — or
deposition buries the surface under transported sediment of a different
provenance — the real surface rock type changes. `weatherability_from:
surface_class` closes this: `_surfaceSourceClass` derives the per-node label each
step from the **dominant provenance class of the top non-empty stratigraphic
layer** (`stratP`), using the same top-layer scan as `_recordInduration`; the
dynamic form is re-resolved every step (not cached). Needs provenance (the
per-layer class record `stratP`); falls back to the bedrock `source_class` where
a column has no sediment. Guard `test_geochem_surface_lithology`.

### 1.9 Testing
- Fixture: 2 species + a lithology split (region A weatherability `[1, 0]`,
  region B `[0, 1]`). Assert `crust_type` == species 0 in A, species 1 in B, and
  that each region's off-species crust is ~0.
- Byte-identity: a uniform map (all-ones) reproduces the scalar run exactly.
- np=2: the field is partition-consistent (halo values match).

### 1.10 Effort
Small — parser (scalar-or-map, ~3 keys), `__init__` resolution (~15 lines
mirroring `_duriWeatherArr`), a one-line use-site change, one fixture + test.

---

## Extension 2 — river dissolved-load coupling

> **STATUS: DONE — incl. all §2.6 refinements.** Opt-in `geochem: river_load:
> true` → `gwRiverLoad`. `_routeRiverSolute` routes the seepage export down the
> flow network **per species** (RHS `gwSoluteFluxSp[:, k]` on the cached `fMati`),
> called at the end of `updateGroundwater` (matrix fresh — no reordering).
> Refinements built: **per-species routing** (`riverSoluteSp` /
> `riverSoluteToOceanSp`); **in-transit reactions** — a first-order per-species
> loss `river_decay = κ` via the shifted operator `(I − Wᵀ + κ I)L = s`
> (`fMati.shift(κ)`), lost mass `riverSoluteLost`; **marine coupling**
> (`marine_coupling`) — delivered coastal flux accumulates into a per-species
> reservoir `marineSolute` + per-node `marineSoluteInput` output. **Exact
> conservation** per species via `Σs = (fMatiᵀ·1)·L + riverSoluteLost`
> (`test_geochem_river_load` conservative, `test_geochem_river_species_reactions_marine`
> full); np=2 confirmed (per-species exact, decay loss, marine accumulation).
> Full `tests/` 155 passed.

### 2.1 Motivation & the current gap
The solute the model dissolves is exported to the ocean via **groundwater
seepage / baseflow** (G3): `gwSoluteFlux` (per-node discharge-to-surface rate,
m³/yr) and the lumped scalar budget `gwOceanFlux`
(`gwplex.py:1006–1008`). But that export is booked as reaching the ocean
**directly** — it never travels **down the river network** the way sediment and
water discharge do. In nature the baseflow-delivered solute rides the surface
drainage to the coast (the dominant pathway for chemical-weathering products to
the sea; the dissolved load), gaining more solute along the way and delivering
it at specific shoreline points. This extension routes the per-node seepage
export along the **same drainage network** goSPL already uses for water and
sediment.

### 2.2 Why it is *much* simpler than sediment
A conservative dissolved load is a **passive tracer** on the flow graph:
- **no deposition / no capacity limit** — solute doesn't build up a bed the way
  sediment fills pits, so the whole `_moveDownstream` / `_distributeSediment`
  pit-overspill cascade (`sedplex.py:190–522`) is **not needed**;
- **no coarse/fine split** — one scalar (or one per species);
- the flow matrix is built on the **filled** topography (`waterFilled`), so a
  single accumulation solve already routes solute *through* filled lakes to
  their spill point and onward to the coast.

So the core is exactly the **sediment-flux accumulation solve** minus the
cascade — the `_getSedFlux` pattern (`sedplex.py:103–188`):
```python
# (I - Wᵀ) L = s   — one implicit solve, reusing the cached flow matrix
self._solve_KSP(False, self.fMati, sourceG, self.riverSoluteG)
```
where the RHS `s` is the per-node seepage export and `L` the accumulated river
dissolved load (m³/yr) at every node — increasing downstream, delivered at the
coast.

### 2.3 Hooks to reuse (from the routing map)
| Role | Existing object | Reference |
|---|---|---|
| Flow matrix `(I − Wᵀ)` (filled topo) | `self.fMat`, cached `self.fMati` | `flowplex.py:814`, `matrixFlow` `:429–475` |
| Accumulation solve | `self._solve_KSP(False, self.fMati, b, x)` | `sedplex.py` `_getSedFlux` `:103–188` |
| Per-node source | `self.gwSoluteFlux` (m³/yr export) | `gwplex.py:1008` |
| Coast / domain exit | `self.seaID`, `self.outletIDs` | `flowplex.py:526–537` |
| Closed basins (trap) | pits with no ocean spill (`pitIDs`, `lFill`) | `flowplex.py:576–710` |

`fMat` is **partition-stable** (the parallel-perf work) so the routed load is
partition-invariant, and `_solve_KSP` is the established MPI-collective pattern.

### 2.4 Design
**Source.** After `_updateSolute` fills `gwSoluteFlux` (per-node, summed over
species) and after `flowAccumulation` has built `fMat`/`fMati` for the step, use
`gwSoluteFlux` as the RHS (already a volumetric rate — no `pointwiseMult` by area
needed, unlike sediment's per-node thickness).

**Solve.** One `(I − Wᵀ) L = s` on `fMati` → `riverSolute` (local/global vecs
`riverSoluteL`/`riverSoluteG`, `(lpoints,)` field). `fatal=False`, `seed=True`
(a substochastic accumulation, like sediment). No cascade loop.

**Delivery & accounting.** `L` at `seaID`/`outletIDs` is the solute leaving the
continent at that point; sum over coastal exits = **total dissolved flux to the
ocean, spatially resolved**. A global reduction gives `riverSoluteToOcean`, which
(absent closed-basin trapping) equals the existing lumped `gwOceanFlux` — a
built-in **conservation cross-check**. Solute routed into a **closed basin** (a
pit that never spills to the sea — an endorheic/evaporite basin) terminates
there rather than reaching the coast: physically correct, and it's exactly what
the filled-topo matrix does (no artificial removal needed).

**Per-species (optional).** To carry each species down the river separately
(dissolved-species delivery to the ocean, e.g. riverine silica vs alkalinity),
route each species' seepage export as its own RHS on the **same** `fMati` —
`n_species` cheap solves reusing one matrix (the G5 provenance trick). v1 routes
the **total**; per-species is a one-loop add.

**Output.** New `riverSolute` field (HDF5+XDMF, gated on the flag); verbose
per-step `riverSoluteToOcean` total. Keep `gwOceanFlux`/`soluteflux` (the budget
scalar and the seepage source) — `riverSolute` is the *routed* field on top.

### 2.5 Sequencing — already favourable, no reordering needed
The solve needs **both** `gwSoluteFlux` and `fMati` for the same step, and the
existing pipeline already provides them in the right order: `model.py`
(`runProcesses`) runs `_FAMesh.flowAccumulation` (`model.py:400`, which caches
`self.fMati` at `flowplex.py:814`) **before** `_GWMesh.updateGroundwater`
(`model.py:407`, which fills `gwSoluteFlux`). So `_routeRiverSolute()` can be
called at the **end of `updateGroundwater`** (right after `_updateSolute`):
`fMati` is fresh and `gwSoluteFlux` is just computed — **no lag, no pipeline
reordering.** (The second `flowAccumulation` before sediment, `model.py:432`,
rebuilds the matrix on post-erosion topo; the river-solute routing intentionally
uses the pre-erosion `fMati` consistent with the groundwater step.)

### 2.6 Invariants & faithfulness caveats
- **Opt-in / byte-identical:** gated on e.g. `geochem: river_load: true`; no
  extra solve, no new field when off.
- **Reuses `fMati`:** no new matrix, no new preconditioner; one extra KSP solve
  per step (or `n_species`), cheap next to the flow/sediment solves.
- **Conservative & partition-invariant** by construction (same operator as water
  discharge).
- **In-transit reactions (BUILT):** a per-species first-order loss
  `river_decay = κ` (in-channel precipitation / uptake) via the shifted operator
  `(I − Wᵀ + κ I)L = s` (`fMati.shift(κ)`); the lost mass `κ·Σ L` is a per-species
  diagnostic `riverSoluteLost`. `κ = 0` (default) is the conservative passive
  routing. NOTE the loss is a **per-node** first-order retention, so it is
  mesh-resolution-sensitive (a distance-weighted form using `distRcv` is the
  further refinement) — the rate is a tunable knob, not a calibrated constant.
- **Marine coupling (BUILT):** the delivered coastal flux (`riverSoluteToOceanSp`)
  accumulates into a per-species reservoir `marineSolute` (m³, integrated over
  time) with a per-node `marineSoluteInput` output. It is a **delivery reservoir**
  (the ocean's weathering-derived inventory), NOT a marine reactive-transport
  model — that (mixing, carbonate system, sedimentation) remains out of scope.

### 2.7 Testing
- Fixture: a ramp draining to one open coast; inject a uniform `gwSoluteFlux`.
  Assert (a) `riverSolute` **increases downstream** (accumulation), (b) the
  coastal-exit sum **== total injected** (conservation, no closed basin), (c)
  `riverSoluteToOcean` matches `gwOceanFlux` to tolerance.
- Closed-basin variant: an interior pit that never spills — assert its solute is
  **trapped** (coastal exit < injected by the trapped amount).
- np=2: `riverSolute` matches the np=1 field (partition invariance of `fMat`).

### 2.8 Effort
Small–moderate — a `_routeRiverSolute()` (~30 lines mirroring `_getSedFlux`), two
scratch vecs + `destroy_DMPlex` registration, one output field, the ordering
hook, a fixture + tests. No new matrix or solver machinery.

---

## Suggested build order
1. **Extension 1** first (self-contained, no pipeline ordering, small): parser →
   `__init__` resolution → one-line use-site → fixture/test. Ships the
   lithology→chemistry link.
2. **Extension 2** next (needs the `model.py` step-order check + a scratch-vec /
   `destroy_DMPlex` registration): `_routeRiverSolute` → output → fixture/tests.

Each is a separate commit on `feat/watertable-geochem` (or a fresh branch),
following the standing invariants and the AGENTS milestone/checklist.

---

## Post-processing & outputs (DONE)

**Per-species outputs.** With `n_species > 1`, `outmesh` writes per-tracer
HDF5+XDMF fields named by species — `solute_<name>` (concentration),
`crust_<name>` (crust contribution), `soluteflux_<name>` (groundwater seepage
export) and, with `river_load`, `riverSolute_<name>` — alongside the aggregated
totals (each per-species set sums to its total). `gwSoluteFluxSp` is allocated in
the general geochem path so `soluteflux_<name>` is available regardless of
`river_load`.

**Basin extraction (`gospl-catchment`).** The geochem fields grid automatically
(`gospl-grid` rasterises every step field). `analyse/catchment.basin_solute_flux`
then extracts, per drainage basin, the **solute outlet** (cell of maximum total
solute flux) with `basin,lon,lat,val` **plus one column per species** (each
species' flux at that outlet; they sum to `val`). The total field auto-detects
`riverSolute` → else `soluteflux`; per-species fields are found by the
`<total>_` prefix. `basin_outflow` adds a `"solute"` entry when present;
`catchment_flux` writes `solute{t}.csv`. This is the dissolved-load analogue of
the water / sediment river-mouth extraction.

