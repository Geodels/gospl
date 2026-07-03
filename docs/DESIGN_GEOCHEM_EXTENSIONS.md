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

> **STATUS: DONE** (built as designed — forms **(a)** direct per-vertex maps +
> **(c)** provenance-driven table; **(b)** the standalone lithology-class map is
> not wired, since **(c)** covers the same need by reusing `source_class`).
> Parser keeps `gwGeoWeather` raw + adds `gwGeoWeatherByClass` / `gwWeatherFrom`;
> `gwplex._resolveGeoWeather` resolves lazily to scalar-or-`(lpoints,)` per
> species; the use-site is one line. Guard `test_geochem_spatial_weatherability`;
> full `tests/` 151 passed. Static-surface-lithology (§1.8) remains the one noted
> refinement.

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

**Recommendation:** implement **(a)** as the primitive (a per-species per-vertex
array), and expose **(c)** as the ergonomic default for provenance runs (it
needs no extra file). **(b)** is (c) generalised to a standalone lithology map;
add it if a lithology map exists independent of provenance.

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

### 1.8 Known limitation — static surface lithology
The map is the **present-day surface** lithology and is **static**: as erosion
exhumes deeper stratigraphy the real surface rock type changes, but v1 keeps the
weatherability field fixed. A faithful refinement would derive the per-node
surface lithology each step from the **top non-empty stratigraphic layer's**
`source_class` (the same top-layer scan `_recordInduration`/`_surfaceComposition`
already do) — deferred, noted here as the natural next step.

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
- **In-transit processes (v1 omits):** real rivers also gain/lose solute in
  transit — extra in-channel weathering, evapoconcentration, carbonate
  precipitation, biological uptake. v1 is a **conservative passive** routing (the
  right level at My steps, where in-channel residence is instantaneous). A
  first-order in-channel loss/gain would add a diagonal reaction term to the
  operator (like the seepage sink in the groundwater solute solve) — deferred.
- **Marine coupling (out of scope):** the delivered coastal flux could seed a
  marine solute / alkalinity state; noted, not built.

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

