# DESIGN — sediment provenance & source-to-sink attribution (Cu prospectivity)

Status: **both approaches implemented.** Standalone post-processor
(`gospl/analyse/provenance.py`); in-model tracers (Approach B, phases B0–B4,
opt-in `provenance:`) carry N source classes through the model's own
erosion/transport/deposition/stratigraphy — conservation-exact for any number of
sources, with exact per-class attribution in **both** sinks (marine and
continental pit/lake) via B2b. See §6 and `docs/tech_guide/provenance.rst`.

## 1. Problem

Given **source-rock polygons** (user shapefiles, each tagged with a rock type)
and **sink-basin polygons**, quantify — from a goSPL run — how much of the
sediment deposited in each sink derives from each source rock type, and where:

1. **per sink basin** — the % contribution of each source class;
2. **per pixel** (mesh node) in the basin — the source-mixture composition;

plus routing diagnostics (transport distance, erosion/deposition pathways) and a
**copper-prospectivity layer** weighting source classes by Cu fertility.

This is a *sediment-provenance* (source-to-sink attribution) problem layered on
goSPL's mass-conserving erosion/transport/deposition.

## 2. What goSPL provides

Per output step (`tout`) the mesh HDF5 holds: `coords`/`cells` (mesh), `elev`,
`erodep` (cumulative ED) and `EDrate`, `FA`, `fillFA`, `waterFill`,
`sedLoad`/`sedLoadF` (total + fine flux), and the stratigraphy
(`stratH`/`stratHf`/`phiS`/`phiF`/`stratK`/`stratZ`). Flow **directions are not
saved**, so the tool re-derives the routing network from the filled elevation
(`fillFA`/`elev`) + the mesh connectivity each step.

**Partitioned output.** goSPL writes one HDF5 per MPI partition per step
(`<base>.<step>.p<rank>.h5`, float32 datasets; partition meshes in
`topology.p<rank>.h5`), each holding only that rank's local nodes. Provenance
routing is global, so the tool (`GosplOutput`) reassembles each field onto the
**global input-mesh** ordering by mapping every partition's local coordinates to
the global vertices with a KDTree — the same local↔global map goSPL builds at
load (`tree.query(lcoords)`). Per-vertex `source_class`/`basin_id` are supplied
in that global ordering; shared ghost nodes carry agreeing (post-halo-sync)
values so the scatter is idempotent.

## 3. Two architectures

### A. Standalone particle/fraction tracking (this design + prototype)
Post-process the saved output: assign each mesh vertex a source class (from the
polygons), re-derive the per-step flow network, route the eroded material
downstream, and accumulate the deposited provenance per node/basin. No re-run
needed; flexible; gives distances and pathways. Approximate: it sees **net**
`erodep` per *output* step (not the compute `dt`), re-derives routing, and is
statistical/heuristic in its deposition rule.

### B. In-model provenance tracers (high fidelity — scoped §6)
Generalise the dual-lithology machinery (which already routes a *second* exact
sub-flux `vSedF` and stores per-layer composition `stratHf`) from 2 fractions to
**N source-provenance fractions**. Exact, mass-conserving, recycling-aware;
requires re-running with tracers on. Reuses proven code.

They are complementary: **A** for exploration on existing runs, **B** when
defensible percentages are needed.

## 4. Standalone tool — data model & algorithm

### Inputs

All per-vertex arrays are in the **global input-mesh** ordering (`npdata` mesh
`v`); build the labels with `classes_from_shapefile` and reassemble model fields
with `GosplOutput`.

| input | shape | required | meaning / how to obtain |
|---|---|---|---|
| mesh `coords` (`v`) | (N, 2\|3) | **yes** | global vertex coordinates; defines ordering + geometry |
| mesh `cells` (`c`) | (M, 3) | **yes\*** | triangular cells → adjacency (\*or pass a prebuilt `neighbours`) |
| `source_class` | (N,) int | **yes** | source-rock class per vertex, in `[0, n_classes)` — point-in-polygon of the source shapefiles; vertices outside every polygon get a **background/"other"** class (not −1) |
| `n_classes` | scalar | **yes** | number of source classes |
| goSPL output series | per step | **yes** | `elev` (routing surface) + `erodep`, fed to `step()` via `GosplOutput` |
| `area` | (N,) | no | Voronoi cell area (m²); default 1. Provide for correct volumes on irregular meshes (percentages are area-independent) |
| `basin_id` | (N,) int | **no** | sink-basin label per vertex, −1 = not a sink. Needed **only** for the per-basin % roll-up (`basin_percentages`); per-pixel provenance is produced everywhere without it. Build with `classes_from_shapefile(..., background=None)` |
| `cu_weight` | (C,) | no | copper fertility per source class; needed only for `cu_fraction` |

`source_class` is validated on construction (must be a full per-vertex array
with values in `[0, n_classes)`), so a mis-sized or −1-containing array fails
fast with a clear message rather than silently mis-attributing.

### State (carried across steps)
- `pile[node, class]` — provenance composition (m³ per class) of the
  *deposited* material currently stored at each node. Enables **recycling**:
  re-eroded sediment carries the pile's mixed provenance, fresh bedrock carries
  the node's `source_class`.
- `dep[node, class]` — cumulative deposited volume per class (the answer).
- `dist[node]` — flux-weighted mean source→deposition transport distance.

### Per output step (Δerodep = erodep(t) − erodep(t−1))
1. **Erode** at nodes with Δerodep < 0: eroded volume `E = −Δerodep·area`. Its
   provenance = the local `pile` composition for the pile portion (recycling),
   `source_class` (bedrock) for the remainder; decrement `pile`.
2. **Route** — build the downhill-flow graph from the filled elevation and sweep
   nodes high→low, pushing the per-class eroded flux to each node's
   receiver(s) (the deterministic continuum of particle advection — no
   Monte-Carlo noise). Routing is **multiple-flow-direction by default** (flux
   splits across all lower neighbours by slope^`flow_exp`, so a source can feed
   several basins partially — the realistic choice for provenance), with single
   steepest descent available (`routing='sfd'`). Carry a distance accumulator
   along the path (split with the flux under MFD).
3. **Deposit** at nodes with Δerodep > 0: capture `D = Δerodep·area` from the
   passing flux with the flux's current composition; add to `pile` and `dep`.
   Flux entering a sink basin / the sea settles there.

Conservation holds within a step by construction (eroded = deposited + outflow).

### Outputs
- **per basin**: `Σ dep[node∈basin, class]` normalised → % per source rock type.
- **per pixel**: `dep[node, :]` normalised → composition map (dominant source or
  full vector), rasterisable.
- **diagnostics**: mean transport distance per source→sink pair; pathway maps.
- **Cu layer**: `Σ_class cu_weight[class]·dep[node, class] / Σ dep` → Cu-sourced
  fraction per pixel/basin.

### Output format (saved per time step)
Results are written **per output step** so the provenance evolution is preserved
and can be viewed in ParaView alongside the goSPL output:
- a single **HDF5** `<prefix>.h5` with `/mesh/{coords,cells}`, `/steps`, and a
  `/step_<s>/` group per step holding `fractions` (npoints × n_classes),
  `dominant` (argmax class), `distance`, and `cu_fraction` (float32, gzipped,
  written incrementally so the full series is never held in memory);
- a **per-basin time-series CSV** (`<prefix>_basins.csv`: `step, basin,
  class0%, …`).
An XDMF temporal wrapper (so the HDF5 overlays the model output through time) is
a P1 add-on.

### Performance
The per-step **flow-graph build** (`downhill_edges`) is fully vectorised in NumPy
(edge slopes/weights via CSR run-length + `bincount`/`lexsort`). The downstream
**accumulation sweep** is inherently sequential (topological order, each node
depends on its donors) — O(N+E) per step.

The sweep is factored into a single Numba-`njit`-compatible kernel (`_sweep_impl`)
used both as the pure-Python reference and, when `numba` is installed, as the
compiled fast path (`method='auto'|'numba'|'python'`, default `auto`). Both give
identical results. Numba removes the Python per-node loop — the dominant cost at
scale — for a several-fold per-step speedup that grows with N (≈2.4× at 40 k
nodes in a quick bench; the vectorised edge build / erosion is shared, so small
meshes gain little).

**A sparse transport-with-loss solve was prototyped and rejected.** Recasting the
deposition as `(I − Wᵀ diag(1−f)) L_c = A_c` per class and factorising with
`splu` is exact only when *no* node is sediment-undersupplied (`f = D/T` ≤ 1) —
but the reconstruction *does* hit under-supply (net per-output-step `erodep` can
exceed the through-flux a coarse-cadence snapshot sees), where it diverges from
the `min(D, flux)` cap; and it was **slower** than the Python sweep at realistic
sizes (LU fill-in on the 2-D mesh matrix + two factorisations per step). The
Numba sweep is exact and faster, so it is the chosen path.

## 5. Copper prospectivity — scope boundary

The tool delivers the **physical provenance and routing** of detritus (which
source, how far, deposited where) — a necessary input to detrital/placer and
sediment-hosted Cu reasoning, **not** an ore-forming model. Real favourability
also needs hydrodynamic concentration (placers: heavy-mineral sorting by
energy/grain size — couple to the dual-lithology coarse/fine split) and, for
sediment-hosted Cu, facies/redox/diagenesis. Output is framed as a
*provenance-and-fertility proxy* (Cu-fertile detritus flux and accumulation),
explicitly upstream of the metallogenic step. The strongest favourability layer
combines **provenance × grain size × depositional setting**.

## 6. In-model provenance tracers (Approach B) — implemented (B0–B4)

Generalise dual lithology (`docs/DESIGN_DUAL_LITHOLOGY.md`) to N provenance
classes:
- **State**: `stratP[node, layer, class]` per-layer provenance fractions
  (mirrors `stratHf`); `provFrac[node, class]` of the routed flux (mirrors
  `fineFrac`).
- **Erosion**: split eroded solid by the consumed layers' provenance (mirrors
  `erodeStrat`'s fine split), bedrock contributing its `source_class`.
- **Transport**: route each class sub-flux through the same flow operator as the
  total (linear, exact) — N extra `_solve_KSP` solves (or one block solve);
  reuse `_getSedFlux`/`_moveDownstream` threading like `vSedF`.
- **Deposition**: deposit with the flux composition; write `stratP`.
- **Conservation**: per-class guard (the analogue of `test_dual_fine_conservation`).
- **Cost/scope**: N is user-set (a handful of source provinces). Heaviest part is
  N sub-fluxes per step; acceptable for the target use. Output `stratP` →
  provenance % per pixel/basin directly, recycling-aware, exact.

Deliver B opt-in (`provenance:` YAML block), inert when absent (byte-identical),
mirroring the dual-lithology delivery discipline. Provenance is a **passive
label** — unlike lithology it has no erodibility/diffusivity/porosity feedback
and no depositional sorting, so there is no `_surfaceLithoK/D` or
`_pitFineFraction` analogue; it simply rides the existing total-sediment routing.

### Phases (in-model, branch per phase, PR into `dev`)

- **B0 — foundation** ✅ *(done)*: `_extraProvenance` parser (`provenance:` →
  `provOn`, `provNb`, source-class `uniform`/`source` map; requires
  `stratNb > 0`; note `cu_weight` is an analysis-tool argument, not a model
  input — see §"Cu layer"); state `stratP[node, layer, class]` (per-class layer thickness,
  Σ = `stratH`) seeded to the bedrock `source_class` in `readStratLayers`; routed
  sub-flux vecs `vSedP[c]` + `provFrac`/`depoProvFrac` + `_provEroded`/
  `_provDeposited` diagnostics in `sedplex`; `destroy_DMPlex` registration.
  Passive (no hooks yet) ⇒ provenance-on is byte-identical to off. Tests:
  `test_provenance_opt_in`, `test_provenance_seeding`.
- **B1 — erosion split** ✅ *(done)*: `erodeStrat` mirrors every `stratHf`
  operation per class — the eroded sediment is split by the consumed layers'
  `stratP` (the bedrock sentinel → the node's `source_class`), producing the
  per-class eroded rate `provEro` (Σ over classes == the total uncompacted
  erosion) and accumulating `_provEroded`; `stratP` is reduced and re-normalised
  so Σ over classes == `stratH`. Tested in isolation (`test_provenance_erosion_split`).
- **B2 — transport** ✅ *(through-flux)*: route the N sub-fluxes `vSedP[c]` through
  the upstream-integration operator `fMati` in `_getSedFlux` (Σ over classes ==
  total flux, exact); `provFrac = vSedP/vSed` is the arriving composition,
  `depoProvFrac` (no sorting bias — passive) is what `deposeStrat` lays down.
- **B3 — deposition** ✅: `deposeStrat` adds `depo · depoProvFrac` to the new
  layer's `stratP` (keeps Σ over classes == `stratH`) and accumulates
  `_provDeposited`. Guarded by `test_provenance_conservation` — with a single
  source, every layer stays 100 % that class after a full run (class-0 leakage ==
  0; `stratP` partitions `stratH` to ~3e-8).
- **B2b — marine composition** ✅ *(done)*: `seaChange._marineProvFraction` sets
  the marine deposit's source composition to the **basin-delivered mix**
  (`ff_mar[c] = Σ provFrac[c]·sedFlux / Σ sedFlux`, uniform — no depth bias, the
  passive-label analogue of `_marineFineFraction`). This was the dominant sink
  (≈79 % of deposition in the test) and was the weakest spot of B2 (marine-only
  nodes had a near-zero through-flux composition). With it the recorded
  composition matches the eroded supply ratio to ~1e-6 and the partition becomes
  machine-exact (Σ over classes == `stratH`, ~1e-24). Domain-uniform, the same
  standard as the dual marine fraction.
- **B2b — continental pit cascade** ✅ *(done)*: pit/lake deposits now carry each
  pit's **cascade-retained** source mix rather than the through-flux composition.
  `_distributeSediment` builds a per-class sub-flux `vSedP[c]·dt` and threads it
  through `_moveDownstream` in lockstep with the total: each pit retains its mix
  proportionally (`ret_frac = pitVol/inV`, so `ret_prov = inVp·ret_frac`) and
  overspills the rest, routed through the same flow matrix as the total (linear),
  so **downstream-lake chains** mix exactly. The per-pit retained provenance
  accumulates in `_pitRetProv`; `_pitProvFraction` (called from `_updateSinks`)
  sets `depoProvFrac[in_pit]` uniformly to `_pitRetProv/depo` — a passive-label
  analogue of `_pitFineFraction` but **with no depocenter bias** (provenance is a
  label, not a grain size), hence simpler than the dual case. Σ over classes of
  `_pitRetProv` == the pit's retained volume, so `depoProvFrac` stays summed-to-1
  and `stratP` partitions `stratH` machine-exactly. Guarded by
  `test_provenance_pit_fraction` (the uniform-mix invariant) plus the existing
  conservation tests (now run with pit attribution active).
- **B4 — advection + I/O + restart** ✅: `stratalRecord` advects each class's
  `stratP[:,:,c]` with the same `strataonesed` interpolation as `stratHf`
  (re-normalised so Σ over classes == `stratH`); `_outputStrat` writes the
  `(lpoints, layers, classes)` `stratP` to the stratal HDF5 and the restart path
  restores it. Provenance per pixel/per layer reads directly off `stratP`.
  Tested via the I/O round-trip (`test_provenance_output_io`); the advection
  branch mirrors the tested `stratHf` path (a horizontal-tectonics fixture would
  exercise it directly).

## 7. Phasing

- **P0** — standalone core (this prototype): neighbour build, per-step receiver
  re-derivation, fraction-routing provenance with recycling, per-basin % +
  per-pixel composition + distance + Cu layer. Synthetic-mesh tests.
- **P1** — I/O: read goSPL HDF5 series; shapefile → `source_class`/`basin_id`
  (lazy `geopandas`); raster/CSV/plot export.
- **P2** — discrete-particle mode for full path-length *distributions* (the
  fraction router gives means); pathway visualisation.
- **P3** — Approach B in-model tracers (§6) for conservation-exact provenance.
- **P4** — Cu favourability: combine provenance × dual-lithology grain size ×
  depositional energy.

## 8. Caveats (carried into the docs)

- Output cadence, not compute `dt`: sub-step erode→redeposit is invisible to A.
- Routing is reconstructed from saved filled elevation, not goSPL's live network.
- A is statistical/heuristic (exact volumes ⇒ B); good for % and pathways.
- Provenance ≠ prospectivity (§5).
