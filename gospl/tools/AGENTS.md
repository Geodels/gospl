# Input parsing, outputs, forcing, diagnostics (gospl/tools/)

Read before editing `inputparser.py` (high-risk: parameter parsing, forcing DataFrames), `outmesh.py`, `addprocess.py` or `zprobe.py`.

This file is part of the goSPL agent guide. Read the root `AGENTS.md` first; it holds
the rules that apply everywhere (MPI contract, solver lifecycles, scratch Vecs,
conventions, commit checklist). The text below was moved verbatim from the root
file in the 2026-10 split, so dates and `## ...` cross-references still refer to
the original section names (most now live in the files listed in the root
`AGENTS.md` > Map).

## The `_extra*` methods are mandatory continuations
NOT optional parsers. Each sets attributes required by other modules. Never delete or rename without following the full call chain in `inputparser._readDomain/_readTime/_readHillslope/_readCompaction/_readFlex/_readOrography/_readIce`.

- `_extraDomain` (inputparser.py:189) → `seaDepo`, `overlap`, `dataFile`, `nodep`, `strataFile`; calls `_extraDomain2`.
- `_extraDomain2` (:229) → `advscheme`, `radius`, `gravity`.
- `_extraHillslope` (:475) → `nlK`, `clinSlp`, `tsStep`, `Gmar`, `offshore`, `nl_pit_volume/depth/K/inlet_bias`, `marineSolver`(`ts`|`picard`)/`picardSub`/`picardIts` (opt-in lagged-diffusivity marine/lake solver).
- `_extraStrata` (called from `_readCompaction`) → `stratLith` (dual-lithology master opt-in), `phi0c/z0c` (coarse porosity curve, defaults to compaction `phi0s/z0s`), `phi0f/z0f` (fine), `fine_k_factor` (fine erodibility multiplier), `bedrock_coarse_frac`, `fine_efficiency`, `pit_inlet_bias_coarse/fine`, `fine_diff_factor` (fine diffusivity multiplier), `bedrock_sentinel`. All the dual-lithology params are inert while `stratLith` is False; forced False if `stratNb == 0`. EXCEPTION: `bedrock_sentinel` is independent of `stratLith` — it inserts a dedicated infinite-bedrock layer (1e6 m, `bedrock_coarse_frac` split) BENEATH file-provided initial layers in `readStratLayers` (off by default; only consulted when an `npstrata` file is supplied — the no-file path always builds a sentinel). `readStratLayers` validates the `npstrata` file up front via `_checkStrataFile` (required `strataH`/`strataZ`/`phiS` present, every layer array shaped `(mpoints, n_layers)`, rank-0 warning if dual-on but `strataHf` absent → all-coarse) and emits a rank-0 `-v` setup summary via `_logStratInit` (both gated by `getattr(self,"verbose",False)` so bare `STRAMesh.__new__` tests stay silent). See `docs/DESIGN_DUAL_LITHOLOGY.md`.
- `_readFlex` → `flexure.method` is `'fem'` (flat, default) or `'global'` (spherical); `'FD'`/`'FFT'`/gFlex were removed. `thick`/`rhoc`/`rhoa`/`young`/`nu`; `ninterp` is global-only (DH-grid KDTree). `regdx` is no longer a key anywhere — the regular grid has been removed entirely (orography is now solved on the mesh too).
- `_extraFlex` (:1656) → `flex_res_deg`, `flex_bcN/S/E/W`, `flex_max_iter`/`flex_tol`/`flex_relax` (varying-Te **global** iterative-solve controls), `flex_interval` (apply flexure every N steps; load accumulates — the eroder snapshots `hOldFlex` only at interval starts, gated by `self.flexCount % flex_interval == 0`, and `model.runProcesses` applies flexure at interval ends).
  - **`'fem'` (flat)** — parallel FV biharmonic on the DMPlex (`_cmptFlexFEM`/`_buildFlexFEM`): single-field `[Lm·diag(D)·Lm + Δρg·I] w = q`; cached operator+factorisation reused each step (serial PETSc LU / parallel **MUMPS** — GMRES+GAMG does NOT converge on the stiff biharmonic, so direct is the default); varying Te = one solve. BCs `0Slope0Shear`/`Mirror` (natural FV) + `0Displacement0Slope` (clamped, Dirichlet `zeroRows`); `0Moment0Shear`/`Periodic` not implemented. Cached `_flexLm/_flexA/_flexKSP` in `destroy_DMPlex`.
  - **`'global'`** — serial on rank 0 (pyshtools SH); the earth-model scaling wall — but with the post-tuning config (`res_deg 0.5`, `interval>1`, capped `maxIter`, warm-start) it is ~0.02–0.07 s/call (no longer a wall). Varying-Te (`temap`) branch is iterative (Anderson), warm-started.
- `_readOrography` → `wind_speed`, `wind_dir` only (the spectral-only `regdx`/`latitude`/`nm`/`hw` were removed). `_extraOrography` → `oro_cw` (= `ref_density`·`moist_lapse_rate`/`env_lapse_rate`), `oro_conv_time`, `oro_fall_time`, `oro_precip_*`, `rainfall_frequency`. Orography is now solved on the mesh (`cptOrography`, see milestone) — no regular grid.
- `_extraIce` → legacy uniform `iceT`/`elaH`/`iceH` time interpolators (scalar or `evol` CSV). `_readIce` also sets the diagnostic params (`ice_slide`, `ice_glen`, `ice_accum_factor`, `ice_accum_max`, `ice_meltfac`, `ice_melt_conserve`, abrasion `ice_Kg`/`ice_abr_l`/`ice_Kl`/`ice_lat_l`, `ice_till_on`, `ice_till_route`) — there is **no `flow_model` selector and no `sia`/`hinit` block** (only the diagnostic model exists). `hela`/`hice`/`hterm` may be uniform scalars, per-vertex maps, or a `glaciers` time series (`_buildIceSeries` → `_iceTimeSeries`, resolved each step by `mesher._updateIce` into `elaMesh`/`iceMesh`/`termMesh`); `hterm` defaults to the `TERMINUS_UNSET` sentinel → sea level. See `## Ice sheet` and `docs/DESIGN_ICE_SHEET.md`.

## YAML parsing helpers (use `_get_param` / `dict.get`)
`tools/inputparser.py` no longer uses bare `try/except KeyError` for default-on-miss. Two patterns are canonical:

- **`self._get_param(*keys, default=None)`** — safe traversal into `self.input`. Returns `default` if any key in the chain is missing. Use when accessing `self.input` directly without pre-extracting a section.
- **`section_dict.get(key, default)`** — Python builtin, used when the section dict has already been pulled into a local variable (e.g. `domainDict = self.input["domain"]`). Preferred over `_get_param` in that case — it makes the data flow more visible.

Both replace the old `try: self.x = dict[key]; except KeyError: self.x = default` boilerplate.

**Do NOT re-introduce that boilerplate.** ~85 blocks were collapsed in 2026-06; the helper docstring (`_get_param` at the top of `ReadYaml`) explains the convention.

**Out of scope for these helpers** (keep the existing `try/except KeyError`, do NOT convert):
1. **Required keys** that print a user-facing diagnostic and raise (16 sites — e.g. `K` in `spl`, `start`/`end`/`dt` in `time`). The except handler is part of the user-facing error contract.
2. **Outer-section blocks** that set multiple defaults + flip a feature flag (e.g. `flexOn=False` when `flexure` is missing — 15 sites). Inverting these to `if "section" in self.input` is a real control-flow refactor, not a 1-line swap.
3. **Other side-effect handlers** (5 sites — the `tout` print, the dual-assignment `soilfile/self.soilFile` pair, the convoluted `sea`/`position`/`curve` nested fallback, the `latitude` bounds-check inside the try).
4. **NPZ archive accesses** (4 sites — `mdata[key]` lookups in `_isKeyinFile` and the rainKey/sedKey/teKey blocks). `mdata` is an `np.load()` result, not `self.input`.

All 33 such sites carry a `# TODO-REFACTOR: complex except, needs manual review` comment with the specific reason.

## Forcing DataFrame layout contract
Consumers use `df.at[nb, col]` named access. **Column names are part of the API.** New columns may be appended in any order; existing consumers reference columns by name, not position, so adding a column no longer breaks anything silently.

| DataFrame | Built in | Columns | Consumers |
|---|---|---|---|
| `self.tecdata` | `_storeTectonics` (inputparser.py) | `start, end, tMap, zMap, hMap` | tectonics.py |
| `self.raindata` | `_defineRain` (inputparser.py) | `start, rUni, rzA, rzB, rMap, rKey` | unstructuredmesh.py |
| `self.evapdata` | `_defineEvap` (inputparser.py) | `start, eUni, eMap, eKey` | unstructuredmesh.py, flowplex.py |
| `self.sedfacdata` | `_defineErofactor` (inputparser.py) | `start, sUni, sMap, sKey` | unstructuredmesh.py |
| `self.tedata` | `_getTe` (inputparser.py) | `start, tUni, tMap, tKey` | addprocess.py |

`evapdata` is parsed from the same `climate:` YAML block as `raindata` (per-row `evap_uniform`/`evap_map` are opt-in extensions to each climate event); it is `None` if no row declares evap. The lake-evap budget is computed in `flowplex._potentialLakeEvap` and applied at `step==0` only inside `_distributeDownstream` so spillover-cascade iterations do NOT debit the same pit twice. Both hooks (channel and lake) accumulate into `self.evapLoss` (m³, running total used by water-balance regression tests).

The previous `iloc[nb, k]` positional pattern was replaced in 2026-06 (30 sites across the 3 consumer files). **Do NOT re-introduce `iloc[nb, k]` on these DataFrames** — it makes the column-order a load-bearing API, and any append silently breaks every consumer.

Out-of-scope iloc uses (kept intact): `pitfilling.py` uses `df["col"].iloc[k]` on a local pit-id DataFrame (not one of the four forcing DataFrames); `inputparser.py` uses `seadata[1].iloc[0]` and `icedata[N].iloc[0]` for first/last-row access on CSV-loaded series (also not forcing DataFrames). Different shape, different objects.

## Debugging an elevation spike (`tools/zprobe.py`)
A recurring failure class in goSPL is a **km-scale elevation needle at one cell** that is only noticed several output steps later, by which time the output alone cannot say which stage produced it (`tout` is usually several `dt`, so one output hides multiple steps of every process). `gospl/tools/zprobe.py` exists to answer that in **one run**: `probeZ(model, tag)` prints the global elevation extremes at a named point in the step and how far the maximum moved since the previous probe (tagging the culprit `<== RAISED THE GLOBAL MAXIMUM`), and `reportVolume(model, tag, thickness)` prints a deposit increment's volume **and peak thickness** — which is what separates "this stage moved a lot of sediment" from "this stage moved sediment badly".

**Verbose-only** (`-v`); grep the log for `zprobe`, `mass rescale`, `fDep at the`.

**Why these hook points.** Most deposition is bounded by a physical envelope — a pit fills toward its spill rim (`_bottomUpDelta` clamps at `lFill`), marine deposition cannot exceed the clinoform surface (`maxDepQs = (clinoH − hl)·larea`, and `clinoH ≤ sealevel`), pit infill is capped at the pre-deposition global max. A **short list of sites is deliberately unbounded**, because they exist to conserve mass when the routing has nowhere else to put the sediment, and those are what the probes cover:
- `sedplex._distributeSediment` closed-sink deposit (`_closedDepo`) — active whenever `_domainHasOutlet` is False, i.e. **always on a global sphere**;
- `seaplex._distOcean` force-deposit at terminal sinks, and the residual drained at cascade exit (both bypass `marVol`);
- plus the two amplifiers/conditioners that turn a normal deposit into a needle: the `_diffuseOcean` global mass rescale (the one place a marine deposit can be multiplied up) and the `fDep` 0.99-cap population in `SPL` (see below);
- plus how much the `_updateSinks` `gmax` guard had to clip — which also reveals when that guard has been rendered useless by an earlier unbounded stage having already raised `gmax`.

**NOT instrumented** (add a probe when you need it): ice/glacial till, soil production and creep, groundwater/duricrust, tectonics and horizontal advection, and the `nlSPL`/`soilSPL` erosion flavours. The facility is intentionally partial; keep it honest by extending the module docstring's scope list when you add a hook.

**MPI**: every call is collective (`Vec.max`/`Vec.min`, `allreduce`) and is gated ONLY on `self.verbose`, a config scalar identical on every rank. Adding a probe behind a rank-local condition is the #1 deadlock class — see `## MPI contract`.

**The `fDep` cap is the first thing to check on a deposition spike.** With `spl: G` non-zero, `fDep = G·larea/PA` is capped at 0.99 to keep the coupled `(I−Wᵀ)Q` + `(1−fDep)h` block non-singular. It reaches the cap where `PA ≤ G·larea/0.99`, which expressed as a multiple of a cell's OWN runoff is **`G/(0.99·rain)` — the cell area cancels, so the criterion is mesh-resolution-independent**. At `G=2` with `rain=2` m/yr that ratio is 1.01, so the entire headwater network pins at the cap with `1−fDep = 0.01`; the solve degrades step after step (a real case: `Solve SPL accounting for sediment deposition` climbing 24→72 s and peaking on the step that produced a 25 km needle at a river mouth) and a trunk-river load can land on a handful of cells. At `G=0.5` the ratio is 0.25, nothing saturates, the solve time plateaus and the needle does not form. **Rule of thumb: keep `G ≲ rain/2` (rain in m/yr).** The probe prints the saturated count every step.
