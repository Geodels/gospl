# AGENTS.md

Last reviewed 2026-06-08 against `v2026.06.08`. Read this at the start of every session. Update it when an invariant here changes. See `REFACTOR_AUDIT.md` for the long-form rationale behind each rule.

## What goSPL does
goSPL is a parallel landscape-evolution model that integrates the stream-power law (river incision), linear and non-linear hillslope diffusion, marine sediment transport, glacial accumulation, flexural isostasy, and horizontal/vertical tectonics on an unstructured Voronoi/Delaunay finite-volume mesh. The mesh is either a 2D flat plane (`self.flatModel == True`) or a global sphere; partitioning, halo exchange, and all linear/non-linear solves run on PETSc DMPlex via petsc4py. Time integration is an explicit outer Euler loop in `Model.runProcesses` with implicit KSP/SNES/TS inner solves for diffusion, flow accumulation, and sediment routing.

## The numpy ↔ PETSc boundary
Every state field exists in two parallel representations.

**Numpy land** (raw arrays indexed by local node ID, dimensional, no halo): `self.lcoords` (m), `self.mCoords` (m), `self.larea` (m²), `self.rainVal` (m/yr), `self.upsub` (m/yr), `self.stratH/stratZ/phiS/stratK`, `self.pitParams`, `self.pitIDs`, `self.lFill`, `self.localFlex`, plus any `vec.getArray().copy()` view.

**PETSc land** (parallel Vec with halo, mutated via `setArray`/`getArray`/`localToGlobal`): `self.hLocal`/`self.hGlobal`, `self.cumED`/`self.cumEDLocal`, `self.FAL`/`self.FAG`, `self.fillFAL`, `self.Eb`/`self.EbLocal`, `self.bL`/`self.bG`, `self.areaLocal`/`self.areaGlobal`, `self.iceHL`/`self.iceFAL`/`self.iceFAG`/`self.iceMeltL`/`self.iceFlex`, `self.Lsoil`/`self.Gsoil`, `self.lHbed`/`self.gHbed`, `self.vSed`/`self.vSedLocal`, `self.fiso`.

Both sides hold physical units; the boundary is about **who owns halo synchronisation**, not units. Cross only via `self.dm.localToGlobal(local, global)`, `self.dm.globalToLocal(global, local)`, `vec.getArray()`, `vec.setArray(arr)`. After mutating a `*Local` array view, you MUST `localToGlobal` before the next collective solve, or ranks see stale halos.

## MPI contract
**Collective** (every rank must call, in the same order): `self.dm.localToGlobal`, `self.dm.globalToLocal`, `MPI.COMM_WORLD.Allreduce/Bcast/bcast/Allgatherv/Reduce/Barrier`, `ksp.solve`, `snes.solve`, `ts.solve`, `vec.sum/max/min`, `vec.assemblyBegin/End`, `mat.assemblyBegin/End`, `vec.duplicate/destroy`, `mat.destroy`, `dm.distribute`. **Rank-local**: `vec.getArray()`, `vec.setArray()`, all numpy ops, anything inside `if MPIrank == 0:`.

**PETSc initialisation happens exactly once**, in `gospl/__init__.py` (`petsc4py.init(sys.argv)`, line 25). Python guarantees the package `__init__` runs before any submodule, so module-level code in submodules (e.g. `MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()` at import time) can rely on PETSc being live. **Do NOT re-introduce `petsc4py.init` in any submodule** — until 2026-06 every submodule called it at import time (15 sites); the call is idempotent so duplicates were harmless but obscured where state was created. Submodules still `import petsc4py` to access `petsc4py.PETSc.X` symbols; that's a separate concern from `init()`.

`MPIcomm` is defined locally in 5 active files. 9 dead-code assignments were removed 2026-06. Active sites already follow the rule below.

| File | `MPIcomm =` | Used for |
|---|---|---|
| `flow/flowplex.py` | `petsc4py.PETSc.COMM_WORLD` | `Mat().create(comm=MPIcomm)` |
| `sed/seaplex.py` | `petsc4py.PETSc.COMM_WORLD` | `Mat/Vec().createNest(comm=MPIcomm)` |
| `eroder/SPL.py` | `petsc4py.PETSc.COMM_WORLD` | `Mat/Vec().createNest(comm=MPIcomm)` |
| `mesher/unstructuredmesh.py` | `MPI.COMM_WORLD` | `bcast`, `Barrier` |
| `tools/outmesh.py` | `MPI.COMM_WORLD` | `bcast`, `gather`, `Barrier` |

Rule: use `MPI.COMM_WORLD` for raw collectives (Allreduce/bcast/Allgatherv); use `petsc4py.PETSc.COMM_WORLD` only when creating PETSc objects (`KSP().create(comm=...)`, `Mat().create(comm=...)`). They wrap the same handle but go through different paths inside PETSc.

## KSP / SNES / TS lifecycle contract
PETSc solvers in goSPL follow two intentional patterns. **Use the right one for new code.**

### CACHED — hot-path solvers (8 sites)
Solvers that run **every timestep** are created lazily on first use, stored as `self._X`, and reused across the entire simulation. Avoids the ~5-10 ms create+destroy churn per call, which compounds over thousands of timesteps. Each site's inline docstring confirms the rationale.

| File | Method | Cached attribute | Solver |
|---|---|---|---|
| `flow/flowplex.py` | `_solve_KSP` | `self._ksp_main` | KSP (main: richardson+bjacobi) |
| `flow/flowplex.py` | `_solve_KSP2` | `self._ksp_fallback` | KSP (fallback: fgmres+asm) |
| `eroder/nlSPL.py` | `_solveNL_ed` | `self._snes_ed` | SNES (transport-limited non-linear SPL) |
| `eroder/nlSPL.py` | `_solveNL` | `self._snes_nl` | SNES (detachment-limited non-linear SPL) |
| `eroder/soilSPL.py` | `_solveSoil` | `self._snes_soil` | SNES (soil-aware SPL) |
| `eroder/soilSPL.py` | `diffuseSoil` | `self._ts_soil` | TS (rosw soil diffusion) |
| `sed/hillslope.py` | `_hillSlopeNL` | `self._snes_hill` | SNES (non-linear hillslope) |
| `sed/hillslope.py` | `_diffuseImplicit` | `self._ts_marine` | TS (rosw marine + lake diffusion) |

**Lifecycle**:
```python
if self._snes_x is None:
    snes = petsc4py.PETSc.SNES().create(comm=petsc4py.PETSc.COMM_WORLD)
    # ...configure...
    self._snes_x = snes
snes = self._snes_x
# ...solve...
# Do NOT call snes.destroy() — destroy_DMPlex handles it at simulation end.
```

**CRITICAL**: any new cached solver MUST be added to the `destroy_DMPlex` loop in `mesher/unstructuredmesh.py:738-799`. The loop iterates over a hardcoded list of attribute names; forgetting to add yours leaks the PETSc object at simulation end. Same applies to cached helper Vecs (`self._snes_X_f`, `self._snes_X_x`, `self._snes_X_J`, etc.).

### AD-HOC — nested-matrix fieldsplit solves (2 sites)
Solvers that build a `Mat().createNest(...)` whose sub-matrices change every call AND configure `pc.setType("fieldsplit")` with IS sets derived from the nested-mat structure. The fieldsplit PC's IS configuration is tied to the specific sysMat instance, so re-using a cached KSP via `setOperators(new_sysMat)` does NOT automatically re-derive the splits. Caching is theoretically possible but requires careful experimentation with `pc.reset()` and the nested-mat lifecycle.

| File | Method | Condition |
|---|---|---|
| `eroder/SPL.py` | `_coupledEDSystem` | `self.fDepa != 0` (transport-limited branch with non-zero `G`) |
| `sed/seaplex.py` | `_depMarineSystem` | `not flatModel AND self.Gmar > 0`, AND only inside the second `_distOcean` pass |

**Lifecycle**: create at the top of the method, configure, solve, then explicitly destroy everything (KSP, sub-KSPs, sub-ISes, PC, sysMat, RHS/solution vectors) at the end. See `SPL.py:191-243` for the canonical pattern. Both sites are COLD-path (conditional, not every step), so the cumulative create+destroy overhead is small.

If a future contributor wants to convert one of these to CACHED, the obstacle is the fieldsplit-PC + nested-mat IS lifecycle, not the KSP object itself. Do it on a focused branch with full regression run.

### Common rules
- All comm arguments use `petsc4py.PETSc.COMM_WORLD` (never `MPI.COMM_WORLD`). Matches the MPI contract above.
- Positional (`KSP().create(PETSc.COMM_WORLD)`) vs keyword (`SNES().create(comm=PETSc.COMM_WORLD)`) is purely cosmetic; both work identically.

## The Model god-class
`gospl/model.py:93-110` declares `Model` as multi-inheritance of 16 mixins. `model.py:126-198` calls each parent's `__init__` **by name, not via `super()`**. The init order is load-bearing and differs from the MRO declaration order — adding `super().__init__()` will break the chain.

| # | Init call (model.py line) | Assumes already populated | Allocates / sets |
|---|---|---|---|
| 1 | `_ReadYaml(filename)` :136 | — | `self.input`, every YAML attr, `self.tNow` |
| 2 | `_STRAMesh()` :139 | ReadYaml (`strataFile`) | `stratH/stratZ/phiS/stratK = None` |
| 3 | `_VoroBuild()` :142 | — | Voronoi cache attrs reset |
| 4 | `_UnstMesh()` :145 | ReadYaml + STRAMesh | `dm`, `hLocal`, `hGlobal`, `locIDs`, `glbIDs`, `lpoints`, `mpoints`, `lcoords`, `mCoords`, `larea`, `FVmesh_ngbID`, `lgmap_*`, `idBorders/idLBounds/ghostIDs`, `bL/bG`, calls `readStratLayers` |
| 5 | `_WriteMesh()` :148 | UnstMesh | `step`, `outputDir`, `upsG/upsL` |
| 6 | `_FAMesh()` :151 | UnstMesh | `iMat`, `fillFAL`, `FAG`, `FAL`, `rtol`; defines `_matrix_build`, `_matrix_build_diag`, `_solve_KSP*` for every downstream class |
| 7 | `_IceMesh()` :154 | FAMesh | `iceHL/iceFAG/iceFAL/iceMeltL/iceFlex` (only if `iceOn`) |
| 8 | `_SPL()` :157 | FAMesh | `hOld`, `hOldLocal`, `hOldFlex`, `Eb`, `EbLocal`, `stepED`, `newH` |
| 9 | `_nlSPL()` :160 | SPL | `snes_rtol/atol/maxit`, lazy `_snes_ed/_snes_nl` |
| 10 | `_soilSPL()` :163 | nlSPL | `Gsoil`, `Lsoil`, `lHbed`, `gHbed`, `prodSoil`, `soil_transition` |
| 11 | `_PITFill()` :166 | UnstMesh | `borders`, `outEdges` |
| 12 | `_SEDMesh()` :169 | FAMesh | **`tmp`, `tmpL`, `tmp1`**, `Qs`, `QsL`, `nQs`, `vSed`, `vSedLocal`, `maxnb` |
| 13 | `_hillSLP()` :172 | SEDMesh | **`h`, `hl`, `dh`**, `mat`, `Dlimit/dexp/minDiff`, lazy `_snes_hill`, `_ts_marine` |
| 14 | `_SEAMesh()` :175 | FAMesh | `zMat` |
| 15 | `_GridProcess()` :178 | ReadYaml + UnstMesh | `localFlex`, regular grid (FD/FFT) or DH grid (global) |
| 16 | `_UnstMesh.applyForces(self)` :181 | everything above | `rainVal`, `upsub`, `sealevel` |
| 17 | `_Tectonics()` :184 | UnstMesh (`hGlobal`) | `fiso`, `tecNb=-1` |

Permuting any line silently uses an unallocated attribute or overwrites one another mixin has already set.

## Scratch vector contract (CRITICAL)
The following PETSc Vecs are **scratch**. Any kernel may overwrite them at any time. They are **NOT persistent state** — they exist only so collective allocations are paid once at init.

| Vec | Allocated at | Local/Global |
|---|---|---|
| `self.tmp` | sedplex.py:30 | global |
| `self.tmpL` | sedplex.py:31 | local |
| `self.tmp1` | sedplex.py:32 | global |
| `self.Qs`, `self.QsL`, `self.nQs` | sedplex.py:33-35 | g, l, l |
| `self.h` | hillslope.py:39 | global |
| `self.hl` | hillslope.py:40 | local |
| `self.dh` | hillslope.py:41 | global |
| `self.newH` | SPL.py:41 | global (also held by SNES) |
| `self.stepED` | SPL.py:39 | global |
| `self.upsG`, `self.upsL` | outmesh.py:64-65 | global, local |

**Reading `self.h` for elevation gives you whatever the last hillslope step wrote.** Use `self.hLocal/hGlobal` for elevation, `self.cumED/cumEDLocal` for cumulative ED, `self.FAL/FAG` for flow accumulation, `self.Eb/EbLocal` for erosion rate, `self.vSed/vSedLocal` for sediment volume.

Rule: any new method that uses a scratch Vec MUST document which ones in its docstring and MUST leave them in a defined state on exit (typically `set(0.0)` or a freshly-written array, never mid-computation).

## The rcvID / rcvIDi convention (CRITICAL)
`flowAccumulation` (`flowplex.py:421-426`) snapshots six arrays **after the first `_buildFlowDirection` on the unfilled topography, before pit filling and before downstream-routing rebuilds the flow matrix**:
```
self.wghtVali = self.wghtVal.copy()
self.rcvIDi   = self.rcvID.copy()
self.distRcvi = self.distRcv.copy()
self.fMati    = self.fMat.copy()
self.lsinki   = self.lsink.copy()
```
The `i` suffix means **initial (pre-fill)**. All SPL kernels (`eroder/SPL.py`, `eroder/nlSPL.py`, `eroder/soilSPL.py`) and `sedplex._getSedFlux` MUST use the `i`-suffix versions. The live `self.rcvID/wghtVal/distRcv/fMat/lsink` are valid **only** inside `flowplex._distributeDownstream` and `sedplex._moveDownstream`, where they are rebuilt against the current filled/sediment-filled topography. Outside those two functions their state is undefined.

## Eb / EbLocal sign convention (unified, thickness-rate)
**Both `self.Eb` (global) and `self.EbLocal` (local) are in the thickness-rate convention: positive for deposition, negative for incision.** Same sign as `cumED`, same sign as the on-disk `EDrate` output field, same sign as the restart loader (`outmesh.py:439-440`). Unified 2026-06.

What each field contains by end-of-step:
- **`self.Eb`** — river-only thickness rate from the most recent SPL flavour (`SPL.py:_getEroDepRate` / `nlSPL.py:_getEroDepRateNL` / `soilSPL.py:_getEroDepRateSoil`). Not re-synced after marine/hillslope contributions, so it reflects ONLY the river step.
- **`self.EbLocal`** — net thickness rate including all axpy contributions from later kernels (`seaplex.py:486`, `hillslope.py:297`, `soilSPL.py:549`). This is what `outmesh.py:272` writes to disk as `EDrate`.

Same convention, different content. `self.Eb` and `self.EbLocal` are NOT local/global views of the same field — the wrapper at the end of each `erodepSPL*` overwrites `EbLocal` with `add_rate = tmp/dt = Eb`, then downstream kernels mutate `EbLocal` only.

Thickness conversion (in case of future refactors):
- Inside `_getEroDepRate*`: `tmp = stepED - hOld` is the elevation change → divide by `dt` directly to get the thickness-rate `Eb` (no inversion).
- Inside `erodepSPL*` wrapper: `tmp = Eb * dt` is the signed thickness change (negative at incising, positive at depositing cells).
- `cumED.axpy(1.0, tmp)` → cumED in thickness convention.
- `hGlobal.axpy(1.0, tmp)` → at incision, h drops correctly.

One quirk worth knowing: **`sedplex._getSedFlux` (sedplex.py:63) negates `self.Eb` before the upstream-integration solve**, because that solve needs an erosion-positive source for `vSed` (m³/yr) to accumulate as positive downstream flux. The `stratNb > 0` branch uses `self.thCoarse` which is already erosion-positive (from `stratplex.erodeStrat`), so no negation there.

## The `_extra*` methods are mandatory continuations
NOT optional parsers. Each sets attributes required by other modules. Never delete or rename without following the full call chain in `inputparser._readDomain/_readTime/_readHillslope/_readFlex/_readOrography/_readIce`.

- `_extraDomain` (inputparser.py:189) → `seaDepo`, `overlap`, `dataFile`, `nodep`, `strataFile`; calls `_extraDomain2`.
- `_extraDomain2` (:229) → `advscheme`, `radius`, `gravity`.
- `_extraHillslope` (:475) → `nlK`, `clinSlp`, `tsStep`, `Gmar`, `offshore`, `nl_pit_volume/depth/K/inlet_bias`.
- `_extraFlex` (:1430) → `nu`, `flex_res_deg`, `flex_bcN/S/E/W`.
- `_extraOrography` (:1517) → `oro_cw`, `oro_conv_time`, `oro_fall_time`, `oro_precip_*`, `rainfall_frequency`.
- `_extraIce` (:1621) → `iceT`, `elaH`, `iceH` interpolators.

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

## Magic numbers
All sentinels and threshold values are defined in `gospl/tools/constants.py` and imported by name. **These literal values MUST NOT be reintroduced inline.** If you add a new constant, add an entry to both this section AND `gospl/tools/constants.py` so the two stay in sync.

| Constant | Value | Role |
|---|---|---|
| `MISSING_DATA_SENTINEL` | `-1.0e8` | Pre-fill before `Allreduce(MAX)` + boundary marker before `fitedges`. |
| `MISSING_LARGE_SENTINEL` | `-1.0e10` | Same as above but for fields whose magnitude can exceed 1e8 (cumED, flexure, soil thickness — used in `tectonics._advectPlates`). Two orders below `MISSING_DATA_SENTINEL` so they never collide. |
| `DISCHARGE_FLOOR` | `1.0e-8` | Minimum FA/discharge/sedLoad value written to HDF5 outputs (avoids `-inf` in log10 viz). |
| `DEPOSIT_FLOOR` | `1.0e-3` (1 mm) | Drop sub-mm marine sediment as numerical noise. |
| `BEDROCK_EXPOSED` | `1.0e-1` (10 cm) | Soil-vs-bedrock threshold in soil-aware SPL. |
| `BEDROCK_SENTINEL` | `1.0e6` | Infinite-bedrock layer-0 sentinel thickness; offset cancels in cumsum arithmetic. |
| `BOUNDARY_FLOW_SENTINEL` | `-1.0e6` | Mfd no-data marker passed to `mfdreceivers`/`mfdrcvrs`; **distinct from `MISSING_DATA_SENTINEL`** — kept two orders apart so a mix-up is obvious in diagnostics. |
| `GRAPH_OUTLIER_CAP` | `1.0e7` | Upper-bound clamp in `pitfilling._performFilling`; entries above this are rewritten to `MISSING_DATA_SENTINEL`. |

**Same value, different role — DO NOT replace these with the listed constants:**
- `1.0e-8` at `flow/iceplex.py:208` (ice-presence threshold), `sed/stratplex.py:219` (thickness numerical-noise floor).
- `1.0e-3` at `flow/flowplex.py:363` (water-routing convergence), `sed/sedplex.py:139` (sediment-routing convergence), `flow/pitfilling.py:568,622` (minh epsilon nudges), `sed/seaplex.py:507,509` (clinoH 1mm offset).
- `1.0e-1` at `flow/iceplex.py:231` (minimum ice thickness for visualization).

Each of these is marked with a permanent `# TODO-REFACTOR: value matches X but distinct role; do not replace` comment so future readers know the coincidence is intentional.

## High-risk modules (do not edit without full regression run)
- **`mesher/unstructuredmesh.py`** — owns `dm`, every shared mesh attribute, forcing dispatch, and `destroy_DMPlex` (lines 738-799) which names every Vec/Mat by hand. Adding a new persistent Vec elsewhere requires adding it to that destroy list or it leaks.
- **`flow/flowplex.py`** — owns `_solve_KSP`, `_solve_KSP2`, `_matrix_build`, `_matrix_build_diag`, `_buildFlowDirection`. Consumed by every downstream module.
- **`tools/inputparser.py`** — owns all parameter parsing; forcing DataFrame column order is API; the `_extra*` chain is mandatory.

## Known bugs (fix before refactoring)
- _(none currently open)_

## Fixed (regression-guarded)
- **`rUni`/`sUni` key/column mismatch** — fixed 2026-06 at `inputparser.py:915` (dict key `'rUni'` → `'sUni'` in the uniform branch of `_defineErofactor`). Before the fix, a YAML with a uniform `sedfactor` event landed `NaN` in `sedfacdata['sUni']`, then crashed `unstructuredmesh._updateEroFactor` on `np.load(None)`. Guarded by `tests/test_regression.py::test_uniform_sedfactor_populates_sUni`.
- **Marine sediment leak at terminal ocean sinks** — fixed 2026-06 at `seaplex._distOcean`. `_matOcean` constructs `dMat1` with a zero column at every terminal sink (cells whose every flow direction is self-referential, typically deep-ocean basin floors); the first `dMat1.mult` therefore annihilated any `sedflux` that landed directly at a sink. `vSed` upstream-integration naturally accumulates large flux at these sinks, so this was the dominant mass-loss path. Diagnostic: 9.748e12 m³ marine sediment supply per step → 4.107e12 m³ landing at sinks (42% of input) → 0 m³ deposited via the routing loop → silently dropped. The fix adds a pre-drain step at the top of `_distOcean` that deposits any sinkVol at terminal sinks directly into `vdep` before the first `dMat1.mult`. The in-loop force-deposit and exit-residual drain (added alongside) handle any sediment that subsequently arrives at sinks during multi-step routing. Guarded by `tests/test_regression.py::test_mass_conservation`, which runs on the `minimal_model` fixture — a **global-sphere mesh** (`flatModel=False`, `len(idBorders)==0`), so the domain is closed by construction with no boundary outflux possible; the test asserts `|dV_surface|/activity < 1e-4` AND `|dV_cumED|/activity < 1e-4`, i.e. that deposited volume equals eroded volume to within ~6e-5 relative (the floor-effect budget from `DEPOSIT_FLOOR=1e-3` and pit-routing residue, both documented in the test's tolerance rationale). The corresponding `is_sink_local` mask is computed in `_matOcean` once per call.
- **`Eb` vs `EbLocal` sign-convention divergence** — fixed 2026-06 by adopting the thickness-rate convention (positive deposition, negative incision) for both fields. Changes: drop the leading minus in `E = -tmp.getArray()` inside `_getEroDepRate` (SPL.py), `_getEroDepRateNL` (nlSPL.py), `_getEroDepRateSoil` (soilSPL.py); drop the minus in `tmp.setArray(-Eb * dt)` inside the three matching `erodepSPL*` wrappers; negate `Eb` inside `sedplex._getSedFlux` (sedplex.py:63) where the upstream-integration solve still wants an erosion-positive source. Before the fix, `Eb` was incision-positive while `EbLocal` (after the wrapper's `add_rate = tmp/dt` overwrite) was deposition-positive, so the two fields disagreed in sign by end-of-step; the on-disk `EDrate` field was always deposition-positive (because it's written from `EbLocal`), and the restart loader at `outmesh.py:439-440` restored `Eb` in deposition-positive convention from `EDrate` — so fresh-run state and restart-state were inconsistent before this fix. Guarded by `tests/test_regression.py::test_erosion_sign_conventions` (which asserts `EbLocal <= 0` at incising nodes; `Eb` follows the same convention now).

## Intentional surprises (do NOT "fix")
- **N/S boundary swap in gFlex**. `addprocess.py:221-222` deliberately assigns `simflex.BC_S = self.flex_bcN` and `simflex.BC_N = self.flex_bcS` to compensate for a coordinate-convention mismatch with gFlex. Verified intentional in project memory.
- **`gid` argument on `mfdreceivers` / `mfdrcvrs`**. The Fortran kernels take an extra `gid(nb)` argument (per-local-node global vertex ID, from `self.gid` set in `mesher/unstructuredmesh.py` immediately after `getLGMap`). Inside the kernels, the stored slope is perturbed by `val * (1 - 1.0e-15 * gid(n))` before the quicksort tie-break. This makes EXACT slope ties resolve deterministically across MPI decompositions (lower global ID wins). The perturbation is well below KSP-solver precision and does not affect physical results; near-ties driven by KSP floating-point noise are a separate problem and are NOT fixed by this pattern. See the relaxation comment on `test_parallel_correctness` in `tests/test_regression.py` and `fortran/functions.F90:mfdreceivers` for the full rationale. **Do NOT remove the `gid` argument** — without it, `mfdreceivers` ordering depends on local iteration order and the parallel-correctness test diverges by ~10x from current baseline.
- **Platform-dependent KSP-noise floor in `test_parallel_correctness`**. With the `gid` exact-tie-break in place, the remaining sum-FA drift between n=1 and n=2 is dominated by KSP-solver floating-point non-determinism at partition boundaries — and the magnitude differs by platform. Measured on `tests/fixtures/minimal.yml`: macOS-14 (arm64, conda-forge OpenMPI) ≈ **0.3%**, ubuntu-latest (x86_64, conda-forge MPICH) ≈ **1.7%**. Same fixture, same goSPL, same Python — the spread is driven by MPI implementation, BLAS variant, and FMA availability at the platform level. The test tolerance for `rel_sum_fa` is set at **5e-2** (5%) to clear both with headroom; a real routing regression would shift it toward 0.5+, not 0.05. If a future contributor sees these numbers and wants to tighten the tolerance, the correct path is to make KSP-precision halo state bitwise-identical across decompositions (hard PETSc work), NOT to lower the bound and hope.

## CI contract
- Workflows: `.github/workflows/tests-pr.yml` (fast tier), `tests-slow.yml` (slow regression + analytical benchmarks), `conda-build.yml` (package build/publish on tag push).
- Matrix in every workflow: `ubuntu-latest + macos-14 × Python 3.11 + 3.12`.
- Conda setup: `miniforge-version: latest` (libmamba is the default solver in modern Miniforge3 — no `use-mamba: true` needed). `channels: conda-forge` with `channel-priority: strict` and `conda-remove-defaults: true`. **No `auto-update-conda: true`** (latest miniforge already gives us latest conda; the extra update step costs minutes for nothing).
- Conda env caching: `actions/cache@v4` step keyed on `runner.os + runner.arch + python + hashFiles('environment.yml')` runs immediately before `setup-miniconda` in `tests-pr.yml` and `tests-slow.yml`. First run after any `environment.yml` change is a full ~5-10 min build (or ~60 min on ubuntu-latest if conda-forge's CDN is throttling); subsequent runs restore the cached env dir and skip most of the install. Bump the `-v<n>-` segment in the cache key to manually flush.
- **Cache path must NOT use `${{ env.CONDA }}`.** On the hosted Ubuntu runner, `env.CONDA` defaults to `/usr/share/miniconda` (the pre-installed system conda), but `setup-miniconda` with `miniforge-version: latest` installs Miniforge3 to `/home/runner/miniconda3`. Caching `${{ env.CONDA }}/envs/gospl` therefore caches an empty directory and the env build runs from scratch every time. The cache step uses `${{ runner.os == 'Linux' && '/home/runner/miniconda3/envs/gospl' || '/Users/runner/miniconda3/envs/gospl' }}` to hit the right location on each OS. macOS-14 happens to coincide with the pre-installed default, but the conditional is still the load-bearing piece for Linux.
- Hard timeout: 240 minutes per `tests-slow.yml` job (benchmarks dominate wall-clock). Fast tier uses the default (no explicit cap).
- goSPL is installed with `pip install --no-deps --no-build-isolation -e . -v`. `--no-deps` because the conda env already supplies every runtime dependency; `--no-build-isolation` because pyshtools rebuilds from source under pip's isolated env and fails on CI runners.
- `environment.yml` conventions — do NOT revert any of these (each one cost a CI iteration to discover):
  - **Python pinned as `python>=3.11,<3.13`, not `=3.11`.** A strict pin conflicts with the workflow's `python-version` matrix override on py3.12 cells, producing `Could not solve for environment specs: python =3.11 vs =3.12` and failing the entire job before any tests run.
  - **`channels:` lists `conda-forge` only — no `defaults`.** With `defaults` in the list, conda pulls older variants of `build` and other packages from `repo.anaconda.com` that don't support newer Pythons, even with `channel-priority: strict` set workflow-side. Also produces libmamba warnings about Anaconda's commercial-channel ToS.
  - **No `pip: - git+...gospl.git` section.** goSPL is installed exclusively by the workflow's `pip install -e .` step. A previous in-env duplicate caused ~30-min CI env builds (full goSPL clone + Fortran rebuild inside the conda phase, then the same build a second time in the workflow's pip step). Local-dev impact: `mamba env create -f environment.yml` now provides deps only; follow with `pip install -e .` as the standard editable-install pattern.
  - **No `build` package in dependencies.** conda-forge has it only as `python-build`, and `meson-python`'s build backend doesn't `import build` anyway. With `pip install --no-build-isolation`, `pyproject.toml`'s `[build-system].requires = ["build"]` is also not consulted — so it's pure dead weight. Adding `build` literally would only resolve via the `defaults` channel (which we don't allow — see above), so a future contributor wondering "why isn't `build` in here?" should stop there.
- **On every PR + push to `master`/`release-candidate`** (`tests-pr.yml`):
    `pytest tests/`
  No marker filter — runs all 13 regression tests (parser-tier + Model-tier). Wall-clock ~5–10 min per cell after the env cache populates.
- **On nightly cron (04:00 UTC), tag push (`v*`), or workflow_dispatch** (`tests-slow.yml`, with `if: github.event_name != 'push'`):
    `pytest tests/ -m 'slow and not benchmark'`
  Gated OFF push events because the same 9 tests already run in `tests-pr.yml` on every push. The non-push triggers (cron / tag / dispatch) keep the slow tier exercised against upstream-package drift.
- **On push to `master`/`release-candidate` only** (gated step inside `tests-slow.yml`):
    `pytest benchmarks/ -m 'benchmark'`
  Skipped on cron / tag / dispatch so release-tagging is not blocked by hours-long benchmark wall-clock.
- Benchmark artifacts (junit XML, Markdown reports, PDF figures under `results/`) uploaded to GitHub on every slow+benchmark run, including failures (`if: always()`). The artefacts land at `<repo_root>/results/<benchmark_name>/` via the teardown copy in `benchmarks/conftest.py`, not directly from `tmp_path`.
- Benchmark failures block merge (`continue-on-error: false`).
- Concurrency: `tests-pr.yml` cancels in-flight runs on the same ref (`cancel-in-progress: true`); `tests-slow.yml` does not (`cancel-in-progress: false`) — slow runs are expensive enough that killing them mid-flight is more wasteful than letting them finish.

## Analytical benchmark suite
Tests in `benchmarks/` validate goSPL physics against exact analytical solutions. They require `scipy` and `matplotlib` and are silently skipped via `pytest.importorskip` in environments without these packages. `matplotlib` is in `environment.yml` for CI; it is NOT a goSPL runtime dependency (intentional — `pyproject.toml` does not list it).

The `benchmark` and `slow` marker names are registered in `pyproject.toml` under `[tool.pytest.ini_options].markers` (authoritative) and mirrored in `tests/conftest.py::pytest_configure` for the existing tests/ tree. Don't expect new markers to work without registering them in pyproject — pytest does not walk sideways into `tests/conftest.py` when invoked from `benchmarks/`.

Each benchmark fixture (`spl_tmp_path`, `hillslope_tmp_path`, `knickpoint_tmp_path` in `benchmarks/conftest.py`) copies the benchmark's input files into a per-test tmp directory and chdirs there. goSPL intermediate output (`sim_output/`, `sims_outputs/`) stays inside `tmp_path` and is reaped by pytest's cleanup. Report artefacts (`*.pdf`, `*.md`, hillslope `*.png`) are copied out on teardown into `<repo_root>/results/<benchmark_name>/` so they survive pytest cleanup AND are picked up by the GitHub Actions `upload-artifact@v4` step in `tests-slow.yml`. `results/` is in `.gitignore` — never commit it.

The SPL fixture's teardown glob is `spl_benchmark*.pdf` / `spl_benchmark*.md` (wildcarded) because `test_spl.py` writes one PDF + one Markdown per benchmark case (`spl_benchmark_100.pdf`, `spl_benchmark_150.pdf`, `spl_benchmark_200.pdf`, plus the combined `spl_benchmark_all.md`). Don't tighten the glob — the wildcard is load-bearing.

Mesh `.npz` files under each benchmark's `boundary_condition[s]/` subfolder are working-tree artefacts (currently untracked in git — see `git status`). Naming convention diverges by benchmark: `spl/boundary_conditions/` and `hillslope/boundary_conditions/` are **plural**, `knickpoint/boundary_condition/` is **singular**. Match the existing folder when adding a new mesh; the YAML configs reference these paths directly. `knickpoint/boundary_condition/drop_baselevel.npz` is regenerated by the test itself during phase 2 of the knickpoint workflow.

`benchmarks/test_hillslope.py` delegates to `scripts/analysis.py` (imported as `anlys` after `hillslope_tmp_path` prepends the copied folder to `sys.path`). The other two benchmarks inline all logic in their `test_*.py` files. `analysis.runFullBenchmark` does NOT expose an `overall_pass` key; use `analysis.print_summary(results)` for the pass/fail bool.

`benchmarks/test_spl.py` runs **three cases per invocation** (`input100.yml`, `input150.yml`, `input200.yml` under `benchmarks/spl/sims/`). `evaluateSPL` accepts a `skip_basin_test` flag — passed `True` for cases 150 and 200 (basin-geometry test is not meaningful for those rates), reducing their pass criteria from 6/6 to 5/5. The test asserts `overall_pass` per-case; failure on any case fails the whole test. Helpers: `_run_spl_benchmarks` (runs all cases + writes per-case + combined Markdown) and `_assert_spl_benchmarks` (the assertion).

| Test file | Process | Analytical basis | Pass criteria |
|---|---|---|---|
| `test_spl.py` | SPL steady state | Braun & Willett 2013; Perron & Royden 2013 | Case 100: 6/6; cases 150 & 200: 5/5 (basin test skipped) |
| `test_hillslope.py` | Hillslope diffusion | `z=(U/2κ)x(L-x)`; Roering, Kirchner & Dietrich 1999 | All `TOL_*` constants met |
| `test_knickpoint.py` | Knickpoint propagation | `c=K·A^m`; Royden & Perron 2013; Tucker & Whipple 2002 | 4/4 sub-tests |

## Conda Package Validation

### Released packages
| Version | Date | Channel | Install |
|---|---|---|---|
| `v2026.06.08` | 2026-06-08 | `geodels` | `mamba install -c geodels -c conda-forge gospl` |

### Local build and smoke-test procedure (osx-arm64)
Run this sequence from the repository root before pushing a release tag to
`geodels`. It validates the conda recipe, builds the package, and exercises
the full fast test suite against the *installed* package (not the source tree).
```bash
mamba install -n base -y conda-build        # one-time; skip if already present
conda build purge                           # clear stale build cache
conda build conda/ \
  -c conda-forge \
  --override-channels \
  --python 3.11 \
  --variants '{"python": ["3.11"]}' \
  2>&1 | tee build.log
mamba create -n gospl-smoke python=3.11 -c conda-forge -y
mamba install -n gospl-smoke \
   $PKG.conda \
  -c local -c conda-forge -y
mamba install -n gospl-smoke pytest -c conda-forge -y
cd /tmp
mamba run -n gospl-smoke python -c \
  "from gospl.model import Model; print('ok')"
mamba run -n gospl-smoke python -m pytest \
  /path/to/gospl/tests/ \
  -v --tb=short \
  --import-mode=importlib
cd -
mamba env remove -n gospl-smoke -y
```

### Known platform constraints (osx-arm64)
- **vtk → vtk-base**: the full `vtk` package on osx-arm64 pulls in `gtk3` /
  `gdk-pixbuf`, whose post-link script fails inside the conda-build sandbox.
  The recipe uses `vtk-base` (headless subset) instead. goSPL only uses VTK
  for unstructured mesh I/O, not rendering.
- **petsc4py ABI mismatch**: conda-forge does not yet publish a `py311`
  osx-arm64 build of `petsc4py`. The `np2py310` build is installed and works
  for all computation, but causes a segfault during MPI finalization when
  spawned as a subprocess via `mpirun`. `test_parallel_correctness` detects
  this at runtime via `_petsc4py_abi_mismatch()` and skips rather than
  failing. The test runs correctly on linux-64 (which has proper `py311`
  builds) and in `gosplenv` (where the mismatch is masked by the environment
  history).
- **Multi-version render**: pass `--variants '{"python": ["3.11"]}'` to
  conda-build to prevent it rendering the recipe for all Python versions in
  the conda-forge global pinnings file (which includes 3.13, incompatible
  with `numpy=1.26`).

## Milestones

| Date | Tag | Description |
|---|---|---|
| 2026-06-08 | `refactor-baseline-2026-06` | Tier 2 AI-readability refactor complete. AGENTS.md written, 6 regression tests passing, 3 scientific bugs fixed (rUni/sUni, marine sediment leak, Eb sign convention), constants.py, _get_param, named DataFrame access, KSP lifecycle documented, HOW_TO_ADD_FORCING.md and HOW_TO_ADD_OUTPUT.md written. |
| 2026-06-08 | `v2026.06.08` | First release from refactored codebase. Analytical benchmark suite integrated and green on all CI cells (ubuntu-latest + macos-14 × Python 3.11 + 3.12). Published to `geodels` conda channel. |

## Checklist before any commit
1. Did you read this file? If invariants here changed, update them.
2. Did you run the full regression test suite
   (pytest tests/ -v)?
3. If you added a column to a forcing DataFrame, did you use named access (df.at[nb, col]) not iloc?
4. If you used a scratch Vec, did you document which ones in the method's docstring?
5. If you changed a method called by `Model.runProcesses` (`model.py:217-286`), did you check every caller AND the mixin init order (`model.py:126-198`)?
6. If you added a new KSP/SNES/TS solver, did you pick the right lifecycle (CACHED for hot-path solvers, AD-HOC for nested-fieldsplit) AND, if cached, add the attribute to the `destroy_DMPlex` list in `unstructuredmesh.py`?
7. If you added a new forcing type, did you follow `docs/HOW_TO_ADD_FORCING.md` including the `destroy_DMPlex` registration?
8. If you added a new output field, did you follow `docs/HOW_TO_ADD_OUTPUT.md` including the `destroy_DMPlex` registration and XDMF entry?
9. If you added a benchmark test, did you apply `pytest.importorskip` for scipy/matplotlib AND mark `@pytest.mark.benchmark @pytest.mark.slow`?
10. If you added a goSPL `Model` init in a test (regression OR benchmark), is `model.destroy()` in a `try/finally` block per AGENTS.md > KSP/SNES/TS lifecycle contract?
