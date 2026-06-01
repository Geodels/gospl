# AGENTS.md

Last reviewed 2026-05-30 against `release-candidate`. Read this at the start of every session. Update it when an invariant here changes. See `REFACTOR_AUDIT.md` for the long-form rationale behind each rule.

## What goSPL does
goSPL is a parallel landscape-evolution model that integrates the stream-power law (river incision), linear and non-linear hillslope diffusion, marine sediment transport, glacial accumulation, flexural isostasy, and horizontal/vertical tectonics on an unstructured Voronoi/Delaunay finite-volume mesh. The mesh is either a 2D flat plane (`self.flatModel == True`) or a global sphere; partitioning, halo exchange, and all linear/non-linear solves run on PETSc DMPlex via petsc4py. Time integration is an explicit outer Euler loop in `Model.runProcesses` with implicit KSP/SNES/TS inner solves for diffusion, flow accumulation, and sediment routing.

## The numpy ↔ PETSc boundary
Every state field exists in two parallel representations.

**Numpy land** (raw arrays indexed by local node ID, dimensional, no halo): `self.lcoords` (m), `self.mCoords` (m), `self.larea` (m²), `self.rainVal` (m/yr), `self.upsub` (m/yr), `self.stratH/stratZ/phiS/stratK`, `self.pitParams`, `self.pitIDs`, `self.lFill`, `self.localFlex`, plus any `vec.getArray().copy()` view.

**PETSc land** (parallel Vec with halo, mutated via `setArray`/`getArray`/`localToGlobal`): `self.hLocal`/`self.hGlobal`, `self.cumED`/`self.cumEDLocal`, `self.FAL`/`self.FAG`, `self.fillFAL`, `self.Eb`/`self.EbLocal`, `self.bL`/`self.bG`, `self.areaLocal`/`self.areaGlobal`, `self.iceHL`/`self.iceFAL`/`self.iceFAG`/`self.iceMeltL`/`self.iceFlex`, `self.Lsoil`/`self.Gsoil`, `self.lHbed`/`self.gHbed`, `self.vSed`/`self.vSedLocal`, `self.fiso`.

Both sides hold physical units; the boundary is about **who owns halo synchronisation**, not units. Cross only via `self.dm.localToGlobal(local, global)`, `self.dm.globalToLocal(global, local)`, `vec.getArray()`, `vec.setArray(arr)`. After mutating a `*Local` array view, you MUST `localToGlobal` before the next collective solve, or ranks see stale halos.

## MPI contract
**Collective** (every rank must call, in the same order): `self.dm.localToGlobal`, `self.dm.globalToLocal`, `MPI.COMM_WORLD.Allreduce/Bcast/bcast/Allgatherv/Reduce/Barrier`, `ksp.solve`, `snes.solve`, `ts.solve`, `vec.sum/max/min`, `vec.assemblyBegin/End`, `mat.assemblyBegin/End`, `vec.duplicate/destroy`, `mat.destroy`, `dm.distribute`. **Rank-local**: `vec.getArray()`, `vec.setArray()`, all numpy ops, anything inside `if MPIrank == 0:`.

Two communicator patterns coexist (until unified):
- `MPIcomm = petsc4py.PETSc.COMM_WORLD` — `flow/*`, `sed/*`, `eroder/*` (10 files).
- `MPIcomm = MPI.COMM_WORLD` — `tools/addprocess.py`, `tools/outmesh.py`, `mesher/unstructuredmesh.py`, `mesher/tectonics.py` (4 files).

Rule: use `MPI.COMM_WORLD` for raw collectives (Allreduce/bcast/Allgatherv); use `petsc4py.PETSc.COMM_WORLD` only when creating PETSc objects (`KSP().create(comm=...)`, `Mat().create(comm=...)`). They wrap the same handle but go through different paths inside PETSc.

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

## Eb / EbLocal sign convention (Eb and EbLocal DISAGREE)
**`Eb` (global) is positive for incision.** Set inside `_getEroDepRate` (SPL.py:313 and the equivalent in nlSPL / soilSPL): `E = -tmp/dt` then `self.Eb.setArray(E)`. After SPL completes, nothing else writes `self.Eb`, so it holds the *river* erosion rate, positive for incision, for the rest of the timestep.

**`EbLocal` (local) is positive for DEPOSITION, negative for incision.** Each SPL `erodepSPL*` wrapper *overwrites* `EbLocal` at the end with `add_rate = tmp/dt = -Eb` (SPL.py:367, nlSPL.py:419, soilSPL.py:341). Subsequent kernels `axpy` thickness rates into it (positive = deposit):
- `seaplex.py:486` — marine deposition (positive contributions).
- `hillslope.py:297` — linear hillslope (positive = deposit, negative = erode).
- `soilSPL.py:549` — soil diffusion.

By end-of-step, `EbLocal` is the **net thickness rate**, same sign convention as `cumED`: positive = net deposition, negative = net erosion. The output field `EDrate` (`outmesh.py:272`) writes `EbLocal`, so Paraview shows "positive = deposition".

**The divergence between `Eb` (incision-positive) and `EbLocal` (deposition-positive) is a long-standing inconsistency, not a documented intentional choice.** Treat `Eb` and `EbLocal` as different fields. If you need the river erosion rate in `m/yr` (positive for incision), read `self.Eb`. If you need the net erosion/deposition thickness rate (positive for deposition), read `self.EbLocal`. Do NOT assume they are local/global views of the same field.

Thickness conversion sanity check, in case of future refactors:
- `tmp = -Eb*dt` → positive for deposition, negative for erosion (thickness convention).
- `cumED.axpy(1.0, tmp)` → cumED in thickness convention (positive = net deposition).
- `hGlobal.axpy(1.0, tmp)` → at incision, h drops.

## The `_extra*` methods are mandatory continuations
NOT optional parsers. Each sets attributes required by other modules. Never delete or rename without following the full call chain in `inputparser._readDomain/_readTime/_readHillslope/_readFlex/_readOrography/_readIce`.

- `_extraDomain` (inputparser.py:189) → `seaDepo`, `overlap`, `dataFile`, `nodep`, `strataFile`; calls `_extraDomain2`.
- `_extraDomain2` (:229) → `advscheme`, `radius`, `gravity`.
- `_extraHillslope` (:475) → `nlK`, `clinSlp`, `tsStep`, `Gmar`, `offshore`, `nl_pit_volume/depth/K/inlet_bias`.
- `_extraFlex` (:1430) → `nu`, `flex_res_deg`, `flex_bcN/S/E/W`.
- `_extraOrography` (:1517) → `oro_cw`, `oro_conv_time`, `oro_fall_time`, `oro_precip_*`, `rainfall_frequency`.
- `_extraIce` (:1621) → `iceT`, `elaH`, `iceH` interpolators.

## Forcing DataFrame layout contract
Consumers use `iloc[nb, k]` positional access. **Column order is part of the API. New columns MUST be appended at the END.**

| DataFrame | Built in | Columns (in order) | Positional consumers |
|---|---|---|---|
| `self.tecdata` | `_storeTectonics` (inputparser.py:788) | `start, end, tMap, zMap, hMap` | tectonics.py:73,78,102,112,116 |
| `self.raindata` | `_defineRain` (inputparser.py:1239) | `start, rUni, rzA, rzB, rMap, rKey` | unstructuredmesh.py:677-688 |
| `self.sedfacdata` | `_defineErofactor` (inputparser.py:929) | `start, sUni, sMap, sKey` | unstructuredmesh.py:725-730 |
| `self.tedata` | `_getTe` (inputparser.py:1079) | `start, tUni, tMap, tKey` | addprocess.py:160-174 |

`tectonics.py:78` reads `iloc[nb, -1]` (== `hMap` today). Appending a column to `tecdata` breaks this line silently.

## Magic numbers
These literals MUST NOT be reintroduced inline. New uses should reference this list (and ideally a shared `constants` module if you add one):
- `MISSING_DATA_SENTINEL = -1.0e8` — pre-fill before `Allreduce(MAX)` (addprocess.py:509,588; tectonics.py:430).
- `DISCHARGE_FLOOR = 1.0e-8` — minimum FA/discharge written to output (outmesh.py:288,299,328,357).
- `DEPOSIT_FLOOR = 1.0e-3` (1 mm) — drop sub-mm sediment deposit (seaplex.py:465; sedplex.py:156).
- `BEDROCK_EXPOSED = 1.0e-1` (10 cm) — soil-vs-bedrock threshold (soilSPL.py:166,243).
- `BEDROCK_SENTINEL = 1.0e6` — infinite-bedrock layer-0 thickness (stratplex.py:99,194,225).

`-1.0e6` (seaplex.py:153; iceplex.py:76) is a **boundary marker** for flow-direction Fortran, NOT a missing-data sentinel — keep distinct.

## High-risk modules (do not edit without full regression run)
- **`mesher/unstructuredmesh.py`** — owns `dm`, every shared mesh attribute, forcing dispatch, and `destroy_DMPlex` (lines 738-799) which names every Vec/Mat by hand. Adding a new persistent Vec elsewhere requires adding it to that destroy list or it leaks.
- **`flow/flowplex.py`** — owns `_solve_KSP`, `_solve_KSP2`, `_matrix_build`, `_matrix_build_diag`, `_buildFlowDirection`. Consumed by every downstream module.
- **`tools/inputparser.py`** — owns all parameter parsing; forcing DataFrame column order is API; the `_extra*` chain is mandatory.

## Known bugs (fix before refactoring)
- **`Eb` (global) and `EbLocal` (local) end the timestep in different sign conventions** (see "Eb / EbLocal sign convention" section). `Eb` is positive-for-incision (set by `_getEroDepRate`, never re-synced after the per-step SPL wrapper overwrites `EbLocal` with `-Eb`). `EbLocal` is positive-for-deposition. Any code that calls `dm.localToGlobal(EbLocal, Eb)` or `globalToLocal(Eb, EbLocal)` after `erodepSPL` will produce a sign-flipped field. Test `test_regression.py::test_erosion_sign_conventions` guards the EbLocal convention. No test guards `Eb` directly. Whoever picks this up: decide which convention is canonical, make them consistent, and add an `Eb`-side test.

## Fixed (regression-guarded)
- **`rUni`/`sUni` key/column mismatch** — fixed 2026-06 at `inputparser.py:915` (dict key `'rUni'` → `'sUni'` in the uniform branch of `_defineErofactor`). Before the fix, a YAML with a uniform `sedfactor` event landed `NaN` in `sedfacdata['sUni']`, then crashed `unstructuredmesh._updateEroFactor` on `np.load(None)`. Guarded by `tests/test_regression.py::test_uniform_sedfactor_populates_sUni`.
- **Marine sediment leak at terminal ocean sinks** — fixed 2026-06 at `seaplex._distOcean`. `_matOcean` constructs `dMat1` with a zero column at every terminal sink (cells whose every flow direction is self-referential, typically deep-ocean basin floors); the first `dMat1.mult` therefore annihilated any `sedflux` that landed directly at a sink. `vSed` upstream-integration naturally accumulates large flux at these sinks, so this was the dominant mass-loss path. Diagnostic: 9.748e12 m³ marine sediment supply per step → 4.107e12 m³ landing at sinks (42% of input) → 0 m³ deposited via the routing loop → silently dropped. The fix adds a pre-drain step at the top of `_distOcean` that deposits any sinkVol at terminal sinks directly into `vdep` before the first `dMat1.mult`. The in-loop force-deposit and exit-residual drain (added alongside) handle any sediment that subsequently arrives at sinks during multi-step routing. Guarded by `tests/test_regression.py::test_mass_conservation`. The corresponding `is_sink_local` mask is computed in `_matOcean` once per call.

## Intentional surprises (do NOT "fix")
- **N/S boundary swap in gFlex**. `addprocess.py:221-222` deliberately assigns `simflex.BC_S = self.flex_bcN` and `simflex.BC_N = self.flex_bcS` to compensate for a coordinate-convention mismatch with gFlex. Verified intentional in project memory.

## Fixed (regression-guarded)
- **`rUni`/`sUni` key/column mismatch** — fixed 2026-06 at `inputparser.py:915` (dict key `'rUni'` → `'sUni'` in the uniform branch of `_defineErofactor`). Before the fix, a YAML with a uniform `sedfactor` event landed `NaN` in `sedfacdata['sUni']`, then crashed `unstructuredmesh._updateEroFactor` on `np.load(None)`. Guarded by `tests/test_regression.py::test_uniform_sedfactor_populates_sUni`.

## Checklist before any commit
1. Did you read this file? If invariants here changed, update them.
2. Did you run the regression tests, including the uniform `sedfactor` case (see Known bugs)?
3. If you touched a forcing DataFrame, did you append columns at the END (not in the middle)?
4. If you used a scratch Vec, did you document which ones in the method's docstring?
5. If you changed a method called by `Model.runProcesses` (`model.py:217-286`), did you check every caller AND the mixin init order (`model.py:126-198`)?
