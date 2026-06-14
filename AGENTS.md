# AGENTS.md

Last reviewed 2026-06-13 against `v2026.6.13` (+ dual-lithology and SIA ice-sheet features on `dev`). Read this at the start of every session. Update it when an invariant here changes. See `REFACTOR_AUDIT.md` for the long-form rationale behind each rule.

## What goSPL does
goSPL is a parallel landscape-evolution model that integrates the stream-power law (river incision), linear and non-linear hillslope diffusion, marine sediment transport, glacial accumulation, flexural isostasy, and horizontal/vertical tectonics on an unstructured Voronoi/Delaunay finite-volume mesh. The mesh is either a 2D flat plane (`self.flatModel == True`) or a global sphere; partitioning, halo exchange, and all linear/non-linear solves run on PETSc DMPlex via petsc4py. Time integration is an explicit outer Euler loop in `Model.runProcesses` with implicit KSP/SNES/TS inner solves for diffusion, flow accumulation, and sediment routing.

## The numpy ↔ PETSc boundary
Every state field exists in two parallel representations.

**Numpy land** (raw arrays indexed by local node ID, dimensional, no halo): `self.lcoords` (m), `self.mCoords` (m), `self.larea` (m²), `self.rainVal` (m/yr), `self.upsub` (m/yr), `self.stratH/stratZ/phiS/stratK`, `self.stratHf/phiF` (dual-lithology fine pile — see `## Dual lithology`), `self.fineFrac/depoFineFrac` (dual), `self.pitParams`, `self.pitIDs`, `self.lFill`, `self.localFlex`, plus any `vec.getArray().copy()` view.

**PETSc land** (parallel Vec with halo, mutated via `setArray`/`getArray`/`localToGlobal`): `self.hLocal`/`self.hGlobal`, `self.cumED`/`self.cumEDLocal`, `self.FAL`/`self.FAG`, `self.fillFAL`, `self.Eb`/`self.EbLocal`, `self.bL`/`self.bG`, `self.areaLocal`/`self.areaGlobal`, `self.iceHL`/`self.iceMeltL`/`self.iceUbL`/`self.iceAbrL`/`self.iceFlex` (SIA ice — see `## Ice sheet`), `self.Lsoil`/`self.Gsoil`, `self.lHbed`/`self.gHbed`, `self.vSed`/`self.vSedLocal`, `self.vSedF`/`self.vSedFLocal` (dual-lithology fine flux), `self.fiso`.

Both sides hold physical units; the boundary is about **who owns halo synchronisation**, not units. Cross only via `self.dm.localToGlobal(local, global)`, `self.dm.globalToLocal(global, local)`, `vec.getArray()`, `vec.setArray(arr)`. After mutating a `*Local` array view, you MUST `localToGlobal` before the next collective solve, or ranks see stale halos.

## MPI contract
**Collective** (every rank must call, in the same order): `self.dm.localToGlobal`, `self.dm.globalToLocal`, `MPI.COMM_WORLD.Allreduce/Bcast/bcast/Allgatherv/Reduce/Barrier`, `ksp.solve`, `snes.solve`, `ts.solve`, `vec.sum/max/min`, `vec.assemblyBegin/End`, `mat.assemblyBegin/End`, `vec.duplicate/destroy`, `mat.destroy`, `dm.distribute`. **Rank-local**: `vec.getArray()`, `vec.setArray()`, all numpy ops, anything inside `if MPIrank == 0:`.

**PETSc initialisation happens exactly once**, in `gospl/__init__.py` (`petsc4py.init(sys.argv)`, line 25). Python guarantees the package `__init__` runs before any submodule, so module-level code in submodules (e.g. `MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()` at import time) can rely on PETSc being live. **Do NOT re-introduce `petsc4py.init` in any submodule** — until 2026-06 every submodule called it at import time (15 sites); the call is idempotent so duplicates were harmless but obscured where state was created. Submodules still `import petsc4py` to access `petsc4py.PETSc.X` symbols; that's a separate concern from `init()`.

**`__version__` is set in `gospl/__init__.py`** via `importlib.metadata` reading the installed package metadata:
```python
try:
    from importlib.metadata import version, PackageNotFoundError
    __version__ = version("gospl")
except PackageNotFoundError:
    __version__ = "unknown"
```
The metadata version is driven by `meson.build` line 4 (`version: '2026.6.13'`). `gospl.__version__` derives from it via `importlib.metadata` — never hardcode the version in `__init__.py`. The `PackageNotFoundError` fallback covers the case where the package is cloned but not installed (e.g. bare `git clone` without `pip install -e .`).

**There is one other version literal that MUST be kept in sync** with `meson.build` at every bump: `conda/meta.yaml` line 2 (`{% set version = "..." %}`). conda-build does NOT introspect `meson.build`; it resolves the jinja `version` literally and embeds it in the `.conda` artefact filename. A drift between the two means the conda channel publishes a different version number than the installed `gospl.__version__`, and the workflow's `anaconda upload --skip-existing` silently no-ops if the conda value matches an already-published release. This bit us once between `v2026.6.12` (PyPI+Docker only) and `v2026.6.13`; the inline comment in `conda/meta.yaml` flags it for future contributors. **Bump checklist: change `meson.build:4` AND `conda/meta.yaml:2` together, always.**

**Version spelling convention (adopted 2026-06-12): no leading zeros on month or day** — e.g. `2026.6.13`, not `2026.06.13`. PyPI auto-normalizes per PEP 440 (strips leading zeros for display and in the wheel/sdist filename), so a `2026.06.13` `meson.build` would show up on PyPI as `2026.6.13` while conda artifacts retained the `2026.06.13` spelling — the two channels would visually diverge for the same release. Writing the no-zero form everywhere keeps PyPI display, conda display, `.conda` filename, `.tar.gz` sdist filename, git tag (`v2026.6.13`), and `gospl.__version__` all bitwise-identical. Past tags (`v2026.06.08`, `v2026.06.11`) stay as historical record; do NOT retroactively re-spell them.

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
| 7 | `_IceMesh()` :154 | FAMesh | `iceHL/iceMeltL/iceUbL/iceAbrL/iceFlex` (only if `iceOn`); seeds `iceHL` from `hinit` on fresh start |
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

## Dual lithology (coarse/fine sediment)
Opt-in, active **only** when stratigraphy is on (`stratNb > 0`). Enabled by `strata: dual: True` in YAML; the master flag is `self.stratLith`. **When off (default), every dual path is gated out and behaviour is byte-identical to the single-fraction code** — guarded by `test_dual_all_coarse_matches_single_fraction` (bitwise) and the regression suite. Full design + rationale: `docs/DESIGN_DUAL_LITHOLOGY.md`.

**State (all gated on `stratLith`):**
- `self.stratH` keeps its meaning = **total** layer thickness. `self.stratHf` = fine-fraction bulk thickness per layer (coarse = `stratH − stratHf`); invariant `0 ≤ stratHf ≤ stratH`.
- `self.phiS` = coarse porosity per layer, `self.phiF` = fine porosity per layer.
- `self.fineFrac` (per node) = `vSedF/vSed` of the routed flux (composition arriving at each node). `self.depoFineFrac` = the fraction `deposeStrat` actually uses; seeded to `fineFrac`, then refined inside pits (`_pitFineFraction`) and the marine domain (`_marineFineFraction`).
- `self.vSedF`/`self.vSedFLocal` = fine sediment flux (registered in `destroy_DMPlex`).
- Parsed params (`_extraStrata`): `phi0c/z0c` (coarse curve, defaults to `phi0s/z0s`), `phi0f/z0f` (fine), `fine_k_factor` (fine erodibility multiplier), `fine_diff_factor` (fine diffusivity multiplier), `bedrock_coarse_frac`, `pit_inlet_bias_*`, `fine_efficiency`.

**One shared hook does the lithology coupling:** `_surfaceComposition()` returns the per-node coarse fraction `fc` of the exposed top layer (all-coarse = 1.0 when off). Built on it: `_surfaceLithoK()` (= `fc + ff·fine_k_factor`, folded into `surfK` in all three SPL flavours) and `_surfaceLithoD()` (= `fc + ff·fine_diff_factor`, applied to `Cd` in `_hillSlope`/`_hillSlopeNL`/`soilSPL.diffuseSoil`). Both are 1.0 (no-op) when single-fraction or no contrast.

**Erosion** (`erodeStrat`) splits eroded solid into `thCoarse`/`thFine` by the consumed layers' composition; the bedrock sentinel is inflated *proportionally to layer-0's fine ratio* so it never alters that layer's composition. **Transport**: `_getSedFlux` routes the TOTAL (`thCoarse+thFine`) into `vSed` and `thFine` into `vSedF` (same operator, exact). **Deposition is composition-only — the deposit geometry (`delta`/`dh`) is never modified, so elevations/flow are unchanged.** `_pitFineFraction` (continental) and `_marineFineFraction` (marine) bias the deposit fine fraction toward the depocenter / deep-distal water by bathymetric depth, conserving the deposited fine volume; coarse stays proximal. **Compaction** (`_depthPorosityDual`) compacts each fraction on its own curve, conserving per-fraction solid. **Advection** (`stratalRecord`) advects the fine pile (`stratHf`,`phiF`) with a SECOND `strataonesed` call — NOT `stratathreesed`, whose extra fields are 0–1 fractions and renormalised (`functions.F90:1968`), wrong for the bulk-thickness representation here.

**Fine-enriched overspill (implemented):** `_moveDownstream` threads the fine sub-volume through the cascade with **coarse-settles-first** retention — a filled pit keeps a coarse-enriched deposit (its `_pitRetFine`, used by `_pitFineFraction`) and the excess overspills fine-enriched, routed downstream and on to the marine domain. `vSedFLocal` accumulates the residual fine *mirroring* `vSedLocal` — **it is NOT zeroed in `_distributeSediment`; that was the bug in the first, reverted attempt.** Per-fraction conservation is guarded by **`test_dual_fine_conservation`** (`_fineEroded` vs `_fineDeposited`, the fine analogue of the total-cumED test, which does NOT catch a fine-only leak). **Any future change to the fine routing MUST keep that test green.**

## Ice sheet (SIA)
Opt-in via an `ice` YAML section (`self.iceOn`). It is a true **Shallow-Ice-Approximation** model — the old MFD flow-routing proxy and the explicit reference scheme were removed, so there is **no `flow_model` selector, no `iceFAL`/`iceFAG`/`iceFA`**. Full design: `docs/DESIGN_ICE_SHEET.md`; change log: `docs/ICE_SHEET_SUMMARY.md`. Lives in `flow/iceplex.py` (`IceMesh`).

**Dynamics.** `iceAccumulation` runs the implicit SIA solve `_iceFlowSIA`: a cached Jacobian-free SNES (`ngmres`, CG/HYPRE) on the backward-Euler residual `F(H)=H−Hold−Δt(ṁ−∇·q)` (`_form_residual_ice`), `q` from the `ice_flux` Fortran kernel on the total surface `s=zbed+H`. Unconditionally stable → one solve per goSPL step. `_iceSIAParams` builds the ELA mass-balance ramp `ṁ = rainVal·clip((z−elaH)/(iceH−elaH),·,1)` — **precipitation-scaled** (so ablation magnitude also scales with `rainVal`; it is not a degree-day melt). `_iceSIAFinalize` clamps `H≥0`, removes ice below `max(hterm, sealevel)`, and sets `iceUbL` (basal speed, `ice_velocity` kernel), `iceAbrL` (= `Kg·|ub|^l`), `iceMeltL` (ablation volume re-injected into the river FA by `flowplex.flowAccumulation`).

**Geometry forcing.** `hela`/`hice`/`hterm` are uniform scalars, per-vertex maps, or a `glaciers` time series (`_buildIceSeries` → `_iceTimeSeries`, stepped each forcing update by `mesher._updateIce` into `elaMesh`/`iceMesh`/`termMesh`, indexed `[locIDs]`). Scalar/`evol` path stays byte-identical (keeps the global no-ice short-circuit); `_iceSIAParams` and the `flowplex` precip split are divide-by-zero hardened for per-node `iceH≤elaH`. `hterm` unset = `TERMINUS_UNSET` sentinel → sea level. `hinit` seeds `iceHL` on a fresh start (restart's `readData` overrides). Helper `scripts/ela_from_temperature.py` builds ELA maps from temperature.

**Glacial erosion / till.** With `till.on` off, abrasion `Eg=Kg·|ub|^l` is added to `Eb` as incision (`SPL._glacialAbrasion`, gated off when `till.on`). With `till.on`, `glacialTill` is a conservative transport: it computes a per-cell deposition weight (`Σ owned = 1`) and applies it bulk (`stratNb==0`, net bed change 0) or through the strata (`_glacialTillStrata`: `erodeStrat`→moraine via `deposeStrat`, split coarse/fine under dual lithology, per-fraction solid conserved; **saves/restores `thCoarse`/`thFine`/`depoFineFrac` around the calls** since `sedChange` consumes them later). Weight = melt-weighted by default, or **catchment-routed** (`till.route`) via `_routeTill` — transport-with-loss on the ice surface (`mfdreceivers` steepest descent, melt-out `f=min(1,ablation·dt/H)` forced 1 at the margin so nothing leaks off-ice), normalised to sum 1 exactly. Conservation guarded by `test_ice_glacial_till_conserves` / `_dual_lithology` / `_routing` / `_routing_conserves`. **Any change to till must keep those green.**

## Analysis tools (`gospl.analyse`)
Post-processing subpackage (NOT part of the `Model.runProcesses` pipeline; safe to edit independently). `gospl.analyse.provenance` attributes sediment deposited in sink basins to source-rock classes (per basin + per pixel) with transport distance and a Cu-fertility layer — design + the (shipped) in-model variant in `docs/DESIGN_PROVENANCE.md`, user/tech docs in `docs/tech_guide/provenance.rst`. Standalone (reads goSPL output); does **not** run inside the model. Key facts: routes the per-step net `erodep` down the re-derived flow graph (MFD default), recycling via a per-node deposited "pile"; `GosplOutput` reassembles the **partitioned** HDF5 (one file per rank per step) onto the global input-mesh ordering by KDTree (the same `tree.query(lcoords)` map goSPL uses); the topological sweep is one Numba-`njit`-compatible kernel (`_sweep_impl`) used as both the Python reference and, via `numba`, the fast path (`method='auto'`, identical results). `numba`/`geopandas` are **optional** (`pip install "gospl[analysis]"`); the tool falls back / errors clearly without them. The package is listed in `pyproject.toml` `[tool.setuptools] packages` (it must be, or it won't install). Tests: `tests/test_provenance.py`. A sparse transport-with-loss solve was evaluated and rejected (inexact under sediment under-supply + slower) — do not re-add it.

**In-model provenance tracers (Approach B — functional, phases B0–B4 done).** Opt-in `provenance:` block (`inputparser._extraProvenance`, requires `stratNb>0`): sets `provOn`/`provNb`/`source_class`/`prov_cu_weight`. State (allocated only when `provOn`, like dual lithology): `stratP[node, layer, class]` per-class layer thickness (Σ over classes == `stratH`, seeded to the bedrock `source_class` in `stratplex._initProvenance`); routed sub-flux `vSedP[c]` + `provFrac`/`depoProvFrac` + `_provEroded`/`_provDeposited` in `sedplex`; registered in `destroy_DMPlex`. Provenance is a **passive label** (no K/D/porosity/sorting feedback) — it rides the total-sediment routing. **Phases B0–B3 done** (functional end-to-end): B1 `erodeStrat` splits eroded sediment by class (`provEro`, Σ == total); B2 `_getSedFlux` routes each `vSedP[c]` through `fMati` (through-flux, Σ == `vSed`) → `provFrac`/`depoProvFrac`; B3 `deposeStrat` adds `depo·depoProvFrac` to `stratP` (keeps Σ over classes == `stratH`, ~3e-8) + `_provDeposited`. Guarded by `test_provenance_conservation` (single source ⇒ every layer stays 100% that class; class-0 leak == 0). **B4 done**: `stratalRecord` advects each `stratP[:,:,c]` like `stratHf` (re-normalised Σ==`stratH`); `_outputStrat` writes/`readData` restores `stratP` (lpoints, layers, classes) in the stratal HDF5 (`test_provenance_output_io`). **Multi-source validated** (`test_provenance_multisource`: 2-class map, both tracked, Σ==`stratH` 3e-8, attribution spatially sensible). **B2b is NOT a conservation fix** — `depoProvFrac` always sums to 1 so `stratP` partitions `stratH` exactly for any N; B2b would only refine per-class *spatial attribution* at pit-internal/marine deposits by threading `vSedP[c]` through `_moveDownstream`+`seaChange` (highest-risk path, buggy-revert history) — deferred unless exact pit/marine attribution is needed. Plan: `docs/DESIGN_PROVENANCE.md` §6. (Bare `STRAMesh.__new__` stub tests don't set `provOn`, so `readStratLayers`/`erodeStrat`/`deposeStrat` read it via `getattr(..., False)`.)

## The `_extra*` methods are mandatory continuations
NOT optional parsers. Each sets attributes required by other modules. Never delete or rename without following the full call chain in `inputparser._readDomain/_readTime/_readHillslope/_readCompaction/_readFlex/_readOrography/_readIce`.

- `_extraDomain` (inputparser.py:189) → `seaDepo`, `overlap`, `dataFile`, `nodep`, `strataFile`; calls `_extraDomain2`.
- `_extraDomain2` (:229) → `advscheme`, `radius`, `gravity`.
- `_extraHillslope` (:475) → `nlK`, `clinSlp`, `tsStep`, `Gmar`, `offshore`, `nl_pit_volume/depth/K/inlet_bias`.
- `_extraStrata` (called from `_readCompaction`) → `stratLith` (dual-lithology master opt-in), `phi0c/z0c` (coarse porosity curve, defaults to compaction `phi0s/z0s`), `phi0f/z0f` (fine), `fine_k_factor` (fine erodibility multiplier), `bedrock_coarse_frac`, `fine_efficiency`, `pit_inlet_bias_coarse/fine`, `fine_diff_factor` (fine diffusivity multiplier). All inert while `stratLith` is False; forced False if `stratNb == 0`. See `docs/DESIGN_DUAL_LITHOLOGY.md`.
- `_extraFlex` (:1430) → `nu`, `flex_res_deg`, `flex_bcN/S/E/W`.
- `_extraOrography` (:1517) → `oro_cw`, `oro_conv_time`, `oro_fall_time`, `oro_precip_*`, `rainfall_frequency`.
- `_extraIce` → legacy uniform `iceT`/`elaH`/`iceH` time interpolators (scalar or `evol` CSV). `_readIce` also sets the SIA params (`sia_Aglen/slide/glen`, `ice_Kg`, `ice_abr_l`, `ice_till_on`, `ice_till_route`). An `ice` section now means **SIA** — there is no `flow_model` selector. `hela`/`hice`/`hterm` may be uniform scalars, per-vertex maps, or a `glaciers` time series (`_buildIceSeries` → `_iceTimeSeries`, resolved each step by `mesher._updateIce` into `elaMesh`/`iceMesh`/`termMesh`); `hterm` defaults to the `TERMINUS_UNSET` sentinel → sea level. `hinit` (scalar/map) seeds initial ice. See `## Ice sheet` and `docs/DESIGN_ICE_SHEET.md`.

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
- `1.0e-8` at `sed/stratplex.py:219` (thickness numerical-noise floor).
- `1.0e-3` at `flow/flowplex.py:363` (water-routing convergence), `sed/sedplex.py:139` (sediment-routing convergence), `flow/pitfilling.py:568,622` (minh epsilon nudges), `sed/seaplex.py:507,509` (clinoH 1mm offset).
- `1.0e-2` in `flow/iceplex.py` (`_iceSIAFinalize`, `_routeTill`) is the SIA ice-presence threshold (m): a distinct physical role, not a tolerance.

Each of these is marked with a permanent `# TODO-REFACTOR: value matches X but distinct role; do not replace` comment so future readers know the coincidence is intentional.

## High-risk modules (do not edit without full regression run)
- **`mesher/unstructuredmesh.py`** — owns `dm`, every shared mesh attribute, forcing dispatch (`applyForces` → `_updateRain`/`_updateIce`/…), and `destroy_DMPlex` (~line 849) which names every Vec/Mat by hand. Adding a new persistent Vec elsewhere requires adding it to that destroy list or it leaks (e.g. the dual-lithology `self.vSedF`/`self.vSedFLocal` and the SIA `self.iceMeltL`/`self.iceUbL`/`self.iceAbrL` are registered there).
- **`flow/flowplex.py`** — owns `_solve_KSP`, `_solve_KSP2`, `_matrix_build`, `_matrix_build_diag`, `_buildFlowDirection`. Consumed by every downstream module.
- **`tools/inputparser.py`** — owns all parameter parsing; forcing DataFrame column order is API; the `_extra*` chain is mandatory.

## Known bugs (fix before refactoring)
- _(none currently open)_

## Fixed (regression-guarded)
- **`rUni`/`sUni` key/column mismatch** — fixed 2026-06 at `inputparser.py:915` (dict key `'rUni'` → `'sUni'` in the uniform branch of `_defineErofactor`). Before the fix, a YAML with a uniform `sedfactor` event landed `NaN` in `sedfacdata['sUni']`, then crashed `unstructuredmesh._updateEroFactor` on `np.load(None)`. Guarded by `tests/test_regression.py::test_uniform_sedfactor_populates_sUni`.
- **Marine sediment leak at terminal ocean sinks** — fixed 2026-06 at `seaplex._distOcean`. `_matOcean` constructs `dMat1` with a zero column at every terminal sink (cells whose every flow direction is self-referential, typically deep-ocean basin floors); the first `dMat1.mult` therefore annihilated any `sedflux` that landed directly at a sink. `vSed` upstream-integration naturally accumulates large flux at these sinks, so this was the dominant mass-loss path. Diagnostic: 9.748e12 m³ marine sediment supply per step → 4.107e12 m³ landing at sinks (42% of input) → 0 m³ deposited via the routing loop → silently dropped. The fix adds a pre-drain step at the top of `_distOcean` that deposits any sinkVol at terminal sinks directly into `vdep` before the first `dMat1.mult`. The in-loop force-deposit and exit-residual drain (added alongside) handle any sediment that subsequently arrives at sinks during multi-step routing. Guarded by `tests/test_regression.py::test_mass_conservation`, which runs on the `minimal_model` fixture — a **global-sphere mesh** (`flatModel=False`, `len(idBorders)==0`), so the domain is closed by construction with no boundary outflux possible; the test asserts `|dV_surface|/activity < 1e-4` AND `|dV_cumED|/activity < 1e-4`, i.e. that deposited volume equals eroded volume to within ~6e-5 relative (the floor-effect budget from `DEPOSIT_FLOOR=1e-3` and pit-routing residue, both documented in the test's tolerance rationale). The corresponding `is_sink_local` mask is computed in `_matOcean` once per call.
- **`fine_diff_factor` design/parser drift** — fixed 2026-06 in `_extraStrata`. The dual-lithology diffusivity contrast was changed from absolute `Dc`/`Df` to a `fine_diff_factor` multiplier in the design doc, but the parser still set the obsolete `self.Dc`/`self.Df` and never parsed `fine_diff_factor`. Surfaced when `_surfaceLithoD` (Phase 7) read `self.fine_diff_factor` → `AttributeError` on any real dual run. Now parses `strata.fine_diff_factor` (multiplier on `Cda`/`Cdm`/`nlK`); guarded by `test_dual_lithology_opt_in`.
- **`Eb` vs `EbLocal` sign-convention divergence** — fixed 2026-06 by adopting the thickness-rate convention (positive deposition, negative incision) for both fields. Changes: drop the leading minus in `E = -tmp.getArray()` inside `_getEroDepRate` (SPL.py), `_getEroDepRateNL` (nlSPL.py), `_getEroDepRateSoil` (soilSPL.py); drop the minus in `tmp.setArray(-Eb * dt)` inside the three matching `erodepSPL*` wrappers; negate `Eb` inside `sedplex._getSedFlux` (sedplex.py:63) where the upstream-integration solve still wants an erosion-positive source. Before the fix, `Eb` was incision-positive while `EbLocal` (after the wrapper's `add_rate = tmp/dt` overwrite) was deposition-positive, so the two fields disagreed in sign by end-of-step; the on-disk `EDrate` field was always deposition-positive (because it's written from `EbLocal`), and the restart loader at `outmesh.py:439-440` restored `Eb` in deposition-positive convention from `EDrate` — so fresh-run state and restart-state were inconsistent before this fix. Guarded by `tests/test_regression.py::test_erosion_sign_conventions` (which asserts `EbLocal <= 0` at incising nodes; `Eb` follows the same convention now).

## Intentional surprises (do NOT "fix")
- **N/S boundary swap in gFlex**. `addprocess.py:221-222` deliberately assigns `simflex.BC_S = self.flex_bcN` and `simflex.BC_N = self.flex_bcS` to compensate for a coordinate-convention mismatch with gFlex. Verified intentional in project memory.
- **`gid` argument on `mfdreceivers` / `mfdrcvrs`**. The Fortran kernels take an extra `gid(nb)` argument (per-local-node global vertex ID, from `self.gid` set in `mesher/unstructuredmesh.py` immediately after `getLGMap`). Inside the kernels, the stored slope is perturbed by `val * (1 - 1.0e-15 * gid(n))` before the quicksort tie-break. This makes EXACT slope ties resolve deterministically across MPI decompositions (lower global ID wins). The perturbation is well below KSP-solver precision and does not affect physical results; near-ties driven by KSP floating-point noise are a separate problem and are NOT fixed by this pattern. See the relaxation comment on `test_parallel_correctness` in `tests/test_regression.py` and `fortran/functions.F90:mfdreceivers` for the full rationale. **Do NOT remove the `gid` argument** — without it, `mfdreceivers` ordering depends on local iteration order and the parallel-correctness test diverges by ~10x from current baseline.
- **Platform-dependent KSP-noise floor in `test_parallel_correctness`**. With the `gid` exact-tie-break in place, the remaining sum-FA drift between n=1 and n=2 is dominated by KSP-solver floating-point non-determinism at partition boundaries — and the magnitude differs by platform. Measured on `tests/fixtures/minimal.yml`: macOS-14 (arm64, conda-forge OpenMPI) ≈ **0.3%**, ubuntu-latest (x86_64, conda-forge MPICH) ≈ **1.7%**. Same fixture, same goSPL, same Python — the spread is driven by MPI implementation, BLAS variant, and FMA availability at the platform level. The test tolerance for `rel_sum_fa` is set at **5e-2** (5%) to clear both with headroom; a real routing regression would shift it toward 0.5+, not 0.05. If a future contributor sees these numbers and wants to tighten the tolerance, the correct path is to make KSP-precision halo state bitwise-identical across decompositions (hard PETSc work), NOT to lower the bound and hope. **Note (2026-06):** `environment.yml` now pins the MPI stack to OpenMPI 4.x on all platforms (see the `environment.yml` conventions in `## CI contract`), so ubuntu CI now also runs OpenMPI rather than the MPICH referenced above — the ubuntu figure is historical and the next CI run will measure the new drift; the 5% tolerance still clears it.

## CI contract
- Workflows: `.github/workflows/tests-pr.yml` (fast tier), `tests-slow.yml` (slow regression + analytical benchmarks), `conda-build.yml` (conda package build/publish on tag push), `docker-build.yml` (HPC container build/publish to Docker Hub on tag push — see `## Docker / HPC Container`).
- Matrix in the test + conda workflows: `ubuntu-latest + macos-14 × Python 3.11 + 3.12`. `docker-build.yml` is the exception — a single `ubuntu-latest` / `linux/amd64` job (the HPC target is amd64 only).
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
  - **MPI stack pinned to the OpenMPI 4.x line** (2026-06): `petsc >=3.21,<3.22`, `petsc4py >=3.21,<3.22`, `openmpi >=4.0,<5.0`, `mpi4py >=4.0`, `h5py * mpi_openmpi*`, `hdf5 * mpi_openmpi*`, `vtk-base` (not `vtk`). Same rationale as the conda recipe (see `## Conda Package Validation > Known platform constraints`): OpenMPI 5.x fails at `MPI_Init` on osx-arm64, petsc4py >=3.22 is mpich-only there, and goSPL needs MPI-linked h5py for collective writes. **Side effect:** because env files can't carry per-platform selectors, pinning `h5py * mpi_openmpi*` + `openmpi <5` forces OpenMPI on **all** platforms — including ubuntu CI, which previously resolved MPICH (see the `test_parallel_correctness` platform note). The 5% `rel_sum_fa` tolerance clears OpenMPI on both ubuntu and macOS, so the test still passes; the historical ubuntu-MPICH drift figure predates this pin.
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
- The issue triage automation is designed to be transparent and idempotent: before posting a new automated response, the bot checks for an identical existing triage comment and skips duplicates when a workflow retry would otherwise repeat the same reply.
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
| `v2026.06.08` | 2026-06-08 | `geodels` | `mamba install -c geodels -c conda-forge gospl=2026.06.08` |
| `v2026.06.11` | 2026-06-11 | `geodels` | `mamba install -c geodels -c conda-forge gospl=2026.06.11` |
| `v2026.6.13` | 2026-06-12 | `geodels` | `mamba install -c geodels -c conda-forge gospl` |

`v2026.6.12` is intentionally absent from this table: that tag fired `pypi-publish` and `docker-build` correctly, but `conda/meta.yaml:2` had not yet been bumped from `2026.06.11`, so the conda-build run on the `v2026.6.12` tag produced a `gospl-2026.06.11-*.conda` artefact, which `anaconda upload --skip-existing` silently no-op'd against the already-published 2026-06-11 release. The recovery was `v2026.6.13` with `meson.build` + `conda/meta.yaml` synced. PyPI users on `gospl==2026.6.12` are functionally equivalent to `2026.6.13` (the only diff is the AGENTS.md / docs sweep below); conda users went straight from `2026.06.11` to `2026.6.13`.

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
- **petsc4py ABI mismatch (historical)**: conda-forge previously published
  only a py310 (`np2py310`) osx-arm64 build of `petsc4py`, which worked for
  all computation but segfaulted during MPI finalization when spawned as a
  subprocess via `mpirun`. `test_parallel_correctness` detects this at runtime
  via `_petsc4py_abi_mismatch()` and skips rather than failing. As of the
  2026-06-11 recipe fix, conda-forge **does** publish py311 and py312
  openmpi-linked `petsc4py 3.21.2` builds (`py311h196a43b_0`,
  `py312ha15fc32_0`), which the recipe now pins — so the ABI-mismatch skip no
  longer fires on osx-arm64 and the test actually runs.
- **`test_parallel_correctness` nested-mpirun env leak (FIXED, regression-guarded)**:
  under OpenMPI, `import gospl` → `petsc4py.init()` → `MPI_Init` in the pytest
  parent exports `OMPI_*`/`PMIX_*`/`PRTE_*` env vars; the test's
  `subprocess.run(["mpirun", ...])` inherited them, so OpenMPI thought it was
  already inside an MPI job and silently refused to launch the nested run
  (`rc=1`, empty output). This was masked in CI until 2026-06 (ubuntu used
  MPICH, which is immune; osx-arm64 *skipped* the test on the old py310-only
  petsc4py). It surfaced on all cells once `environment.yml` pinned OpenMPI 4.x
  + petsc4py 3.21.x (real py311/py312 builds → no skip). **Fix:** `run_at_rank`
  now scrubs `OMPI_*`/`PMIX_*`/`PRTE_*`/`OPAL_*` from the child env (preserving
  `OPAL_PREFIX`) before spawning `mpirun`. Harmless under MPICH. The same leak
  affects the published conda package and HPC container (both OpenMPI), so this
  is a real fix, not just a CI patch.
- **Multi-version render**: pass `--variants '{"python": ["3.11"]}'` to
  conda-build to prevent it rendering the recipe for all Python versions in
  the conda-forge global pinnings file (which includes 3.13, incompatible
  with `numpy=1.26`).
- **OpenMPI 4.x pin (osx-arm64)**: the recipe pins `openmpi >=4.0,<5.0` in
  both `host` and `run`. `openmpi` is a pure transitive dependency (pulled in
  via `mpi4py`/`petsc4py`), so without an explicit pin conda-forge resolves to
  `openmpi 5.x`, which fails at `MPI_Init` on macOS with `PML add procs
  failed / Not found (-13)`. The 4.x line initialises cleanly. Because
  `openmpi` is otherwise invisible in the recipe, the pin must be listed
  explicitly — relying on a downstream package to constrain it does not work.
- **petsc4py 3.21.x ceiling (osx-arm64)**: the recipe pins both `petsc
  >=3.21,<3.22` and `petsc4py >=3.21,<3.22` in `host` and `run`. conda-forge
  `petsc4py >=3.22` is **mpich-only** on osx-arm64 — no openmpi-linked variant
  exists — so allowing >=3.22 silently drags in mpich and conflicts with the
  openmpi pin above. `3.21.2` is the last openmpi-linked `petsc4py` build
  published for osx-arm64.
- **h5py must be MPI-linked**: the recipe pins `h5py * mpi_openmpi*` and
  `hdf5 * mpi_openmpi*` in `run`. goSPL performs parallel (collective) HDF5
  writes via h5py; the conda-forge `nompi` variant of h5py is otherwise a
  valid solve and the solver will pick it, but it silently fails on collective
  writes under MPI. Pinning the `mpi_openmpi*` build string forces the
  MPI-linked variant and keeps it consistent with the openmpi 4.x line.

## Docker / HPC Container

This section covers the Singularity/Apptainer container that packages goSPL for
NCI Gadi and Pawsey Setonix. Read it when working on anything under `docker/`.

**Note the MPI-stack difference from the conda package above.** The conda
package (osx-arm64) links **OpenMPI 4.x** and pins `h5py * mpi_openmpi*`. The
HPC container links **MPICH** (hybrid-MPI bind-mount, below) and builds h5py
against that MPICH. The `mpi_openmpi*` build string is OpenMPI-specific and
MUST NOT be carried into the container — OpenMPI is not MPICH-ABI-compatible.

### Target systems
| System | Scheduler | Container runtime | MPI stack | Notes |
|---|---|---|---|---|
| NCI Gadi | PBS | Singularity (module: `singularity`) | Intel MPI (MPICH ABI) | Cannot build images on Gadi — root required |
| Pawsey Setonix | Slurm | Singularity 4.1.0 (module: `singularity/4.1.0-mpi`) | Cray MPI (MPICH ABI) | Ubuntu 24.04 base mandatory (CPE 25.03+) |

### File layout
```
docker/
├── Dockerfile            Multi-stage: MPICH → PETSc → goSPL venv
├── Singularity.def       Native Apptainer definition (alternative build path)
├── build.sh              docker build → .sif conversion; run on local Linux
└── slurm/
    ├── gadi.pbs          PBS job script (NCI Gadi)
    └── setonix.slurm     Slurm job script (Pawsey Setonix)
```
A repo-root `.dockerignore` keeps `.git`/`results/`/caches out of the build
context (the Dockerfile `COPY . /opt/gospl-src`s the whole tree).

### MPI contract inside the container
The container ships MPICH as a build-time placeholder. At runtime, the cluster's
native high-performance MPI (Cray on Setonix, Intel MPI on Gadi) is bind-mounted
over it by the `-mpi` Singularity module flavour. This is the **hybrid MPI**
pattern — do not deviate from it. Consequences:

- `mpi4py`, `petsc4py`, **and parallel HDF5 + h5py** inside the container **must**
  be compiled from source against the container's MPICH (not pre-built
  conda/pip binary wheels). The Dockerfile and Singularity.def both do this
  explicitly (`--no-binary=mpi4py`, `--no-binary=petsc4py`, source HDF5 +
  `HDF5_MPI=ON --no-binary=h5py`, and `pip install --no-deps` for goSPL so the
  binary wheels are never pulled). Do not swap to binary wheels.
- The MPICH version in the Dockerfile (`ARG MPICH_VERSION`, currently `4.2.3`)
  must remain ABI-compatible with the cluster MPI. As of 2026-06, `mpich 4.2.x`
  is compatible with both Cray MPI and Intel MPI. If you bump it, verify ABI
  compatibility with Setonix's Cray MPI via Pawsey docs before merging.
- On Setonix, load `singularity/4.1.0-mpi`; this module sets `SINGULARITY_BINDPATH`
  and `SINGULARITYENV_LD_LIBRARY_PATH` automatically — no manual `--bind` for MPI
  paths is needed.

### goSPL-specific container invariants
These are required for correctness inside the container and must not be changed:

1. **`petsc4py.init` is called exactly once**, in `gospl/__init__.py`. The container
   does not alter this. Do not add `petsc4py.init()` calls to entrypoint scripts,
   Slurm wrappers, or any other layer around goSPL.

2. **Threading must be disabled inside each MPI rank.** Both the Dockerfile/def
   (`ENV`) and the job scripts export:
   ```bash
   export OMP_NUM_THREADS=1
   export OPENBLAS_NUM_THREADS=1
   export MKL_NUM_THREADS=1
   ```
   Do not remove these.

3. **Collective calls must span all ranks.** Nothing in the job scripts or
   entrypoint may cause one rank to diverge from the MPI contract described in
   `## MPI contract` above. In particular, do not redirect `stdout` of a subset
   of ranks before a `Barrier` or `Allreduce`.

4. **`destroy_DMPlex` must run on clean exit.** goSPL model runs should always
   call `model.destroy()` (the job scripts' `-c` driver does:
   `Model(...).runProcesses(); m.destroy()`). Do not `SIGKILL` the job before
   PETSc finalization unless debugging — leaked Vecs cause HPC sysadmin
   complaints and false leak reports.

5. **The scratch Vecs (`self.tmp`, `self.tmpL`, `self.tmp1`, `self.h`, `self.hl`,
   `self.dh`, etc.) are not safe to read after the step that wrote them.** If a
   post-processing or checkpointing script runs inside the container alongside
   goSPL, do not read these Vecs for elevation or cumulative ED — use
   `self.hLocal/hGlobal` and `self.cumED/cumEDLocal` respectively (see
   `## Scratch vector contract`).

### Known Setonix issue: parallel I/O inside Singularity
Pawsey has an active known issue: MPI-parallel HDF5 collective I/O inside
Singularity on Setonix can fail. goSPL uses collective HDF5 writes via the
MPI-linked h5py built into the container (against the container's MPICH — see
the MPI contract above). If you hit this on Setonix:
- First check Pawsey's Known Issues page for resolution status.
- Workaround: switch goSPL output to independent I/O mode if supported, or
  post-process outside the container after the simulation.
- Do NOT work around this by switching to the `nompi` h5py variant — goSPL
  requires MPI-linked h5py for collective writes. (Inside the container that
  means the MPICH-linked build, NOT the conda `mpi_openmpi*` build — see the
  note at the top of this section.)

### Build workflow (local Linux machine)
```bash
cd docker/
./build.sh v2026.6.13           # build + smoke-test + convert to .sif (no push)
./build.sh v2026.6.13 --push    # build + push Docker Hub + convert to .sif
```
The `build.sh` script:
1. Runs `docker build --platform linux/amd64` from the repo root with
   `--build-arg GOSPL_VERSION=<tag without leading v>`.
2. Smoke-tests the image (`import gospl; import petsc4py; import mpi4py` →
   `goSPL container build OK`) and asserts `gospl.__version__` matches the tag.
3. Optionally pushes to `docker.io/geodels/gospl-hpc` (override the repo with
   `REGISTRY=... ./build.sh ...`).
4. Converts to `.sif` via `apptainer build` / `singularity build` from the
   local Docker daemon (`docker-daemon://`), or from the registry (`docker://`)
   when `--push` was given.

**CI**: `.github/workflows/docker-build.yml` builds the same image on every `v*`
tag push and pushes `docker.io/geodels/gospl-hpc:<tag>` + `:latest` to Docker
Hub (needs `DOCKERHUB_USERNAME` / `DOCKERHUB_TOKEN` secrets). `workflow_dispatch`
builds + smoke-tests only (no push). The from-source MPICH/PETSc/HDF5 build is
slow (~30-60 min cold); the workflow uses the gha build cache.

**You cannot build on Gadi or Setonix.** Build the `.sif` locally and `scp` it,
or — once the image is on Docker Hub — `singularity pull docker://geodels/gospl-hpc:v2026.6.13`
on a login node (a conversion, not a root build).

numpy (`1.26.4`) and cython (`<3.1`) are pinned for the whole venv via a
constraints file wired to both `PIP_CONSTRAINT` and `PIP_BUILD_CONSTRAINT`
(both env vars because pip ≥26.2 stops honouring `PIP_CONSTRAINT` for build
deps): numpy stays on the conda/CI baseline (no transitive drift to 2.x), and
petsc4py 3.21.x cannot be cythonized by Cython ≥3.1 (`cyautodoc`
"ExpressionWriter" crash; h5py is fine with <3.1 too).

**setuptools is deliberately NOT pinned globally** — petsc4py 3.21.x and h5py
have *opposite* requirements: petsc4py's `confpetsc.py` needs the classic
`distutils.util.execute(dry_run=...)` that setuptools ≥74 dropped (py3.12 has no
stdlib distutils), while h5py's build requires setuptools ≥77.0.1. They're
resolved per-package: **petsc4py gets a dedicated build-constraint file** pinning
`setuptools<74` (passed as `PIP_BUILD_CONSTRAINT` for just that `pip install`),
while **h5py keeps build isolation** under the global (setuptools-free) constraint
and pulls its own ≥77. (Note: `PIP_BUILD_CONSTRAINT`/`--build-constraint` cannot
be combined with `--no-build-isolation` — pip rejects it — so petsc4py stays in
build isolation.) These build-toolchain contortions are a consequence of building
the older petsc4py 3.21.x from source; a future PETSc bump may let them be relaxed.

The base layers are split into cached stages — `mpich-build` → `petsc-build` →
`hdf5-build` (HDF5 chained `FROM petsc-build`). `docker-build.yml` warms them
with `--target hdf5-build` into a dedicated gha cache scope so a failure in the
python layers never triggers a cold ~45-min base rebuild.

### Version bumping
When bumping the goSPL version in the container, change **only** the
`GOSPL_VERSION` build arg (`build.sh` derives it from the `<tag>` you pass). Do
not hardcode version strings in the Dockerfile/def — goSPL's version derives
from `gospl.__version__`, set via `importlib.metadata` reading `meson.build`
line 4 (see `## MPI contract > __version__`). The build asserts that the
installed `gospl.__version__` equals `GOSPL_VERSION`, so the build arg and
`meson.build` must match on any published `.sif`.

### Checklist before publishing a new .sif
1. `./build.sh <tag>` completes without errors and the smoke test prints
   `goSPL container build OK`.
2. The MPICH version inside the container is ABI-compatible with the cluster MPI
   on both Gadi and Setonix.
3. `OMP_NUM_THREADS=1` / `OPENBLAS_NUM_THREADS=1` are exported in both job
   scripts.
4. The `.sif` filename includes the version tag (e.g. `gospl-hpc-v2026.6.13.sif`).
5. Update the `## Milestones` table in this file with the new tag and `.sif`
   publication date.
6. The `docker/slurm/gadi.pbs` and `docker/slurm/setonix.slurm` `CONTAINER=`
   paths reference the new `.sif` filename.

## Milestones

| Date | Tag | Description |
|---|---|---|
| 2026-06-08 | `refactor-baseline-2026-06` | Tier 2 AI-readability refactor complete. AGENTS.md written, 6 regression tests passing, 3 scientific bugs fixed (rUni/sUni, marine sediment leak, Eb sign convention), constants.py, _get_param, named DataFrame access, KSP lifecycle documented, HOW_TO_ADD_FORCING.md and HOW_TO_ADD_OUTPUT.md written. |
| 2026-06-08 | `v2026.06.08` | First release from refactored codebase. Analytical benchmark suite integrated and green on all CI cells (ubuntu-latest + macos-14 × Python 3.11 + 3.12). Published to `geodels` conda channel. |
| 2026-06-11 | — | `gospl.__version__` added via `importlib.metadata`; single source of truth is `meson.build`. |
| 2026-06-11 | `v2026.06.11` | Conda recipe fixes osx-arm64 OpenMPI 5.x `MPI_Init` failure: explicit `openmpi >=4.0,<5.0`, `petsc`/`petsc4py >=3.21,<3.22` (last openmpi-linked builds), `h5py`/`hdf5 * mpi_openmpi*`, relaxed `mpi4py >=4.0`. Version bumped 2026.06.08 → 2026.06.11 for the re-publish. |
| 2026-06-11 | — | `docker/` HPC container added (Dockerfile, Singularity.def, build.sh, slurm/gadi.pbs, slurm/setonix.slurm) for NCI Gadi + Pawsey Setonix: hybrid-MPI MPICH-ABI build with mpi4py/petsc4py/parallel-HDF5+h5py compiled from source. Not yet built/published (`.sif` TBD — build is local-Linux only). |
| 2026-06-12 | `v2026.6.12` | First PyPI sdist release via `.github/workflows/pypi-publish.yml` (Trusted Publishing / OIDC, sdist-only, no wheels). No-leading-zero version-spelling convention adopted (`2026.6.12` not `2026.06.12`) for PyPI ↔ conda artefact-name parity. **conda channel not bumped on this tag** — `conda/meta.yaml:2` still referenced `2026.06.11`, so the conda-build job silently skipped the upload (`anaconda upload --skip-existing` matched the existing 2026-06-11 release). PyPI + Docker Hub `:v2026.6.12` published correctly. |
| 2026-06-12 | `v2026.6.13` | Sync release: bumped `conda/meta.yaml:2` to track `meson.build:4`; added inline cross-reference comment; AGENTS.md > `__version__` corrected to acknowledge `conda/meta.yaml` as a second version literal that must be bumped in lockstep. Docs sweep: every `v2026.06.11` example in `docs/getting_started/install{HPC,Conda,Docker}.rst`, `docker/build.sh`, `docker/Singularity.def`, and AGENTS.md (example commands + `.sif` filename + Released-packages table) bumped to `v2026.6.13`. New `docs/getting_started/installPip.rst` documents the PyPI-sdist install path; wired into the `getting_started/index.rst` grid + toctree. |
| 2026-06-12 | — (on `dev`) | **Dual-lithology (coarse/fine) sediment** — core merged into `dev` (PRs #415 co-transport, #416 separate transport/diffusivity, #417 docs). Opt-in via `strata: dual`; single-fraction byte-identical. Phases 0–7: `_extraStrata` parser, `stratHf`/`phiF` state + I/O/restart, `_surfaceComposition`/`_surfaceLithoK`/`_surfaceLithoD` hooks, `erodeStrat` split, separate fine flux (`vSedF`), composition-only depth-biased lake (`_pitFineFraction`) + marine (`_marineFineFraction`) deposition, per-fraction compaction (`_depthPorosityDual`), fine-pile advection, composition-weighted hillslope + soil diffusivity. Bitwise parity + closed-sphere mass-conservation tests. Full plan: `docs/DESIGN_DUAL_LITHOLOGY.md`. |
| 2026-06-13 | — (on `dev`) | **Dual-lithology completion** — both retrospective limitations resolved + I/O round-trip (PRs #418–#420). #418: fine-specific conservation guard (`test_dual_fine_conservation`, `_fineEroded`/`_fineDeposited`) THEN fine-enriched overspill — `_moveDownstream` threads the fine sub-volume with coarse-settles-first retention (`_pitRetFine`), excess overspills fine-enriched to the distal marine basin (~98% of deposited fine reaches the sea; conserved to 3e-5). The prior reverted attempt's bug: zeroing `vSedFLocal` instead of mirroring `vSedLocal`. #419: write the fine sediment load to output (`sedLoadF`, XMF `SLf`; coarse = `sedLoad − sedLoadF`). #420: robust + documented per-layer initial composition via npstrata `strataHf`/`phiF` (clamp `0 ≤ strataHf ≤ strataH`, `phiF` defaults to `phi0f`, warn if dual off). Suite: 25 passed, 1 skipped. No correctness gaps remain — only enhancements (true two-field differential-transport solve; quantitative basin validation). |
| 2026-06-13 | — (on `dev`) | **SIA ice-sheet model** — true Shallow-Ice-Approximation replacing the MFD proxy (PRs #423–#425 dynamics+velocity / abrasion+till / loading+output, then #426–#429). Implicit non-linear-diffusion thickness solve (`ice_flux` kernel, cached SNES), velocity-based abrasion, glacial till→moraine (bulk + stratigraphic/dual-lithology, conservative), meltwater re-injected into the river FA, ice loading via existing flexure. #426: MFD + explicit schemes removed → SIA is the sole model; till coupled to dual lithology; `iceMelt`/`iceAbr` outputs; spatial + time-series ELA maps (`glaciers`); `hinit` seed-and-evolve. #427: terminus floor `max(hterm, sea level)`, default sea level. #428: opt-in catchment-aware till routing (`till.route`, `_routeTill`). #429: `scripts/ela_from_temperature.py` (ELA from paleo temperature). Opt-in via `ice:`; ice-off byte-identical. Full design: `docs/DESIGN_ICE_SHEET.md`; change log: `docs/ICE_SHEET_SUMMARY.md`. Suite: 43 passed, 1 skipped. |
| 2026-06-14 | — (on `dev`) | **Sediment-provenance tool** (`gospl.analyse.provenance`, PRs #431/#432 + follow-up) — standalone source-to-sink attribution: deposited sediment → source-rock classes per basin/pixel, transport distance, Cu-fertility layer. MFD fraction routing with recycling; partition-aware HDF5 reassembly (`GosplOutput`, KDTree to global mesh); per-step HDF5 + XDMF (ParaView) + per-basin CSV; Numba-JIT sweep fast path (`method='auto'`); `classes_from_shapefile` (geopandas). `numba`/`geopandas` optional (`gospl[analysis]`); `gospl.analyse` added to `pyproject.toml` packages. Sparse transport-with-loss solve evaluated and rejected (inexact under under-supply + slower). Not in the model pipeline. Design: `docs/DESIGN_PROVENANCE.md`; tech doc: `docs/tech_guide/provenance.rst`. Suite: 51 passed, 1 skipped. |
| 2026-06-14 | — (on `dev`) | **In-model provenance tracers** (Approach B, PRs #434–#438) — opt-in `provenance:` carries N source-rock classes through the model's own erosion (`erodeStrat` split) → transport (`_getSedFlux` routes each `vSedP[c]` through `fMati`) → deposition (`deposeStrat` writes `stratP`) → stratigraphy (advected like `stratHf`, written/restored in the stratal HDF5 as `stratP`). Passive label (no K/D/sorting feedback); requires stratigraphy; ice-off byte-identical. **Conservation exact for any N sources** (`stratP` partitions `stratH` ~3e-8 since `depoProvFrac` sums to 1); single-source + 2-class-map tests. Optional B2b (per-class pit/marine *attribution* refinement — not conservation) deferred. Design/tech/user docs updated. Suite: 57 passed, 1 skipped. |

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
11. If you bumped the version, did you change **both** `meson.build:4` AND `conda/meta.yaml:2` in the same commit, with the same no-leading-zero spelling (e.g. `2026.6.13`, not `2026.06.13`)? (`gospl/__init__.py`, `pyproject.toml`, and `gospl.__version__` derive from `meson.build` via `importlib.metadata`; conda-build does NOT — it reads the literal in `meta.yaml`. A drift between the two causes a silent no-op on `anaconda upload --skip-existing`. See AGENTS.md > `__version__`.)