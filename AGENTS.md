# AGENTS.md

Rules for anyone (human or AI agent) changing goSPL. Read this file at the start of
every session, and the subsystem file in the **Map** below before you edit that part
of the tree. Update the relevant file when an invariant changes. Long-form rationale:
`REFACTOR_AUDIT.md`; past bugs and their guards: `docs/dev/FIXED_BUGS.md`; per-change
history: `docs/dev/CHANGELOG_DEV.md`.

Last reviewed 2026-10-03 against `dev` (after `v2026.7.14`). In October 2026 this file
was split: it keeps only the rules that apply everywhere, and the subsystem detail
moved verbatim into the files below.

## Map

| Before you edit… | Read |
|---|---|
| `gospl/flow/` (flow KSP, pits, cascade, ice, groundwater) | `gospl/flow/AGENTS.md` |
| `gospl/sed/` (sediment, marine, hillslope, strata, dual lithology, provenance tracers) | `gospl/sed/AGENTS.md` |
| `gospl/eroder/` (SPL, nlSPL, soilSPL) | `gospl/eroder/AGENTS.md` |
| `gospl/tools/` (`inputparser`, forcing DataFrames, `zprobe`, outputs) | `gospl/tools/AGENTS.md` |
| `gospl/analyse/` (post-processing tools) | `gospl/analyse/AGENTS.md` |
| `fortran/` (kernels, neighbour tables, memory safety) | `fortran/AGENTS.md` |
| `tests/` | `tests/README.md` |
| `benchmarks/` | `benchmarks/AGENTS.md` |
| `.github/`, `environment.yml` | `.github/AGENTS.md` |
| `conda/` | `conda/AGENTS.md` |
| `docker/` | `docker/AGENTS.md` |
| a version bump / release | `docs/dev/RELEASE.md` and the `release` skill |
| a new process, output or forcing | `docs/HOW_TO_ADD_FORCING.md`, `docs/HOW_TO_ADD_OUTPUT.md`, the `add-feature` skill |

Recurring procedures are written up as skills in `.claude/skills/` (`release`,
`debug-spike`, `parallel-check`, `add-feature`); they are plain Markdown, usable by any
agent.

## Working loop (commands)

| Task | Command |
|---|---|
| Quick check (166 tests, ~25 s) | `pytest tests/ -n 4 -m "not mpi"` |
| One subsystem | `pytest tests/ -m flow` (markers: `tests/README.md`) |
| Multi-rank guards (~70 s) | `pytest tests/ -m mpi` |
| Full suite before a commit | `pytest tests/ -n 4` (CI runs it serially) |
| Analytical benchmarks (minutes) | `pytest benchmarks/ -m benchmark` |
| Collective-under-rank-guard lint | `python scripts/lint_mpi_collectives.py` (CI runs it) |
| np=1 vs np=N field comparison | `python scripts/ab_partition.py <input.yml> -n 4` |
| Mass budget of a run | `python scripts/budget.py <output-dir>` |
| Machine-readable per-step summary | `gospl -i input.yml --summary run.jsonl` |
| Inspect / compare output | `gospl-inspect <output-dir> [--step N] [--compare <other-dir>] [--json]` (or `python -m gospl.analyse.runinspect` before a `pip install -e .` registers the command) |
| Fortran bounds-checked suite | `scripts/fcheck.sh` |
| goSPL-examples smoke run | `python scripts/examples_smoke.py <examples-dir> --only Local-examples` |
| CI runs / logs | `gh run list`, `gh run view <id> --log-failed` (`gh` is authenticated) |

## What goSPL does
goSPL is a parallel landscape-evolution model that integrates the stream-power law (river incision), linear and non-linear hillslope diffusion, marine sediment transport, glacial accumulation, flexural isostasy, and horizontal/vertical tectonics on an unstructured Voronoi/Delaunay finite-volume mesh. The mesh is either a 2D flat plane (`self.flatModel == True`) or a global sphere; partitioning, halo exchange, and all linear/non-linear solves run on PETSc DMPlex via petsc4py. Time integration is an explicit outer Euler loop in `Model.runProcesses` with implicit KSP/SNES/TS inner solves for diffusion, flow accumulation, and sediment routing.

## The numpy ↔ PETSc boundary
Every state field exists in two parallel representations.

**Numpy land** (raw arrays indexed by local node ID, dimensional, no halo): `self.lcoords` (m), `self.mCoords` (m), `self.larea` (m²), `self.rainVal` (m/yr), `self.upsub` (m/yr), `self.stratH/stratZ/phiS/stratK`, `self.stratHf/phiF` (dual-lithology fine pile — see `## Dual lithology`), `self.fineFrac/depoFineFrac` (dual), `self.pitParams`, `self.pitIDs`, `self.lFill`, `self.localFlex`, plus any `vec.getArray().copy()` view.

**PETSc land** (parallel Vec with halo, mutated via `setArray`/`getArray`/`localToGlobal`): `self.hLocal`/`self.hGlobal`, `self.cumED`/`self.cumEDLocal`, `self.FAL`/`self.FAG`, `self.fillFAL`, `self.Eb`/`self.EbLocal`, `self.bL`/`self.bG`, `self.areaLocal`/`self.areaGlobal`, `self.iceHL`/`self.iceMeltL`/`self.iceUbL`/`self.iceAbrL`/`self.iceFlex` (diagnostic glacial model — see `## Ice sheet`), `self.Lsoil`/`self.Gsoil`, `self.lHbed`/`self.gHbed`, `self.vSed`/`self.vSedLocal`, `self.vSedF`/`self.vSedFLocal` (dual-lithology fine flux), `self.fiso`.

Both sides hold physical units; the boundary is about **who owns halo synchronisation**, not units. Cross only via `self.dm.localToGlobal(local, global)`, `self.dm.globalToLocal(global, local)`, `vec.getArray()`, `vec.setArray(arr)`. After mutating a `*Local` array view, you MUST `localToGlobal` before the next collective solve, or ranks see stale halos.

## MPI contract
**Collective** (every rank must call, in the same order): `self.dm.localToGlobal`, `self.dm.globalToLocal`, `MPI.COMM_WORLD.Allreduce/Bcast/bcast/Allgatherv/Reduce/Barrier`, `ksp.solve`, `snes.solve`, `ts.solve`, `vec.sum/max/min/norm/dot`, `mat.norm`, `vec.assemblyBegin/End`, `mat.assemblyBegin/End`, `vec.duplicate/destroy`, `mat.destroy`, `dm.distribute`, `petsc4py.PETSc.garbage_cleanup`. **Rank-local**: `vec.getArray()`, `vec.setArray()`, all numpy ops, anything inside `if MPIrank == 0:`.

**THE #1 parallel deadlock class — a collective MUST NOT be gated on a rank-local condition.** If a collective (anything in the list above) is reached on some ranks but not others, the ranks that called it block forever in the underlying `MPI_Allreduce`/scatter waiting for the ranks that skipped — the run **hangs at np>1** (and serial ALWAYS passes, because there is no partner to wait for: the classic serial-passes/parallel-hangs signature). The trap is subtle because the gating condition often *looks* fine:
- **Rank-local guards that desync**: `if MPIrank == 0:` / `!= 0`, `if local_mask.any()`, `if len(self.xPts) > 0`, `if seaID...`, an early `return`/`continue`/`break` on a per-rank value, or a `for`/`while` whose trip-count differs per rank. A `vec.sum()`/`vec.norm()`/`mat.norm()` looks innocent but is collective — putting one inside `if MPIrank == 0` deadlocks (this was the 2026-06-24 `flowplex._solve_KSP2` hang).
- **Global guards that are SAFE**: a YAML/config scalar (`stratNb`, `fDepa`, `flexOn`, …), a global per-pit array (`pitParams`-derived `pitVol`/`eV`; note `eV=inV−pitVol` is safe only because `inV` is `Allreduce`d first), an already-`Allreduce`d value, or a quantity that is *itself* collective (`vec.sum()` returns the same global number on every rank, so `if vec.sum() > x:` is consistent). The same `.any()`/`.size`/`for` is safe iff its operand is one of these.

**The fix recipe**: compute the collective on **every** rank, then gate only the rank-local part (the `print`/IO). For a genuine rank-local decision that must gate a collective, **reduce it first**: `do_it = MPI.COMM_WORLD.allreduce(local_flag, op=MPI.LOR)` (or `LAND`), then branch on `do_it`. The `bcast(..., root=0)` pattern (rank 0 sends, the `else` branch receives) is fine — both sides call `bcast`.

**THE #2 partition-dependence class — matrix/vector assembly MUST set OWNED rows only, never ghost rows.** `Mat.setValuesLocal(row, …)` maps a *local* row through the lgmap; a GHOST local row maps to a global row owned by **another** rank, and its stencil computed on the ghosting rank uses an **incomplete** neighbour set (the ghost's full neighbourhood is not all present locally). Setting it ships that wrong value off-process, where under `INSERT_VALUES` it **collides** with the owner's correct value — the result is **undefined / partition-dependent** (owner vs ghoster race). So loop only over owned rows (`self.glIDs`), exactly as `_evalJacobianMardDiff` and every residual (`F.setArray(res[self.glIDs])`) do. **This does NOT hang** — it silently corrupts the operator's boundary rows, so it passes serial (one partition → no ghost rows) and produces a wrong-but-finite result at np>1 (the 2026-06-24 soil-diffusion `_evalJacobianSoil` spike: it looped `range(self.lpoints)`, incl. ghost rows → partition-dependent Jacobian → an isolated elevation spike at sub-domain boundaries by the 2nd output; np=1 clean). Diagnostic tell: an **exact** solve (swap the PC for `lu`/MUMPS) still spikes ⇒ the *matrix*, not the solver, is partition-dependent. Don't chase the linear solver (KSP/PC/tolerance) until you've confirmed assembly is owned-rows-only.

**Both classes were audited clean on 2026-06-24** (every `setValues*` site, and every collective reached under a condition or inside a loop; full records in `docs/dev/FIXED_BUGS.md` > MPI audits). Safe assembly patterns: (a) **diagonal-only** matrices; (b) **receiver-based** stencils built from the partition-invariant drainage arrays (`rcvIDi`/`wghtVali`); (c) **additive FV-Laplacian** assembly with `ADD_VALUES`. Never `INSERT` a neighbour-summing row over all `lpoints`. The live reductions that make gated collectives safe: the coastline rebuild `allreduce(rebuild, MPI.LOR)` (`hillslope.py`), the closed-sink deposit `allreduce(has_closed, MPI.LOR)` (`sedplex.py`), the per-pit `inV`/`inlet_count` `Allreduce`s. **Run `python scripts/lint_mpi_collectives.py` whenever you add or move a collective**; end-to-end guards: `pytest -m mpi`.

**PETSc initialisation happens exactly once**, in `gospl/__init__.py` (`petsc4py.init(sys.argv)`, line 25). Python guarantees the package `__init__` runs before any submodule, so module-level code in submodules (e.g. `MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()` at import time) can rely on PETSc being live. **Do NOT re-introduce `petsc4py.init` in any submodule** — until 2026-06 every submodule called it at import time (15 sites); the call is idempotent so duplicates were harmless but obscured where state was created. Submodules still `import petsc4py` to access `petsc4py.PETSc.X` symbols; that's a separate concern from `init()`.

**Versions** live in four places that move together: `meson.build:4`, `conda/meta.yaml:2`, `docs/conf.py` and `docs/_static/version_switch.json`. CalVer with no leading zeros (`2026.6.13`, not `2026.06.13`). `gospl.__version__` derives from `meson.build` via `importlib.metadata`; never hardcode it. Background: `docs/dev/RELEASE.md`.

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
PETSc solvers follow two intentional patterns. **Use the right one for new code.** Full per-solver notes (prefixes, why each PC, failure modes) are in the subsystem files.

### CACHED — hot-path solvers (14 sites)
Created lazily on first use, stored as `self._X`, reused for the whole run (avoids ~5–10 ms create/destroy per call).

| File | Method | Cached attribute | Solver |
|---|---|---|---|
| `flow/flowplex.py` | `_solve_KSP` | `_ksp_main` | KSP `fgmres`+`bjacobi` (`flowacc_`); fatal solve full budget, cascade capped |
| `flow/flowplex.py` | `_solveIDAExact` | `_ksp_exact` | KSP `fgmres`+`bjacobi`/exact `lu` (`flowaccx_`); routing-solve rescue, sticky |
| `flow/flowplex.py` | `_solve_KSP2` | `_ksp_fallback` | KSP `richardson`+`none` bounded fallback (`flowaccfb_`) |
| `flow/gwplex.py` | `_solveHead` | `_ksp_gw` | KSP `fgmres`+hypre BoomerAMG (`gw_`); AMG is required |
| `eroder/nlSPL.py` | `_solveNL_ed` | `_snes_ed` (+`_fb`) | SNES `qn`, `ngmres` fallback |
| `eroder/nlSPL.py` | `_solveNL` | `_snes_nl` | SNES nrichardson + analytic Jacobian |
| `eroder/soilSPL.py` | `_solveSoil` | `_snes_soil` (+`_fb`) | SNES `qn`, `ngmres` fallback, adaptive sub-stepping |
| `eroder/soilSPL.py` | `diffuseSoil` | `_ts_soil` | TS rosw |
| `sed/hillslope.py` | `_hillSlopeNL` | `_snes_hill` | SNES |
| `sed/hillslope.py` | `_diffuseImplicit` | `_ts_marine` | TS rosw (marine + lake); stage solve gmres+gasm (NOT preonly: rosw has no Newton correction) |
| `sed/hillslope.py` | `_solveSmooth` | `_ksp_smooth` (+`_smoothMat`) | KSP, operator rebuilt on coastline move |
| `sed/hillslope.py` | `_hillSlope(smooth=0)` | `_ksp_hill_lin` (+`_hillMat`) | KSP, operator rebuilt on coastline move |
| `sed/hillslope.py` | `_diffuseImplicitPicard` | `_ksp_picard` | KSP (opt-in `marineSolver: picard`) |
| `sed/hillslope.py` | `_diffuseProvTracers` | `_ksp_provdiff` | KSP (only when `provOn`) |

**CRITICAL — collective rebuild decision.** Both coastline-gated caches above (`_smoothMat`/`_hillMat`) rebuild the operator (and redo `PCSetUp`) only when `seaID` moves. `seaID` is **rank-local**, but `_buildDiffMat` assembly and `PCSetUp` are **collective** — so the "did the coastline move?" test MUST be reduced across ranks (`rebuild = MPI.COMM_WORLD.allreduce(local_changed, op=MPI.LOR)`) before it gates those calls. Without the reduce, one rank rebuilds while another reuses, the two take different collective paths, and the run **deadlocks at np>1** (it ran fine until the per-partition coastlines drifted apart, then hung). A rank forced to rebuild with an unchanged mask reproduces its own cached matrix, so the reduce stays bit-faithful. Any future cached operator gated on a rank-local condition needs the same reduce.

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

**CRITICAL**: any new cached solver MUST be added to the `destroy_DMPlex` loop in `mesher/unstructuredmesh.py` (`def destroy_DMPlex`, runs to EOF). The loop iterates over a hardcoded list of attribute names; forgetting to add yours leaks the PETSc object at simulation end. Same applies to cached helper Vecs (`self._snes_X_f`, `self._snes_X_x`, `self._snes_X_J`, etc.).

### AD-HOC — nested-matrix fieldsplit solves (2 sites)
Solvers that build a `Mat().createNest(...)` whose sub-matrices change every call AND configure `pc.setType("fieldsplit")` with IS sets derived from the nested-mat structure. The fieldsplit PC's IS configuration is tied to the specific sysMat instance, so re-using a cached KSP via `setOperators(new_sysMat)` does NOT automatically re-derive the splits. Caching is theoretically possible but requires careful experimentation with `pc.reset()` and the nested-mat lifecycle.

| File | Method | Condition |
|---|---|---|
| `eroder/SPL.py` | `_coupledEDSystem` | `self.fDepa != 0` (transport-limited branch with non-zero `G`) |
| `sed/seaplex.py` | `_depMarineSystem` | `not flatModel AND self.Gmar > 0`, AND only inside the second `_distOcean` pass |

**Lifecycle**: create at the top of the method, configure, solve, then explicitly destroy everything (KSP, sub-KSPs, sub-ISes, PC, sysMat, RHS/solution vectors) at the end. See `SPL.py:116-259` (`_coupledEDSystem`) for the canonical pattern. Both sites are COLD-path (conditional, not every step), so the cumulative create+destroy overhead is small.

If a future contributor wants to convert one of these to CACHED, the obstacle is the fieldsplit-PC + nested-mat IS lifecycle, not the KSP object itself. Do it on a focused branch with full regression run.

### Common rules
- All comm arguments use `petsc4py.PETSc.COMM_WORLD` (never `MPI.COMM_WORLD`). Matches the MPI contract above.
- Positional (`KSP().create(PETSc.COMM_WORLD)`) vs keyword (`SNES().create(comm=PETSc.COMM_WORLD)`) is purely cosmetic; both work identically.

## Flow accumulation (summary)
The IDA system `(I − Wᵀ) q = b` drives discharge and every water/sediment cascade. Its non-convergence is the most common symptom of trouble at scale. Know before touching it: the primary is `fgmres` (stationary Richardson diverges partition-dependently); a failed routing solve is first retried with exact block factors (`_solveIDAExact`, sticky once it works: long drainage chains are well posed, just beyond ILU/Richardson); the main discharge solve is `fatal=True` but **ponds** a small finite un-drained region instead of aborting (cap `GOSPL_UNDRAINED_CAP`); a failed iterate is never accepted; the outer cascade loop is bounded (stagnation break, relative floor, `_cascade_max_steps`); a reachability pin was tried and rejected. All of it, with the reasons: `gospl/flow/AGENTS.md`.

## The Model god-class
`gospl/model.py` declares `Model` (`class Model(`) as multi-inheritance of 17 mixins. `Model.__init__` calls each parent's `__init__` **by name, not via `super()`**. The init order is load-bearing and differs from the MRO declaration order — adding `super().__init__()` will break the chain. Line numbers below are as of 2026-10-03; search for the call if they have drifted.

| # | Init call (model.py line) | Assumes already populated | Allocates / sets |
|---|---|---|---|
| 1 | `_ReadYaml(filename)` :247 | — | `self.input`, every YAML attr, `self.tNow` |
| 2 | `_STRAMesh()` :256 | ReadYaml (`strataFile`) | `stratH/stratZ/phiS/stratK = None` |
| 3 | `_VoroBuild()` :259 | — | Voronoi cache attrs reset |
| 4 | `_UnstMesh()` :262 | ReadYaml + STRAMesh | `dm`, `hLocal`, `hGlobal`, `locIDs`, `glbIDs`, `lpoints`, `mpoints`, `lcoords`, `mCoords`, `larea`, `FVmesh_ngbID`, `lgmap_*`, `idBorders/idLBounds/ghostIDs`, `bL/bG`, calls `readStratLayers` |
| 5 | `_WriteMesh()` :265 | UnstMesh | `step`, `outputDir`, `upsG/upsL` |
| 6 | `_FAMesh()` :268 | UnstMesh | `iMat`, `fillFAL`, `FAG`, `FAL`, `rtol`; defines `_matrix_build`, `_matrix_build_diag`, `_solve_KSP*` for every downstream class |
| 7 | `_IceMesh()` :271 | FAMesh | `iceHL/iceMeltL/iceMeltRiverL/iceUbL/iceAbrL/iceFAL/iceFAG/iceFlex` + cached `iceMat` (only if `iceOn`) |
| 8 | `_GWMesh()` :276 | FAMesh | groundwater head/recharge/duricrust/solute state + cached `_ksp_gw`/`_gwMat` (only if `gwOn`); runs before `_soilSPL`, so `lHbed` is initialised in `soilSPL.__init__` |
| 9 | `_SPL()` :279 | FAMesh | `hOld`, `hOldLocal`, `hOldFlex`, `Eb`, `EbLocal`, `stepED`, `newH` |
| 10 | `_nlSPL()` :282 | SPL | `snes_rtol/atol/maxit`, lazy `_snes_ed/_snes_nl` |
| 11 | `_soilSPL()` :285 | nlSPL | `Gsoil`, `Lsoil`, `lHbed`, `gHbed`, `prodSoil`, `soil_transition` |
| 12 | `_PITFill()` :288 | UnstMesh | `borders`, `outEdges` |
| 13 | `_SEDMesh()` :291 | FAMesh | **`tmp`, `tmpL`, `tmp1`**, `Qs`, `QsL`, `nQs`, `vSed`, `vSedLocal`, `maxnb` |
| 14 | `_hillSLP()` :294 | SEDMesh | **`h`, `hl`, `dh`**, `mat`, `Dlimit/dexp/minDiff`, lazy `_snes_hill`, `_ts_marine` |
| 15 | `_SEAMesh()` :297 | FAMesh | `zMat` |
| 16 | `_GridProcess()` :300 | ReadYaml + UnstMesh | `localFlex`, cached FEM operator/KSP (flat `fem`) or DH grid (`global`); cached orography advection operators if orography on |
| 17 | `_UnstMesh.applyForces(self)` :303 | everything above | `rainVal`, `upsub`, `sealevel` |
| 18 | `_Tectonics()` :306 | UnstMesh (`hGlobal`) | `fiso`, `tecNb=-1` |

Permuting any line silently uses an unallocated attribute or overwrites one another mixin has already set.

## Scratch vector contract (CRITICAL)
The following PETSc Vecs are **scratch**. Any kernel may overwrite them at any time. They are **NOT persistent state** — they exist only so collective allocations are paid once at init.

| Vec | Allocated at | Local/Global |
|---|---|---|
| `self.tmp` | sedplex.py:29 | global |
| `self.tmpL` | sedplex.py:30 | local |
| `self.tmp1` | sedplex.py:31 | global |
| `self.Qs`, `self.QsL`, `self.nQs` | sedplex.py:32-34 | g, l, l |
| `self.h` | hillslope.py:40 | global |
| `self.hl` | hillslope.py:41 | local |
| `self.dh` | hillslope.py:42 | global |
| `self.newH` | SPL.py:40 | global (also held by SNES) |
| `self.stepED` | SPL.py:38 | global |
| `self.upsG`, `self.upsL` | outmesh.py:64-65 | global, local |

**Reading `self.h` for elevation gives you whatever the last hillslope step wrote.** Use `self.hLocal/hGlobal` for elevation, `self.cumED/cumEDLocal` for cumulative ED, `self.FAL/FAG` for flow accumulation, `self.Eb/EbLocal` for erosion rate, `self.vSed/vSedLocal` for sediment volume.

Rule: any new method that uses a scratch Vec MUST document which ones in its docstring and MUST leave them in a defined state on exit (typically `set(0.0)` or a freshly-written array, never mid-computation).

## The rcvID / rcvIDi convention (CRITICAL)
`flowAccumulation` (`flowplex.py:787-791`) snapshots five arrays **after the first `_buildFlowDirection` on the unfilled topography, before pit filling and before downstream-routing rebuilds the flow matrix**:
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

## Opt-in features (summaries)
Each is gated by a flag and must stay **byte-identical when off** (guard that with a test when you touch it).
- **Dual lithology** (`strata: dual`, `stratLith`): coarse/fine sediment; deposition is composition-only, so it never changes elevation. Detail: `gospl/sed/AGENTS.md`, `docs/DESIGN_DUAL_LITHOLOGY.md`.
- **In-model provenance tracers** (`provenance:`, `provOn`): a passive label carried through erosion, transport, deposition and stratigraphy; conservation is machine-exact. Detail: `gospl/sed/AGENTS.md`, `docs/DESIGN_PROVENANCE.md`.
- **Diagnostic glacial model** (`ice:`, `iceOn`): routed discharge → Bahr thickness → sliding → abrasion/till/meltwater. There is no SIA solver; do not reintroduce a serial epsfill for the ice routing surface. Detail: `gospl/flow/AGENTS.md`, `docs/DESIGN_ICE_SHEET.md`.
- **Groundwater, duricrust, Level-B geochemistry** (`groundwater:`, `gwOn`/`duriOn`/`gwGeochemOn`): implicit Dupuit head (hypre), duricrust armoring via `_surfaceLithoK`, conservative multi-species solutes (direct LU). Detail: `docs/DESIGN_WATERTABLE_DURICRUST.md`, `docs/DESIGN_WATERTABLE_GEOCHEM.md`, `docs/tech_guide/groundwater.rst`.
- **Soil / regolith** (`soil:`): `mode: regolith` keeps deposits in the stratigraphy. Detail: `docs/DESIGN_SOIL_REGOLITH.md`.

## Input parsing and forcing (summary)
The `_extra*` parser methods are **mandatory continuations**, not optional parsers. Use `self._get_param(...)` / `dict.get` for defaults (do not reintroduce `try/except KeyError` boilerplate; the 33 remaining ones are deliberate and marked `# TODO-REFACTOR`). Forcing DataFrames are read with **named** `df.at[nb, col]`, never `iloc[nb, k]`. Detail: `gospl/tools/AGENTS.md`.

## Input mesh (summary)
Maximum vertex degree is **12** (fixed-width Fortran neighbour tables). `epsfill` seeds only on `elev < cut-off` and its 1-ULP increment does not survive float32. Detail: `fortran/AGENTS.md`.

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
| `MARINE_SMOOTH_N_LAND` | `1.0` | Dimensionless smoothing strength for emergent land in the marine flow-direction smoother (`hillslope._hillSlope(smooth=2)`). |
| `MARINE_SMOOTH_N_SEA` | `5.0` | Dimensionless smoothing strength for the seafloor in the same smoother (heavier than land). Per-node `Kd = N·cell_area`, so the smoothing is timestep- and resolution-independent (replaces the old `Cd·dt`, `Cd∈{1e5,5e6}` m²/yr, whose effective strength `v≈Cd·dt/Δx²` silently varied with both `dt` and mesh resolution — `~1e-3` on the coarse fixtures, `~55` on a 30 km/`dt=1e4` run). Used only to derive marine flow directions in `seaplex._matOcean`; never alters elevation. |
| `ICE_COVER_MIN` | `1.0e-2` (1 cm) | Ice thickness above which a land cell is treated as ice-covered for the subaerial soil gate (`soilSPL._iceFrozenMask`): pedogenic soil production is suppressed and the regolith is held **frozen inert** (preserved, NOT zeroed — unlike the subaqueous marine/lake case). Matches the ice-presence threshold used inline in `flow/iceplex.py` (`_glacialMeltwater`/`_routeTill`), which is a **distinct role** (kept inline there with a `# TODO-REFACTOR` note). See `DESIGN_SOIL_REGOLITH.md` §3. |
| `ICE_DISCHARGE_REL_FLOOR` | `1.0e-6` | Fraction of the maximum smoothed ice discharge below which `_iceFlowMFD` sets the discharge to zero before the Bahr thickness `H ∝ Q^0.3`. The implicit smoothing leaves an unresolved tail on every cell (below ~100× the solve's `rtol` of 1e-8); `Q^0.3` turned it into a thin ice fringe set by solver noise, so by the partition. Removes only cells with zero routed discharge (0.7% of ice volume on glacial_erosion). |

**Same value, different role — DO NOT replace these with the listed constants:**
- `1.0e-8` at `sed/stratplex.py:219` (thickness numerical-noise floor).
- `1.0e-3` at `flow/flowplex.py:363` (water-routing convergence), `sed/sedplex.py:139` (sediment-routing convergence), `flow/pitfilling.py:568,622` (minh epsilon nudges), `sed/seaplex.py:507,509` (clinoH 1mm offset).
- `1.0e-2` in `flow/iceplex.py` (`_glacialMeltwater`, `_routeTill`) is the ice-presence threshold (m): a distinct physical role, not a tolerance.

Each of these is marked with a permanent `# TODO-REFACTOR: value matches X but distinct role; do not replace` comment so future readers know the coincidence is intentional.

## High-risk modules (do not edit without full regression run)
- **`mesher/unstructuredmesh.py`** — owns `dm`, every shared mesh attribute, forcing dispatch (`applyForces` → `_updateRain`/`_updateIce`/…), and `destroy_DMPlex` (~line 1096, runs to EOF) which names every Vec/Mat by hand. Adding a new persistent Vec elsewhere requires adding it to that destroy list or it leaks (e.g. the dual-lithology `self.vSedF`/`self.vSedFLocal` and the diagnostic-glacial `self.iceMeltL`/`self.iceUbL`/`self.iceAbrL` are registered there).
- **`flow/flowplex.py`** — owns `_solve_KSP`, `_solve_KSP2`, `_matrix_build`, `_matrix_build_diag`, `_buildFlowDirection`. Consumed by every downstream module.
- **`tools/inputparser.py`** — owns all parameter parsing; forcing DataFrame column order is API; the `_extra*` chain is mandatory.

## Debugging (summary)
- **An elevation spike:** follow the `debug-spike` skill. `gospl/tools/zprobe.py` (`-v`) names the stage that raised the maximum; check the `fDep` 0.99 cap first (keep `spl: G ≲ rain/2`). Detail: `gospl/tools/AGENTS.md`.
- **A result that changes with rank count, or a hang at np>1:** follow the `parallel-check` skill (`scripts/lint_mpi_collectives.py`, `scripts/ab_partition.py`).
- **A wandering `SIGABRT` on Linux that macOS does not reproduce:** a Fortran out-of-bounds write. Run `scripts/fcheck.sh`. Detail: `fortran/AGENTS.md`.
- **A mass imbalance:** `scripts/budget.py <output-dir>`, or `--summary run.jsonl` for per-step budgets.
- **The model aborted on a large un-drained region:** first run `scripts/fill_mesh_pits.py --check` on the input mesh. A topography problem is not a solver problem; `GOSPL_FLOW_KSP=richardson` does not fix a singular block.

## Known bugs (fix before refactoring)
All issues found by the 2026-10 analytical benchmarks and partition checks are fixed (see `docs/dev/FIXED_BUGS.md`). When you find one, record it here with how to reproduce it, and pin it with a benchmark or test at the current behaviour.
- _(none currently open)_

## Lessons from fixed bugs (full list: `docs/dev/FIXED_BUGS.md`)
- Any `Vec`/`Mat` reduction, scatter or `Allreduce` under `if MPIrank == 0` or a rank-local `.any()` deadlocks at np>1, and serial always passes.
- A scratch Vec that holds a result you keep must not be reused for an intermediate reduction on the skip path (IIOE2 zeroed fields this way). The converse holds too: a scratch Vec passed as a KSP solution with `guess=True` is a STALE starting guess; `_solve_KSP` now resets one that is worse than zero, but seed routing solves explicitly (`b`, `seed=True`).
- A diffusion on absolute elevation that clamps its increment to one sign is not conservative on relief: rescale the increment to the input volume.
- One YAML key must have one meaning at every site that reads it (`oFill`).
- A fatal solve distinguishes a genuinely broken state (large or non-finite: abort) from a knife-edge local singularity (small and finite: pond and continue).
- A serial rank-0 step inside `Model.__init__` presents as a hang at np>1 if it is slow (the `domain: radius` DH-grid query).
- A linearly-implicit integrator (Rosenbrock `rosw`) needs each stage SOLVED; `ksp preonly` turns it into an inexact, partition-dependent scheme its error estimator cannot see.
- A ghost node's Fortran `FVarea` (and any stencil quantity computed over its truncated local neighbourhood: range, outflow count) is NOT the owner's. Harmless while ghost rows are dropped at assembly; wrong as soon as an owned row reads a ghost's derived value (the IIOE2 `thetain = 1 - thetaout(ghost)` bug). Pass halo-synced inputs.
- `idBorders`, `outletIDs` and `advectBorders` list OWNED nodes only (ghost copies of edge nodes are not in them), while `northPts`/`southPts`/`eastPts`/`westPts` are geometric and include ghosts. Used as an exclusion mask for a neighbour average/min, an owned-only set lets a ghost edge neighbour's raw value in, which is partition-dependent (the advection edge-reset bug). Sync the flag to the ghosts first (see `tectonics._advEdgeHalo`).
- Test the parameter regimes the fixtures do not: every stratigraphy fixture had `G: 0`, which hid a transport-limited double count for months. When a code path branches on a parameter, a guard needs a fixture on each side.
- A regression guard must fail without its fix. Verify that before committing.

## Intentional surprises (do NOT "fix")
- **`gid` argument on `mfdreceivers` / `mfdrcvrs`**. The Fortran kernels take an extra `gid(nb)` argument (per-local-node vertex ID `self.gid`, set in `mesher/unstructuredmesh.py`; since 2026-06-21 it is `self.locIDs`, the partition-INVARIANT input-mesh id, no longer PETSc's per-partition numbering; see `gospl/flow/AGENTS.md`). Inside the kernels, the stored slope is perturbed by `val * (1 - 1.0e-15 * gid(n))` before the quicksort tie-break. This makes EXACT slope ties resolve deterministically across MPI decompositions (lower global ID wins). The perturbation is well below KSP-solver precision and does not affect physical results; near-ties driven by KSP floating-point noise are a separate problem and are NOT fixed by this pattern. See the relaxation comment on `test_parallel_correctness` in `tests/test_parallel.py` and `fortran/functions.F90:mfdreceivers` for the full rationale. **Do NOT remove the `gid` argument** — without it, `mfdreceivers` ordering depends on local iteration order and the parallel-correctness test diverges by ~10x from current baseline.
- **Loose tolerances in `test_parallel_correctness`** (`rel_sum_fa` 5e-2, `rel_mean_h` 5e-9): the platform-dependent KSP floating-point floor, not a regression. Tightening them requires bitwise-identical halo state across decompositions. Detail: `docs/dev/FIXED_BUGS.md` > Test tolerance notes.
- **`filllabel`** is uncalled but kept on purpose; do not remove it in a dead-code sweep (`fortran/AGENTS.md`).
- **Fortran `meshparams` allocatables** must be freed per mesh in `definetin`; see `fortran/AGENTS.md`.

## CI contract (summary)
- `tests-pr.yml`: `pytest tests/` + the MPI lint on PRs and pushes to `master`/`release-candidate`. **Pushes to `dev` run no CI**; the nightly cron (`tests-slow.yml`, `examples-smoke.yml`) runs from the default branch, which is `dev`. So `dev` is tested nightly, not per push: run the suite locally before pushing.
- `tests-slow.yml`: slow tier nightly/tags/dispatch, benchmarks on `master`/`release-candidate` pushes. `examples-smoke.yml`: every local goSPL-examples input for 2 steps at np=2 (`scripts/examples_smoke.py`), nightly + dispatch (`only: Global-examples` for the heavy ones).
- On a `v*` tag, `conda-build`, `pypi-publish` and `docker-build` each `needs:` the `_release-gate.yml` tests job, so nothing publishes unless the suite passes.
- Matrix: `ubuntu-latest + macos-15 × Python 3.11 + 3.12`. macOS is pinned (not `macos-latest`) because the osx-arm64 OpenMPI 4.x / petsc4py 3.21 stack is validated per image.
- `environment.yml` conventions, cache-path pitfalls, MPI pins: `.github/AGENTS.md`.

## Docs / Read the Docs (autodoc)
The API reference (`docs/api_ref/*.rst`) is Sphinx **autodoc** — it imports each module to extract docstrings, so the modules MUST import under the docs build. Two invariants keep the API pages from rendering empty (both regressed once and produced blank pages for *most* classes):
- **`gospl` must be importable.** It is NOT pip-installed on Read the Docs, so `docs/conf.py` puts BOTH the repo root (`..`, so `from gospl.tools.constants import …` and other `from gospl.X import` resolve) AND `../gospl/` (so the legacy top-level `.. autoclass:: sed.hillslope.hillSLP` paths resolve) on `sys.path`. Don't remove either.
- **Mock only the genuinely-absent compiled/MPI deps via `autodoc_mock_imports`** (`h5py`, `mpi4py`, `petsc4py`, `vtk`, `pyshtools`, `gospl._fortran`). Do NOT mock packages `docs/requirements.txt` installs (`numpy`, `scipy`, `pandas`, `numpy-indexed`, `ruamel.yaml`) — the old manual `sys.modules[m] = Mock()` shadowed the real `scipy` and broke `from scipy.special import …` (and didn't cover submodules), blanking the pages. `autodoc_mock_imports` mocks submodules automatically; the `READTHEDOCS` env-guard around `from gospl._fortran import …` is the belt-and-braces for the compiled extension.
- When you add a public/private method that should appear in the API, add it to BOTH the `.. autosummary::` and `.. automethod::` lists in the relevant `docs/api_ref/*_ref.rst` (they are hand-maintained, not auto-generated).

## Analytical benchmark suite (summary)
`benchmarks/` validates physics against exact solutions (`pytest benchmarks/ -m benchmark`; scipy/matplotlib, skipped when absent). Mark new ones `@pytest.mark.benchmark @pytest.mark.slow`, build the mesh at runtime in `tmp_path` (see `benchmarks/test_dupuit.py`), and keep each under ~2 minutes. Conventions and per-benchmark notes: `benchmarks/AGENTS.md`.

| Test file | Process | Analytical basis | Pass criteria |
|---|---|---|---|
| `test_spl.py` | SPL steady state | Braun & Willett 2013; Perron & Royden 2013 | Cases 150 & 200: 5/5 each (basin test skipped — no case meets >80%); case 100 excluded (finest-mesh basin/R² artifacts + ~34 min cost) |
| `test_hillslope.py` | Hillslope diffusion | `z=(U/2κ)x(L-x)`; Roering, Kirchner & Dietrich 1999 | All `TOL_*` constants met |
| `test_knickpoint.py` | Knickpoint propagation | `c=K·A^m`; Royden & Perron 2013; Tucker & Whipple 2002 | 4/4 sub-tests |
| `test_dupuit.py` | Groundwater head (Dupuit-Boussinesq) | `s²=s_d²+(R/K)(2Wx−x²)` | RMSE 0.00 %, R²=1 |
| `test_flexure.py` | Flat FEM + global flexure | cosine modes `w=q/(Dk⁴+Δρg)` (exact for the 5-point stencil); infinite-plate line load with image loads; degree-l spherical-harmonic response `q/(Δρg+D·P_l)` + flat-plate limit | FEM: discrete round-off, continuum <1%, 2nd-order convergence; line load RMSE <0.5%; global spectral 1e-9, mesh pipeline RMSE <1% |
| `test_orography.py` | Orographic rain (Smith-Barstad, no mountain wave) | `P̂ = Cw·iku·ĥ/((1+ikuτc)(1+ikuτf))` on a 1-D ridge, with the background/floor clip reproduced | central-row RMSE <1.5% of peak, wind-parallel edge rows <1.5% (inflow-only BC since 2026-10); first-order convergence |
| `test_advection.py` | Horizontal advection (upwind / iioe1 / iioe2) | translated Gaussian; uniform field (iioe2 zero-excess path) | mass < 1e-4 for every scheme (iioe2 mass-neutral since 2026-10), centroid, peak, RMSE; uniform field preserved |
| `test_marine_diffusion.py` | Marine diffusion (`ts` and `picard`) | 2-D heat kernel on a constant-`Cd` deposit | volume 1e-12; D err, peak law, RMSE < 2e-3 for both solvers (the `ts` stage-solve fix, 2026-10) |
| `test_geochem_balance.py` | Level-B geochemistry | per-step `dissolved = precipitated + exported`; 1-D strip `c=D0/(R+p)·(Ls/u)^(1+p/R)` | closure ~1e-15; strip vs exact discrete 2e-5 |

## Checklist before any commit
1. Did you read this file? If invariants here changed, update them.
2. Did you run the full regression suite (`pytest tests/ -n 4`, or at least the markers you touched plus `-m mpi` if you touched a collective, a solver or an assembly)? Did a new guard fail without your fix?
3. If you added a column to a forcing DataFrame, did you use named access (df.at[nb, col]) not iloc?
4. If you used a scratch Vec, did you document which ones in the method's docstring?
5. If you changed a method called by `Model.runProcesses` (`model.py`, `def runProcesses`), did you check every caller AND the mixin init order (`Model.__init__`)?
6. If you added a new KSP/SNES/TS solver, did you pick the right lifecycle (CACHED for hot-path solvers, AD-HOC for nested-fieldsplit) AND, if cached, add the attribute to the `destroy_DMPlex` list in `unstructuredmesh.py`?
6b. If you added or moved a **collective** call (`localToGlobal`/`Allreduce`/`ksp.solve`/`vec.norm`/`vec.sum`/`garbage_cleanup`/…), is it reached by ALL ranks — NOT inside `if MPIrank==0`, a rank-local `.any()`/`.size`/early-`return`, or a loop with a per-rank trip count? If a rank-local decision must gate it, reduce the flag first (`allreduce(..., op=MPI.LOR)`). Run `python scripts/lint_mpi_collectives.py` and `pytest tests/ -m mpi`.
7. If you added a new forcing type, did you follow `docs/HOW_TO_ADD_FORCING.md` including the `destroy_DMPlex` registration?
8. If you added a new output field, did you follow `docs/HOW_TO_ADD_OUTPUT.md` including the `destroy_DMPlex` registration and XDMF entry?
9. If you added a benchmark test, did you apply `pytest.importorskip` for scipy/matplotlib AND mark `@pytest.mark.benchmark @pytest.mark.slow`?
10. If you added a goSPL `Model` init in a test (regression OR benchmark), is `model.destroy()` in a `try/finally` block per AGENTS.md > KSP/SNES/TS lifecycle contract?
11. If you bumped the version, did you change all four locations (`meson.build:4`, `conda/meta.yaml:2`, `docs/conf.py`, `docs/_static/version_switch.json`) in the same commit, with the no-leading-zero spelling (`2026.6.13`)? See `docs/dev/RELEASE.md`.
12. If you changed a Fortran kernel, did you edit `fortran/functions.F90` AND `fortran/functions.pyf` together, bound-check every write into a fixed-width `(npoints, 12)` table, and run `scripts/fcheck.sh`? An out-of-bounds write passes on macOS and surfaces as a wandering `SIGABRT` on glibc. See `fortran/AGENTS.md`.
13. Did you add a row to `docs/dev/CHANGELOG_DEV.md` for a feature or fix (and an entry in `docs/dev/FIXED_BUGS.md` for a bug fix, naming its guard test)?
