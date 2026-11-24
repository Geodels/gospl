# AGENTS.md

Rules for anyone (human or AI agent) changing goSPL. Read this file at the start of
every session, and the subsystem file in the **Map** below before you edit that part
of the tree. Update the relevant file when an invariant changes. Long-form rationale:
`REFACTOR_AUDIT.md`; past bugs and their guards: `docs/dev/FIXED_BUGS.md`; per-change
history: `docs/dev/CHANGELOG_DEV.md`.

Last reviewed 2026-11-23 on `feat/gravity-currents` (gravity currents phases G1-G3, not merged into `dev`; after merging `dev` @ `01ca0d2`, the base-model partition fixes). In October 2026 this file
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
| Quick check (255 of 278 tests, ~25 s) | `pytest tests/ -n 4 -m "not mpi"` |
| One subsystem | `pytest tests/ -m flow` (markers: `tests/README.md`) |
| Multi-rank guards (~70 s) | `pytest tests/ -m mpi` |
| Full suite before a commit | `pytest tests/ -n 4` (CI runs it serially) |
| Analytical benchmarks (minutes) | `pytest benchmarks/ -m benchmark` |
| Collective-under-rank-guard lint | `python scripts/lint_mpi_collectives.py` (CI runs it) |
| np=1 vs np=N field comparison | `python scripts/ab_partition.py <input.yml> -n 4` |
| Mass budget of a run | `python scripts/budget.py <output-dir>` |
| Wall time per sub-stage (marine smoother, epsfill, cascades, sediment load, SPL), rank 0 vs max rank | `mpirun -n N python scripts/stage_profile.py <input.yml>` |
| Per-rank memory by setup stage / step phase (rank 0 vs others) | `mpirun -n N python scripts/mem_profile.py <input.yml> [--json out.json]` |
| Machine-readable per-step summary | `gospl -i input.yml --summary run.jsonl` |
| Inspect / compare output | `gospl-inspect <output-dir> [--step N] [--compare <other-dir>] [--json]` (or `python -m gospl.analyse.runinspect` before a `pip install -e .` registers the command) |
| Fortran bounds-checked suite | `scripts/fcheck.sh` |
| goSPL-examples smoke run | `python scripts/examples_smoke.py <examples-dir> --only Local-examples` |
| CI runs / logs | `gh run list`, `gh run view <id> --log-failed` (`gh` is authenticated) |

## What goSPL does
goSPL is a parallel landscape-evolution model that integrates the stream-power law (river incision), linear and non-linear hillslope diffusion, marine sediment transport, glacial accumulation, flexural isostasy, and horizontal/vertical tectonics on an unstructured Voronoi/Delaunay finite-volume mesh. The mesh is either a 2D flat plane (`self.flatModel == True`) or a global sphere; partitioning, halo exchange, and all linear/non-linear solves run on PETSc DMPlex via petsc4py. Time integration is an explicit outer Euler loop in `Model.runProcesses` with implicit KSP/SNES/TS inner solves for diffusion, flow accumulation, and sediment routing.

## The numpy ↔ PETSc boundary
Every state field exists in two parallel representations.

**Numpy land** (raw arrays indexed by local node ID, dimensional, no halo): `self.lcoords` (m), `self.mCoords` (m; whole input mesh, a lazily loaded property: see below), `self.larea` (m²), `self.rainVal` (m/yr), `self.upsub` (m/yr), `self.stratH/stratZ/phiS/stratK`, `self.stratHf/phiF` (dual-lithology fine pile — see `## Dual lithology`), `self.fineFrac/depoFineFrac` (dual), `self.pitParams`, `self.pitIDs`, `self.lFill`, `self.localFlex`, plus any `vec.getArray().copy()` view.

**PETSc land** (parallel Vec with halo, mutated via `setArray`/`getArray`/`localToGlobal`): `self.hLocal`/`self.hGlobal`, `self.cumED`/`self.cumEDLocal`, `self.FAL`/`self.FAG`, `self.fillFAL`, `self.Eb`/`self.EbLocal`, `self.bL`/`self.bG`, `self.areaLocal`/`self.areaGlobal`, `self.iceHL`/`self.iceMeltL`/`self.iceUbL`/`self.iceAbrL`/`self.iceFlex` (diagnostic glacial model — see `## Ice sheet`), `self.Lsoil`/`self.Gsoil`, `self.lHbed`/`self.gHbed`, `self.vSed`/`self.vSedLocal`, `self.vSedF`/`self.vSedFLocal` (dual-lithology fine flux), `self.fiso`.

Both sides hold physical units; the boundary is about **who owns halo synchronisation**, not units. Cross only via `self.dm.localToGlobal(local, global)`, `self.dm.globalToLocal(global, local)`, `vec.getArray()`, `vec.setArray(arr)`. After mutating a `*Local` array view, you MUST `localToGlobal` before the next collective solve, or ranks see stale halos.

**No whole-mesh arrays on every rank.** Memory per rank must scale with the partition (`lpoints`), not the mesh (`mpoints`): at 192 ranks a per-rank whole-mesh copy costs more than the rank's own data, and 48 ranks share a node. What holds this (2026-11):
- **Mesh build (`_buildMesh`).** Only rank 0 reads the input mesh. Every rank gets its vertices' input ids (`locIDs`, exact: `distributeField` of an index field, no KDTree) and elevations through the DMPlex distribution; `glbIDs == locIDs[glIDs]` (checked collectively). `self.mCoords` is a property that loads the whole-mesh coordinates on first use: only the global-flexure DH grid (rank 0) and plate advection read it, and the latter releases it. `VoroBuild.releaseVoronoi()` frees the Voronoi intermediates after `_meshStructure`.
- **Forcing maps** (`rainMesh`, `evapMesh`, `sedFacMesh`, `elaMesh`/`iceMesh`/`termMesh`) are local slices (`[locIDs]` at load). Follow `docs/HOW_TO_ADD_FORCING.md` for new ones.
- **Plate interpolation** (`advect: interp`) builds the KDTree on the whole advected mesh once per plate event, then drops it; it keeps only local `tec_*` rows (local `tec_onIDs`) plus a cached owner-fetch plan (`_tecNeeded`/`_tecRows`/`_tecWant`). With `interp_remap` the barycentric tables read the cell list once as int32 and keep only the cells touching the displaced vertices near the rank's nodes (no whole-mesh centroid tree). `_advectPlates` fetches the surface fields at the source nodes from their owners (no whole-mesh `Allreduce`) and scatters owned rows with `[glIDs]`; `_advectStrati` fetches the layer rows the same way (`_fetchOwnerRows`: memory-scalable AO + one blocked scatter) and applies `stratalRecord`.
- **Advection schemes** (`upwind`/`iioe1`/`iioe2`) move the stratigraphy semi-Lagrangian (`_advectStratiScheme`: departure point `x - v*dt`, tables built once per tectonic event, then `_fetchOwnerRows` + `stratalRecord` every step; ~3% of run time with a layer every step on 1M nodes). On a cyclic model plate interpolation and the departure points displace nodes exactly around the cylinder (`_cylinderDisplace`).
- **Serial rank-0 steps** (marine and gravity depression fills, global flexure) gather with `_gatherGlobalOnRoot` and return the result with `_scatterFromRoot` (one `Scatterv` of each rank's `locIDs` rows); never `bcast(full)[self.locIDs]`, which put an `mpoints` array on every rank (2026-11-23).
- **Guards:** `test_scatter_from_root`, `test_mesh_not_replicated_per_rank`, `test_plate_interp_tables_local`, `test_plate_strata_partition_invariant` (interp + iioe2), `test_scheme_advection_moves_stratigraphy`, the cyclic/sphere advection tests.

**Memory benchmark** (`scripts/mem_profile.py`, flat synthetic meshes, macOS, peak RSS in MB, `dev` @ `3e17a13` vs this branch):

| mesh, ranks | other ranks, init | other ranks, run | rank 0, run |
|---|---|---|---|
| 1M, np=4 | 1475 → 1414 | 2754 → 2494 (−9%) | 3334 → 3150 (−6%) |
| 2M, np=16 | 1288 → 1096 (−15%) | 1670 → 1527 (−9%) | 3839 → 3822 |

Read it as: the import baseline is ~0.55 GB/rank; the per-rank local cost is ~5–8 kB per local node; the removed whole-mesh part was ~100 B per MESH node per rank at the setup peak (+~30 B held all run), so the saving grows with mesh size and rank count. RSS rarely drops after a free, so attribute memory by per-stage growth, not absolute values. **Rank 0 is now the hotspot:** it builds the serial DMPlex (`createFromCellList` + `distribute`, ~1.8 kB per mesh node; 3.8 GB at 2M) and holds `FVgnID` (`epsfill`) and the pit-spillover graph. Two ways to shrink it were measured:
- *Parallel construction* via PETSc's HDF5 load (`dm.load` of a `dmplex_storage_version` 3.0.0 file, each rank reading a chunk, then `interpolate()` and a ParMETIS repartition; petsc4py 3.21 does not wrap `DMPlexCreateFromCellListParallelPetsc`): it loads correctly, but every rank pays ~190 B per mesh node during the load (4M nodes: 3.9/2.1/1.5 GB per rank at np=2/4/8 vs 2.3 GB on rank 0 alone serially), i.e. worse per node, and the local conda PETSc cannot ParMETIS-repartition a distributed graph (fails even on a known-good mesh; PT-Scotch shim). Not adopted.
- *Uninterpolated serial build + parallel interpolate* (`createFromCellList(interpolate=False)`, distribute, then `dm.interpolate()` per partition; vertex points and coordinate order are preserved): rank 0 run peak 3822 → 2881 MB (−25%) at 2M/np=16, others unchanged, partition quality equal (balance, ghost counts). **Rejected (2026-11-11, maintainer decision):** ParMETIS gets a different cell graph, so on irregular meshes the partition differs from `dev` and np>1 results are no longer bitwise equal to earlier runs (np=1 is); bitwise reproducibility with existing runs was preferred over the rank-0 saving. Do not re-propose without that trade-off.


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
| `eroder/soilSPL.py` | `diffuseSoil` | `_ts_soil` | TS rosw; stage solve gmres+gasm (NOT preonly), atol 1e-3 m, rtol 0 |
| `sed/hillslope.py` | `_hillSlopeNL` | `_snes_hill` (+`_snes_hill_sG/sL` slope halo) | SNES `qn` (`hillnl_`; bare `ngmres` aborted on fine meshes) |
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
The IDA system `(I − Wᵀ) q = b` drives discharge and every water/sediment cascade. Its non-convergence is the most common symptom of trouble at scale. Know before touching it: the primary is `fgmres` (stationary Richardson diverges partition-dependently); a failed routing solve is first retried with exact block factors (`_solveIDAExact`, sticky once it works: long drainage chains are well posed, just beyond ILU/Richardson); the main discharge solve is `fatal=True` but **ponds** a small finite un-drained region instead of aborting (cap `GOSPL_UNDRAINED_CAP`); a failed iterate is never accepted; the outer cascade loop is bounded (stagnation break, relative floor, `_cascade_max_steps`); a reachability pin was tried and rejected. All of it, with the reasons: `gospl/flow/AGENTS.md`. Since 2026-11-12: routing reads `lFill` for cells within 1 nm of it (`FAMesh._routingSurface`; the shifted-frame fill leaves 1-ulp noise between `h` and `lFill`, and flat routing on `lFill` vs MFD on `h` formed 2-cycles; clamping `lFill >= h` instead split DEM plateaus into ~1e5 isolated flat sinks, do not); a routing solve's primary is accepted only on its true residual (`_trueResidualOK`, fgmres reported convergence on a 2-cycle with x = 2e13); and a non-fatal routing solve that fails on a small region is repaired, not zeroed (`_pondRoutingRegion`: those cells become sinks, the rest is re-solved, the trapped volume is deposited), with a `[flow] routing repair region` diagnostic listing the cells.

## The Model god-class
`gospl/model.py` declares `Model` (`class Model(`) as multi-inheritance of 17 mixins. `Model.__init__` calls each parent's `__init__` **by name, not via `super()`**. The init order is load-bearing and differs from the MRO declaration order — adding `super().__init__()` will break the chain. Line numbers below are as of 2026-10-03; search for the call if they have drifted.

| # | Init call (model.py line) | Assumes already populated | Allocates / sets |
|---|---|---|---|
| 1 | `_ReadYaml(filename)` :247 | — | `self.input`, every YAML attr, `self.tNow` |
| 2 | `_STRAMesh()` :256 | ReadYaml (`strataFile`) | `stratH/stratZ/phiS/stratK = None` |
| 3 | `_VoroBuild()` :259 | — | Voronoi cache attrs reset |
| 4 | `_UnstMesh()` :262 | ReadYaml + STRAMesh | `dm` (built on rank 0, distributed), `hLocal`, `hGlobal`, `locIDs`/`glbIDs` and the initial elevation (delivered by `distributeField`, not a KDTree), `lpoints`, `mpoints` (bcast), `lcoords`, `larea`, `FVmesh_ngbID`, `lgmap_*`, `idBorders/idLBounds/ghostIDs`, `bL/bG`, calls `readStratLayers`. `mCoords` is a lazy property (whole mesh, loaded on first use) |
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
- **Channel cells out of the hillslope diffusion** (`diffusion: channelQ`, m³/yr): cells with discharge ≥ channelQ are identity rows of the diffusion and their hillslope supply joins the same step's routed load (`hillslopeBeforeRouting` moves the hillslope step between incision and routing when deposition is on). Removes the resolution dependence of the channel steepness. Detail: `gospl/sed/AGENTS.md`, `docs/tech_guide/hill.rst`; example `goSPL-examples/Local-examples/hillslope_channels`.
- **Marine gravity currents** (`gravity:`, `gravOn`, phases G1-G3): part of the river load reaching the sea (fixed or flood fraction) is routed down the marine network with a run-out length, `(I - W^T diag(1-Lambda)) q = s`, and deposited beyond the clinoforms; coarse and fine have their own run-out under dual lithology; provenance classes ride the same solve (phase G2). Phase G3b: own network on the eps-filled raw seafloor (`network: raw`, default) and optional linear slope erosion (`erosion_length`, erode-or-deposit per node, `max_erosion`). Detail: `gospl/sed/AGENTS.md`, `docs/DESIGN_GRAVITY_CURRENTS.md`. Channel belts and lobes at regional resolution: separate coarse/fine networks (`receivers_fine`, `flow_exp_fine`), smoothing (`network_smooth`, `slope_smooth`), conservative deposit spreading (`deposit_spread`, default two cell sizes), `coarse_factor` (default 0.3). Phase G3: slope-failure trigger (`failure_slope`). Slope erosion is an expert option (warning on cells coarser than `GRAVITY_EROSION_MAX_CELL`). Benchmark `benchmarks/test_gravity.py`. Branch `feat/gravity-currents`, not merged into `dev`.
- **Plate interpolation remaps** (`domain: interp_remap`, mesh spacings, `advect: interp` only): each event is applied in remaps of about that many spacings (whole number of steps set by the fastest node), with barycentric weights in the displaced triangles (`tectonics._barycentricTables`, IDW only where no triangle covers a node). Plain `interp` keeps IDW. Detail: `docs/tech_guide/tecto.rst` (lateral-motion test table).

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
| `DEPOSIT_FLOOR` | `1.0e-9` (1 nm) | Round-off floor of the routed marine deposit before the diffusion. Was 1 mm until 2026-11-25: it discarded real sub-mm deposits without putting the volume back (4e-5 of the activity on the closed test sphere, ~5 % of thinly spread coastal-erosion material). |
| `BEDROCK_EXPOSED` | `1.0e-1` (10 cm) | Soil-vs-bedrock threshold in soil-aware SPL. |
| `BEDROCK_SENTINEL` | `1.0e6` | Infinite-bedrock layer-0 sentinel thickness; offset cancels in cumsum arithmetic. |
| `BOUNDARY_FLOW_SENTINEL` | `-1.0e6` | Mfd no-data marker passed to `mfdreceivers`/`mfdrcvrs`; **distinct from `MISSING_DATA_SENTINEL`** — kept two orders apart so a mix-up is obvious in diagnostics. |
| `GRAPH_OUTLIER_CAP` | `1.0e7` | Upper-bound clamp in `pitfilling._performFilling`; entries above this are rewritten to `MISSING_DATA_SENTINEL`. |
| `MARINE_SMOOTH_N_LAND` | `1.0` | Dimensionless smoothing strength for emergent land in the marine flow-direction smoother (`hillslope._hillSlope(smooth=2)`). |
| `MARINE_SMOOTH_N_SEA` | `5.0` | Dimensionless smoothing strength for the seafloor in the same smoother (heavier than land). Per-node `Kd = N·cell_area`, so the smoothing is timestep- and resolution-independent (replaces the old `Cd·dt`, `Cd∈{1e5,5e6}` m²/yr, whose effective strength `v≈Cd·dt/Δx²` silently varied with both `dt` and mesh resolution — `~1e-3` on the coarse fixtures, `~55` on a 30 km/`dt=1e4` run). Used only to derive marine flow directions in `seaplex._matOcean`; never alters elevation. |
| `ICE_COVER_MIN` | `1.0e-2` (1 cm) | Ice thickness above which a land cell is treated as ice-covered for the subaerial soil gate (`soilSPL._iceFrozenMask`): pedogenic soil production is suppressed and the regolith is held **frozen inert** (preserved, NOT zeroed — unlike the subaqueous marine/lake case). Matches the ice-presence threshold used inline in `flow/iceplex.py` (`_glacialMeltwater`/`_routeTill`), which is a **distinct role** (kept inline there with a `# TODO-REFACTOR` note). See `DESIGN_SOIL_REGOLITH.md` §3. |
| `ICE_DISCHARGE_REL_FLOOR` | `1.0e-6` | Fraction of the maximum smoothed ice discharge below which `_iceFlowMFD` sets the discharge to zero before the Bahr thickness `H ∝ Q^0.3`. The implicit smoothing leaves an unresolved tail on every cell (below ~100× the solve's `rtol` of 1e-8); `Q^0.3` turned it into a thin ice fringe set by solver noise, so by the partition. Removes only cells with zero routed discharge (0.7% of ice volume on glacial_erosion). |
| `GRAVITY_EROSION_MAX_CELL` | `200.0` (m) | Mean cell size above which goSPL warns that gravity-current slope erosion (`gravity: erosion_length`) is not meaningful: submarine channels need several cells across. An estimate (615 m tested, finer not), used only for the warning in `SEAMesh.__init__`. |
| `FILL_LEVEL_TOL` | `1.0e-9` (1 nm) | Distance (m) within which a cell is AT a depression-filled surface (`lFill`, `waterFilled`): routing reads `lFill` there (`FAMesh._routingSurface`), an overflow pass does not treat it as under water, the end of `flowAccumulation` counts it in the lake, and the groundwater recharge/seepage/subaerial masks do not treat it as lake (`gwplex._underLake`). The fill leaves 1-ulp noise between the surfaces, and that noise depends on the partition: a pit's spill node kept or lost its whole discharge on it (2026-11-24). |

**Same value, different role — DO NOT replace these with the listed constants:**
- `1.0e-8` at `sed/stratplex.py:219` (thickness numerical-noise floor).
- `1.0e-3` at `flow/flowplex.py:363` (water-routing convergence), `sed/sedplex.py:139` (sediment-routing convergence), `flow/pitfilling.py:568,622` (minh epsilon nudges), `sed/seaplex.py:507,509` (clinoH 1mm offset): they matched the former 1 mm `DEPOSIT_FLOOR` and are distinct roles.
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
- **A job killed for memory (OOM, `SIGKILL`, "exceeded memory allocation"):** run `scripts/mem_profile.py` at two rank counts; a stage whose peak on the non-root ranks does not shrink with more ranks holds a whole-mesh array per rank (see the memory rule above). On Gadi request full-node memory (`docs/getting_started/installHPC.rst`).
- **`[flow] routing repair …` or `non-fatal routing solve did not converge` in the log:** a routing system is singular somewhere. The `[flow] routing repair region` block lists the cells (input-mesh id, lon/lat, z, fill, sedfill, pit, sink, b, x, upstream donors); look for two cells that are each other's donors (a cycle) and check why the routing rules disagree there before touching the solver (that is how the 1-ulp fill bug was found).
- **The model aborted on a large un-drained region:** first run `scripts/fill_mesh_pits.py --check` on the input mesh. A topography problem is not a solver problem; `GOSPL_FLOW_KSP=richardson` does not fix a singular block.

## Known bugs (fix before refactoring)
All issues found by the 2026-10 analytical benchmarks and partition checks are fixed (see `docs/dev/FIXED_BUGS.md`). When you find one, record it here with how to reproduce it, and pin it with a benchmark or test at the current behaviour.
- _(none currently open)_

## Known limitations (documented, not bugs)
Stated in the GMD model-description paper (revision of Nov 2026) and in the user and technical guides. Do not "fix" them silently: a change of behaviour must update the docs and be discussed with the maintainer.
- **Hillslope diffusion acts on channel cells**, so the hillslope/fluvial balance depends on resolution: with `hillslopeKa: 0.2` the steady channel steepness of an uplifted block is 30.5 / 41.6 / 66.4 at 200 / 100 / 50 m spacing instead of U/K = 25 (25.1-25.3 without diffusion). **Addressed (2026-11) by `diffusion: channelQ`** (m³/yr): cells whose discharge reaches it are identity rows of the diffusion operator (`hillslope._channelMask`, cached operator rebuilt when the mask changes, collective LOR decision like the coastline), and `_channelDiffusion` hands what the hillslopes deliver to them (`d = hOld - (A h)` with the unmasked rows `_hillCoeffs`, or `dt * flux` for the non-linear laws) to the same step's routed source through `_hillDelivG` (the hillslopes then run between the river incision and the sediment routing, `hillslopeBeforeRouting`, which restores the river step's `thCoarse`/`thFine`/`provEro`) (consumed in `sedplex._getSedFlux`); a negative `d` is taken from the channel cell, so the volume budget is exact to the solve tolerance and nothing is left pending at the end of a run (no one-step lag since 2026-11-19: the supply bypasses only the SPL in-river deposition with `G > 0`). k_s 25.1/25.2/25.4 at 200/100/50 m. With stratigraphy, `_channelComposition` (after the hillslope `erodeStrat`) converts the bulk supply to solid with the donors' eroded solid fraction and splits it into fine (`_hillDelivFine`) and per-class (`_hillDelivProv`) parts for the fine and provenance routing. With the soil model (`soilSPL._soilChannelSupply`): channel rows fixed in the TS residual and Jacobian (`_soilChan`); the total supply is exact (volume lost by non-channel cells), its split among channel cells uses the final-state inflow. The soil TS itself was then fixed to solve each Rosenbrock stage (see FIXED_BUGS). Guards: `test_channel_cells_excluded_from_diffusion[linear|series]`, `test_channel_cells_reach_the_river_load`, `test_channel_q_option_checks`, `test_channel_cells_dual_fine_conservation`, `test_channel_cells_provenance_conservation`, `test_channel_cells_total_sediment_conserved`, `test_channel_hillslope_partition_invariant[linear|series]` (np=3 with deposition and `G > 0`; its threshold is off the regular mesh's discharge lattice, since an exact tie flips cells by round-off). Without the option, users must choose `hillslopeKa` for their resolution.
- **Eulerian advection (`upwind`/`iioe1`/`iioe2`) does not carry the drainage across a lateral motion**: the valley pattern lags and then reorganises (Hergarten & Fraters 2026). **Addressed (2026-11) by `domain: interp_remap`**: `advect: interp` applies each event's motion in remaps of about that many mesh spacings (`tectonics._remapInterval`/`_nextRemap`; tables built once per event and reused, rebuilt only for a shorter last or first remap). Lateral-motion test (40 x 20 km at 200 m, 4 mm/yr = 0.1 cell/step, 1 Myr; table in `docs/tech_guide/tecto.rst`): barycentric remaps of 0.5/1/2/5 cells move the valleys 4.0 of 4 km; 0.1-cell remaps move them at the right speed but smooth them away by 1 Myr (IDW 0.1-cell remaps moved them only 0.7 km: IDW copies the nearest node). With the option the interpolation weights are ALWAYS barycentric (no switch) in the displaced triangle containing each node (`tectonics._barycentricTables`, IDW kept where none does; candidate cells only those touching the 12 displaced vertices nearest a local node, so memory scales with the partition), which moves slow regions of a spatially varying velocity at their own speed (IDW is biased towards the nearest node at sub-cell displacements); plain `interp` keeps IDW (byte-identical). Guards: `test_interp_remap_carries_the_surface` (bitwise equal to hand-split events), `test_interp_remap_variable_velocity` (shear, slow bump 80.9 m of 80 m; IDW 70 m). The Eulerian schemes themselves are unchanged.
- **Knickpoint benchmark front**: the first-order implicit scheme smears a knickpoint, so its 5 % and 95 % levels lead and lag the mid-point by ~12 %; only the mid-point follows `∫K A^m dt` (slope 0.996-1.000). This is discretisation, not a celerity error.

The paper's scripts and review tests live outside the repository (paper folder, `review_tests/`). When a change alters a number quoted in the docs or the paper (benchmark results, example outputs), update both.

## Lessons from fixed bugs (full list: `docs/dev/FIXED_BUGS.md`)
- Any `Vec`/`Mat` reduction, scatter or `Allreduce` under `if MPIrank == 0` or a rank-local `.any()` deadlocks at np>1, and serial always passes.
- A scratch Vec that holds a result you keep must not be reused for an intermediate reduction on the skip path (IIOE2 zeroed fields this way). The converse holds too: a scratch Vec passed as a KSP solution with `guess=True` is a STALE starting guess; `_solve_KSP` now resets one that is worse than zero, but seed routing solves explicitly (`b`, `seed=True`).
- A diffusion on absolute elevation that clamps its increment to one sign is not conservative on relief: rescale the increment to the input volume.
- One YAML key must have one meaning at every site that reads it (`oFill`).
- A fatal solve distinguishes a genuinely broken state (large or non-finite: abort) from a knife-edge local singularity (small and finite: pond and continue).
- A serial rank-0 step inside `Model.__init__` presents as a hang at np>1 if it is slow (the `domain: radius` DH-grid query).
- A linearly-implicit integrator (Rosenbrock `rosw`) needs each stage SOLVED; `ksp preonly` turns it into an inexact, partition-dependent scheme its error estimator cannot see. When fixing one solver, check every other solver of the same kind: the soil TS kept `preonly` for a year after the marine TS was fixed.
- A ghost node's Fortran `FVarea` (and any stencil quantity computed over its truncated local neighbourhood: range, outflow count) is NOT the owner's. Harmless while ghost rows are dropped at assembly; wrong as soon as an owned row reads a ghost's derived value (the IIOE2 `thetain = 1 - thetaout(ghost)` bug). Pass halo-synced inputs.
- `idBorders`, `outletIDs` and `advectBorders` list OWNED nodes only (ghost copies of edge nodes are not in them), while `northPts`/`southPts`/`eastPts`/`westPts` are geometric and include ghosts. Used as an exclusion mask for a neighbour average/min, an owned-only set lets a ghost edge neighbour's raw value in, which is partition-dependent (the advection edge-reset bug). Sync the flag to the ghosts first (see `tectonics._advEdgeHalo`).
- Test the parameter regimes the fixtures do not: every stratigraphy fixture had `G: 0`, which hid a transport-limited double count for months. When a code path branches on a parameter, a guard needs a fixture on each side.
- A regression guard must fail without its fix. Verify that before committing.
- A code path with no test rots silently: plate advection with stratigraphy crashed on an attribute removed in a 2024 refactor (`idNbghs`) for two years, and hid a second fault behind it. Give every configurable path at least one smoke run.
- A function whose only caller is a serial test hides np>1 bugs (`stratalRecord` set full local columns into a global scratch Vec: right only at np=1). Any helper that touches a global Vec needs an np>1 guard.
- A cyclic model's mesh is 3-D (a cylinder) but its inputs are in the flat frame: any geometric operation (displacement, departure point) must map through `_cylinderVelocity`/`_cylinderDisplace`, never add the flat vector to the 3-D coordinates.
- Memory that scales with the MESH on every rank is invisible at low rank counts (local data dominates) and decisive at high ones; measure the non-root ranks at two rank counts.
- `flowAccumulation` runs twice per step, so it must not overwrite the per-step forcing it reads: `bL`/`bG` are read-only there and the net source (evaporation, ice split) goes to `_srcL`/`_srcG`. Writing it back evaporated the runoff twice (2026-10-27). `UnstMesh._setRunoffSource` is the only writer of `bL` (rain and orographic paths alike).
- A face coefficient set by a state on one side (deposit thickness) must come from the cell material leaves (the donor), not the average of the two: averaging lets a bare cell shed bed into a deposit, and the `>= 0` clip after the solve keeps only the gain. A global mass rescale hid this for months by halving every marine deposit (2026-11-12, `_donorDiffCoeffs`). A rescale factor far from 1 is a symptom, not a fix. A face-by-face fix is not enough when an implicit step can move more than a cell holds: the same clip created volume again at a large diffusion number until a supply limiter capped each cell's export (it only showed at 10 km resolution, not on the 28 km test mesh).
- A threshold derived from one surface must be computed on the surface the algorithm actually uses: the marine eps-fill cut-off came from the raw elevation while the fill ran on the smoothed one, so it never had a seed on a global mesh (no fill, ~660 terminal sinks, km-thick force-deposits).
- PETSc fieldsplit sub-solver settings made with `subksp.setType()` before `setFromOptions()` do not survive the switch to Schur: set them through the options database and check with `-<prefix>ksp_view`. Both coupled deposition solves silently ran GMRES in their blocks (SPL solve ~5,000 s/step at 5 km).
- A non-fatal solve that is discarded on failure must not be one whose result carries mass: the sediment routing solves were zeroed every step on fine global meshes until they got the exact-factor rescue (`route=True`).
- Two rules that route the same cells must read the same surface: flat routing on `lFill` and MFD on `h` differed by 1-ulp round-off and formed closed flow cycles. Make them consistent at the routing step, not by clamping the fill (that fragmented DEM plateaus into ~1e5 sinks). Find the cycle's cells (the routing-repair diagnostic) before adding more solver rescue, and check a fix on the production case, not only on the guard test.
- A composition field read by a shared routine (`depoFineFrac` in `deposeStrat`) must be set by EVERY caller for its own deposits: hillslope creep inherited the stale routed value and created fine sediment in proportion to the diffusivity (+19 % at 1e7), invisible at the fixture's tiny diffusivity. Test conservation at the parameter extremes, not only at the fixture values.
- A failure localised to a few cells must be repaired locally, not by discarding the whole solve: a 3-cell flow cycle zeroed a full sediment-cascade pass (`_pondRoutingRegion` now makes those cells sinks and re-solves).
- An adaptive solver that can give up (TS `DIVERGED_STEP_REJECTED`) needs a fallback, not an abort: a 144-rank job should not die on an option that has a robust alternative.
- A non-linear flux law is implemented as its flux factor `q/(D S)`, not its derivative `dq/dS`, and the factor is evaluated on the gradient magnitude, not the edge slope (`S cos(theta)` on an unstructured mesh). Both errors survived two years because no option had a benchmark against its exact solution (`hillslp_nl`, 2026-11).
- Operator splitting decides what an output shows: uplift applied after the processes left one un-eroded `U*dt` in every output and made the steady relief depend on `dt`. The uplift runs at the START of a step (`applyVerticalForcing`, 2026-11); keep it there. The open-edge reset stays at the END (`applyForces`): run before the processes, its flat edge runs made outlet capture partition-dependent (`test_advection_partition_invariant`). Check partition invariance at several rank counts (`scripts/ab_partition.py -n 2/3/4`): one rank count can pass by luck. The hillslope step has two positions: after the deposition (default), or between the river incision/glacial till and the sediment routing when `diffusion: channelQ` is set with deposition on (`hillslopeBeforeRouting`), so its channel supply is routed in the same step; a kernel that writes the river step's routing arrays (`thCoarse`/`thFine`/`provEro`) between SPL and `sedChange` must save and restore them.
- When a tolerance fixes a knife edge, apply it at EVERY test of the same quantity: the flow fix gave `FILL_LEVEL_TOL` to `flowplex` but the groundwater model kept `lFill > z` in three places, and a lake cell at its own surface flipped recharge with the rank count (head 44 m off at np=3, 2026-11-24). `grep` for the comparison, not the function.
- A solver that reports convergence has met one of its criteria, not necessarily yours: the soil and non-linear SPL SNES stopped on PETSc's default step-size test (`stol = 1e-8`) with residuals 1e3 times `rtol`, a partition-dependent accuracy (2026-11-24). Check `getConvergedReason` and switch off the defaults you do not mean.
- A cached solver keeps every setting of its previous calls: `_solve_KSP(guess=False)` only set the initial-guess flag in its `guess=True` branch, so zero-guess solves warm-started from stale Vecs (2026-11-24). Set per-call options on every branch.
- A threshold on a field that takes discrete values is a tie generator: on a regular test mesh with uniform rain the discharge only takes multiples of cell area × rain, so `channelQ` set exactly on one flips cells by round-off and looks like a ~1 m partition bug. Put test thresholds off the lattice before suspecting the code.
- Inverse-distance weighting of the nearest displaced nodes is biased towards the nearest one when the displacement is a fraction of a cell (a 0.2-cell remap moved a bump 70 m of 80 m): use weights that reproduce a linear field (barycentric) wherever displacements can be sub-cell.
- A node can have a physical role (a deep fill seed linked to the open sea) and a parallel role (partition boundary) at once; a flag that encodes only one of them makes the physics depend on where the boundaries fall. The depression fill lost the sea link of deep seeds on partition boundaries, which only the ubuntu partition hit (seafloor filled from -5590 m to +317 m, 2026-11-24). A partition-exact result on one platform says little about another: ParMETIS splits differently, so check `pytest -m mpi` on the Linux CI too.
- A restart must reproduce the continuous run, and only a comparison with one shows it: counters that are "the next value" (layer index, next layer time, the forcing of the next step) were one step ahead after a restart for years, and the outputs a restart reads were single precision (the 1-ulp depression-fill tilts were lost). Any state carried between steps must be saved and restored; `tests/test_restart.py` compares a restarted run with a continuous one (2026-11-24).

## Intentional surprises (do NOT "fix")
- **`gid` argument on `mfdreceivers` / `mfdrcvrs`**. The Fortran kernels take an extra `gid(nb)` argument (per-local-node vertex ID `self.gid`, set in `mesher/unstructuredmesh.py`; since 2026-06-21 it is `self.locIDs`, the partition-INVARIANT input-mesh id, no longer PETSc's per-partition numbering; see `gospl/flow/AGENTS.md`). Inside the kernels, the stored slope is perturbed by `val * (1 - 1.0e-15 * gid(n))` before the quicksort tie-break. This makes EXACT slope ties resolve deterministically across MPI decompositions (lower global ID wins). The perturbation is well below KSP-solver precision and does not affect physical results; near-ties driven by KSP floating-point noise are a separate problem and are NOT fixed by this pattern. See the relaxation comment on `test_parallel_correctness` in `tests/test_parallel.py` and `fortran/functions.F90:mfdreceivers` for the full rationale. **Do NOT remove the `gid` argument** — without it, `mfdreceivers` ordering depends on local iteration order and the parallel-correctness test diverges by ~10x from current baseline.
- **Loose tolerances in `test_parallel_correctness`** (`rel_sum_fa` 5e-2, `rel_mean_h` 5e-9): the platform-dependent KSP floating-point floor, not a regression. Tightening them requires bitwise-identical halo state across decompositions. Detail: `docs/dev/FIXED_BUGS.md` > Test tolerance notes.
- **`filllabel`** is uncalled but kept on purpose; do not remove it in a dead-code sweep (`fortran/AGENTS.md`).
- **Fortran `meshparams` allocatables** must be freed per mesh in `definetin`; see `fortran/AGENTS.md`.

## CI contract (summary)
- `tests-pr.yml`: `pytest tests/` + the MPI lint on PRs and pushes to `master`/`release-candidate`, and on manual dispatch (`gh workflow run tests-pr.yml --ref <branch>`), the only way to run the FULL suite on `dev` or a feature branch in CI (a pull request into `dev` starts NO tests): the nightly `tests-slow` runs just the `slow` subset (69 of 278 tests) and never the benchmarks. **Pushes to `dev` run no CI**; the nightly cron (`tests-slow.yml`, `examples-smoke.yml`) runs from the default branch, which is `dev`. So `dev` is tested nightly, not per push: run the suite locally before pushing.
- `tests-slow.yml`: slow tier nightly/tags/dispatch, benchmarks on `master`/`release-candidate` pushes. `examples-smoke.yml`: every local goSPL-examples input for 2 steps at np=2 (`scripts/examples_smoke.py`), nightly + dispatch. The Global-examples inputs are notebook products, gitignored in goSPL-examples, so CI cannot run them (a dispatch skips them all, and the script exits 2 when nothing ran): run `python scripts/examples_smoke.py <goSPL-examples> --only Global-examples` locally.
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
| `test_hillslope_nonlinear.py` | Non-linear hillslope (Taylor series, Andrews-Bucknam/Roering, power law) | steady ridge `D S F(S) = U x'` integrated exactly, hexagonal lattice | each law no more than 1 % (divide, RMSE) less accurate than the linear law on the same mesh (2026-11) |
| `test_marine_diffusion.py` | Marine diffusion (`ts` and `picard`) | 2-D heat kernel on a constant-`Cd` deposit | volume 1e-12; D err, peak law, RMSE < 2e-3 for both solvers (the `ts` stage-solve fix, 2026-10) |
| `test_geochem_balance.py` | Level-B geochemistry | per-step `dissolved = precipitated + exported`; 1-D strip `c=D0/(R+p)·(Ls/u)^(1+p/R)` | closure ~1e-15; strip vs exact discrete 2e-5 |
| `test_gravity.py` | Gravity-current run-out (`_gravRoute`) | uniform slope, square lattice: `D_k = Q0 (1 - e^(-dx/L)) e^(-k dx/L)` | round-off at dx 1000/500 m, 1 and 6 receivers; fitted L = L |

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
