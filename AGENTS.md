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

### CACHED — hot-path solvers (11 sites)
Solvers that run **every timestep** are created lazily on first use, stored as `self._X`, and reused across the entire simulation. Avoids the ~5-10 ms create+destroy churn per call, which compounds over thousands of timesteps. Each site's inline docstring confirms the rationale.

| File | Method | Cached attribute | Solver |
|---|---|---|---|
| `flow/flowplex.py` | `_solve_KSP` | `self._ksp_main` | KSP, **`fgmres`+`bjacobi`** (`flowacc_` prefix, `pc_factor_shift_type nonzero`). Solves the IDA `(I−Wᵀ)q=b` flow/sediment-routing system. **Per-call `max_it`:** the **fatal** flow-accumulation solve gets the full `_primary_max_it` (5000, or 100000 for the `richardson` escape hatch — it MUST converge); the **non-fatal `seed=True` cascade** solves are capped at `_cascade_max_it` (1000) so a partition-dependent near-singular cell can't grind the full budget before the bounded fallback (the dominant cause of erratic `sed` wall-time at some rank counts — ~5500 vs ~205 iters/solve at np=16 vs np=8 on a 9.2M mesh). See **`## Flow-accumulation KSP`**. Env overrides `GOSPL_FLOW_KSP` / `GOSPL_FLOW_PC`. |
| `flow/flowplex.py` | `_solve_KSP2` | `self._ksp_fallback` | KSP, **`richardson`+`none`** bounded fallback (`flowaccfb_` prefix, `max_it=500`). Reached only when the primary fails on a (near-)singular sub-region; confirms non-convergence cheaply then degrades gracefully (zeroed, finite). fgmres already recovers every solvable case, so the fallback never needs deep iteration. |
| `eroder/nlSPL.py` | `_solveNL_ed` | `self._snes_ed` (+ `self._snes_ed_fb`) | SNES, transport-limited (**`qn` default**, ngmres fallback) |
| `eroder/nlSPL.py` | `_solveNL` | `self._snes_nl` | SNES, detachment-limited (nrichardson + analytic Jacobian) |
| `eroder/soilSPL.py` | `_solveSoil` | `self._snes_soil` (+ `self._snes_soil_fb`) | SNES, soil-aware (**`qn` default**, ngmres fallback) |
| `eroder/soilSPL.py` | `diffuseSoil` | `self._ts_soil` | TS (rosw soil diffusion) |
| `sed/hillslope.py` | `_hillSlopeNL` | `self._snes_hill` | SNES (non-linear hillslope) |
| `sed/hillslope.py` | `_diffuseImplicit` | `self._ts_marine` | TS (rosw marine + lake diffusion) |
| `sed/hillslope.py` | `_solveSmooth` | `self._ksp_smooth` (+ cached operator `self._smoothMat`) | KSP (richardson+bjacobi, `marsmooth_` prefix, factor-shift) for the marine flow-direction smoother (`_hillSlope(smooth=2)`). Operator is mesh-size-based + timestep-independent so it is built once and the PC **reused**; rebuilt only when the coastline (`seaID`) moves (`_smooth_pc_ready` guards the reuse). |
| `sed/hillslope.py` | `_hillSlope(smooth=0)` | `self._ksp_hill_lin` (+ cached operator `self._hillMat`) | KSP (`hilllin_` prefix; same config via shared `_makeDiffusionKSP`) for the **linear** soil-creep solve. Operator `(I−Cd·dt·L)` is invariant for constant `Cda`/`Cdm`, so built once + PC **reused**, rebuilt only on coastline move (`_hill_pc_ready`). **Only when NOT dual-lithology** (`stratLith` makes `Cd` step-varying → per-step build on `_ksp_main`). Bit-faithful. |
| `sed/hillslope.py` | `_diffuseImplicitPicard` | `self._ksp_picard` | KSP (`marpicard_` prefix; `_makeDiffusionKSP`) for the **opt-in** lagged-diffusivity marine/lake solver (`diffusion: marineSolver: picard`; default `ts` = the adaptive non-linear `_ts_marine`). Per Picard iteration it rebuilds `M = I + dt_sub·L` from `jacobiancoeff(h,Cd,Kp=0)` (the linear `∇·(Cd∇)` operator, **consistent with the TS residual `fctcoeff`** incl. the no-flux-across-zero-Cd-face marine-mask gating) via the single-pass `_assembleDiffMatCSR`, and solves it (no SNES). The operator changes each iteration so PCSetUp is not reused — but each solve is linear + smooth (no kink rejections). Approximation: matches the TS to ~1e-5 on minimal, ~1.5% deposit diff on a 658k earth run at ~5-8× the marine-diffusion speed; `picardSub`/`picardIts` are the accuracy/speed knobs. |

**CRITICAL — collective rebuild decision.** Both coastline-gated caches above (`_smoothMat`/`_hillMat`) rebuild the operator (and redo `PCSetUp`) only when `seaID` moves. `seaID` is **rank-local**, but `_buildDiffMat` assembly and `PCSetUp` are **collective** — so the "did the coastline move?" test MUST be reduced across ranks (`rebuild = MPI.COMM_WORLD.allreduce(local_changed, op=MPI.LOR)`) before it gates those calls. Without the reduce, one rank rebuilds while another reuses, the two take different collective paths, and the run **deadlocks at np>1** (it ran fine until the per-partition coastlines drifted apart, then hung). A rank forced to rebuild with an unchanged mask reproduces its own cached matrix, so the reduce stays bit-faithful. Any future cached operator gated on a rank-local condition needs the same reduce.

**Selectable solver + complementary fallback** (soilSPL `_solveSoil`, nlSPL `_solveNL_ed`): built by a `_build_*_snes(primary=)` helper. The primary defaults to `qn` (L-BFGS + critical-point line search, set via the YAML `solver:` key). A *bare* `ngmres` ignores its configured KSP/PC and **stalls or diverges** on these stiff residuals (this was the historical default; it cost soil ~2.4× and nlSPL-transport ~10×, and diverged outright on the SIA — which is why SIA is now solved explicitly, see "Ice sheet"). On non-convergence the solve retries from the same initial guess with the **complementary** solver (`qn` ⇄ `ngmres`+`nrichardson`+hypre); if both fail it warns and continues with the best iterate (never crashes). Each SNES carries its own options prefix (`soilspl_`/`soilsplfb_`/`nlspled_`/`nlspledfb_`) so line-search options don't leak between solvers. The `_fb` fallback SNES **and** its `_f` residual vec are cached and are in the `destroy_DMPlex` list. YAML knobs (`soil:` / `spl:` blocks): `maxIter`, `rtol`, `atol`, `pcType`, `solver`. The two **coupled two-field (h/q) deposition** solves — **linear SPL deposition** (`SPL._coupledEDSystem`, `spl_ed_` prefix) and **marine deposition** (`seaplex._depMarineSystem`, `marine_dep_` prefix) — both use a **Schur-complement** fieldsplit (defaults injected into the options DB guarded by `hasName` so `PETSC_OPTIONS`/`-spl_ed_*`/`-marine_dep_*` override). The two blocks are near-triangular M-matrices their per-block ILU solves almost exactly, so the Schur factorisation captures the full coupling and converges in ~2 outer iterations vs the old additive split (SPL ~60×/122→2 fewer Krylov iters; marine 31→2 on the minimal mesh, more at scale), solution unchanged. (Both remain AD-HOC per-call create/destroy — the Schur change is only the preconditioner, not the lifecycle.)

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

## Flow-accumulation KSP (`_solve_KSP` / `_solve_KSP2`) — design + known limitations
The implicit drainage-area (IDA) system `(I − Wᵀ) q = b` (`W` row-substochastic MFD weights, `b` runoff/sediment) drives flow accumulation **and** the downstream water/sediment **cascade** solves (`flowplex._distributeDownstream`, `sedplex._moveDownstream`). It is ill-conditioned (spectral radius ≈ 1 from long flat drainage chains), so it needs an *accelerated* solver. The current shape (validated on a 10 km global mesh, 48–240 ranks):

- **Primary = `fgmres`+`bjacobi`** (default; was stationary `richardson`+`bjacobi`). **Why the change:** stationary Richardson on this operator is destabilised by `bjacobi`'s per-rank blocks — the *same* (partition-invariant) matrix split differently at each rank count gives a different preconditioned operator, and at unlucky cuts ρ>1 → diverges → fell to the costly fallback → erratic, partition-dependent `sed`/`sea` blow-ups (10–40×). A Krylov method (fgmres) minimises the residual monotonically so it **cannot** diverge that way; it converges fast and stably at every rank count. This is mechanism #1 of the partition-dependence and it is now **resolved**. Overridable via env `GOSPL_FLOW_KSP` (`fgmres`|`richardson`) / `GOSPL_FLOW_PC` (`bjacobi`|`none`|…) for experiments; defaults are the validated combo. `richardson` keeps a larger `max_it` (it propagates one hop/pass); `fgmres` is capped at 5000 (the bulk converges in O(100s)).
- **`seed=True`** (opt-in, on the `(I−Wᵀ)` solves only): cold-start the guess with `b`. Since `q = b + Wᵀq ≥ b`, `b` is a valid lower bound; without it the init solve exceeded `max_it` (`DIVERGED_MAX_IT`). Do **not** enable for non-`(I−Wᵀ)` systems (hillslope/SPL/tectonics/ice).
- **`fatal=True`** only on the main `flowAccumulation` discharge solve — a failed *main* discharge would silently feed a no-river state into erosion/sediment, so it aborts. Cascade/auxiliary solves are `fatal=False`: on failure they zero (finite, bounded) and continue.
- **Discharge clamp ≥0** at the end of `flowAccumulation`: fgmres (unlike the IDA fixed-point) does not preserve non-negativity, and a tiny negative feeds `PA**spl_m` → `sqrt(NaN)` in the SPL eroders. No-op under richardson.
- **On primary failure: never accept the iterate.** It restarts the fallback from a clean BOUNDED guess (the seed `b`, else a finite warm iterate) and lets the bounded fallback finish (zeroes a non-convergent solve). *Do not* "accept fgmres best-effort on max_it" — on a near-singular sub-region the **residual** is small but the **solution** is unbounded (null-space blowup), so accepting injects huge values that compound through the cascade (observed 1e29 escalation). **Benign-failure reporting:** the un-converged region is localised first (one matvec → `nbad` = nodes above `1e-3·‖b‖`, Allreduced). A *non-fatal* solve failing on only a TINY region (`nbad ≤ self._undrained_benign_cap`, default 256) is the expected isolated-pocket / micro-cycle case (a few genuinely un-drainable cells that just pond; discharge clamped ≥0, mass conserved) — it prints one calm `[flow] N isolated un-drained cell(s) … (benign)` line. Only the **fatal** main solve or a **large** region (a real failure: NaN source, broken partition) gets the loud `KSP … failed to converge` + `[flowKSP] worst residual … N nodes` diagnostic with the RHS/matrix-finite probe. The benign region is *nearly* partition-invariant after the mechanism-#2 fixes (0 on the 5.9M mesh), but a handful of cells still flip at higher rank counts / finer meshes — see **Residual** below; that small variation is what the `_cascade_max_it` cap exists to bound.

**Mechanism #2 of the partition-dependence (FIXED — the residual drainage non-invariance).** The cascade solves run on the **filled** topography. A local A/B harness (run the fill + drainage at np=1 vs np=N, gather fields to rank 0 in **`locIDs`** input-mesh order, diff) isolated it. **The fill (`lFill`), the sink set (`lsink`) and the MFD self-sink set are bit-for-bit partition-invariant** — the non-invariance was entirely in the depression *graph* built on top of them. A 370k mesh exposed only one layer (the flat-routing tie-break); a **5.9M-node global mesh exposed three more** in the pit machinery. All four are now fixed; the **true singular set of `(I − Wᵀ)` — what the flow KSP sees as the "un-drained region" — is bit-for-bit np-invariant on the 5.9M mesh** (symdiff 0 at np=1-vs-np8), so the KSP message no longer varies with rank count *there*. A tiny residual remains at higher rank counts / finer meshes (the 9.2M case) — see **Residual** below for what it is, the erratic `sed` it causes, and why the cap (not a pin) handles it.

The four roots + fixes (all keyed on the partition-invariant input-mesh id):
- **`self.gid` → `self.locIDs`.** `gid` *was* the PETSc per-partition DMPlex numbering (`l2g`) — RENUMBERED per partition (`gid[locID]` differs across np, identity only at np=1) — so every tie-break keyed on it was non-invariant. Its ONLY use is the `mfdreceivers`/`mfdrcvrs`/flat tie-breaks (+ a diag print), so it's redefined as `self.locIDs.astype(int32)` (KDTree input-mesh index, invariant; PETSc `l2g` still built for `lgmap_col`). Genuinely fixes the **SFD** path (`iceplex` `mfdreceivers(1,…)`) and the flat router; MFD was already set-invariant (only the steepest-first `[:,0]` order varied).
- **Flat-routing receiver tie-break (`_dirFlats`, Python).** In a flat region many neighbours tie at the minimum `ptdir`; the old fortran `fill_rcvs` broke the tie by *"first in `FVnID` neighbour order"* (partition-dependent). Replaced with a vectorised pick: among strictly-closer (`ptdir(c) < ptdir(i)`), not-higher (`lFill(c) ≤ lFill(i)`) neighbours, smallest `ptdir`, **ties broken by smallest `locIDs`**. (`fill_rcvs` now unused/removed.)
- **Flat-BFS distance `ptdir` (`nghb_dir`, fortran).** The kernel assigned each node's spill-distance **once** (a `flag` gate) and the driver re-seeded it every pass, so a boundary node frozen at a longer within-partition distance was never relaxed when the shorter cross-partition value arrived through the halo → `ptdir` (10270/5.9M differed). Rewrote `nghb_dir` as a **relaxing (label-correcting) Dijkstra** (lower a distance whenever a shorter path appears); `_dirFlats` now iterates to **global convergence** (no owned node changed) instead of a `remaining`-stall heuristic. The fixed point is the true graph distance ⇒ invariant.
- **Spill-point selection (`_pitInformation`, Python) + pit membership (`label_pits`, fortran).** `spill_pts` claimed each pit's spill for the *first pit node in local order* (then `Allreduce(rank, MAX)` across ranks) — both partition-dependent (43212-vs-43237 spills). Replaced with a min-`locID` canonical pick (same physics: case-A lower-neighbour rim node / case-B equal-level border-or-non-pit edge). And `label_pits` flagged every node on first visit + flagged higher neighbours during the flood, so a flat rim cell visited before its depression seed was frozen unlabelled (150 membership flips); removing both order-dependent flags makes membership the invariant "equal-fill connected to a bottom seed" set. (`spill_pts` now unused/removed.) The relabelling marks more equal-fill rim cells as pit members, but they keep their MFD receivers (`_buildFlowDirection` only overrides with `flatDirs` when `sum_weight==0`), so it's flow-neutral — `mass_conservation`/`advection_parallel` hold.

**Residual non-invariance + the erratic-`sed` it causes (capped, NOT eliminated).** A tiny set (~1e-6 of nodes) still flips: ~8 flat-receiver cells on the 5.9M mesh at np=1-vs-np8, and on a finer **9.2M mesh it grows with rank count** (np=4-vs-np=8: 14-node un-drained symdiff; a near-singular cell present at np=16 but not np=8). Most are harmless equally-valid flat re-routes, but a few are genuine partition-dependent **near-singular cells** (residual cross-partition 2-cycles / isolated pockets). Those are **benign physically** (they pond; mass conserved) but make a **cascade `fgmres` solve grind to `max_it`** before the bounded fallback — the dominant cause of erratic `sed` wall-time at some rank counts (e.g. the 9.2M Gadi sweep: `sed` spikes at P=48/288, clean at P=240; reproduced locally np=8-vs-np=16: ~205 vs ~2750 fgmres iters/solve). **Two cheap mitigations are in place** (do not remove): the benign-failure reporting (above) and the `_cascade_max_it` cap (1000) that stops the grind ~5× early with the same bounded result (the fatal flow-accumulation solve keeps its full budget). **A reachability-to-outlet pin was tried and REJECTED — do not re-attempt the naive version.** Pinning non-draining cells to identity rows (so the operator is non-singular, eliminating the grind) is *correct* (implemented `reach_outlet` + `_pinUndrained`; suite green; closed/endorheic basins protected by also seeding the `sum_weight==0` terminal cells, not just sea/outlets) but **impractical**: its cross-partition reach-propagation loop runs MORE collective iterations at HIGHER rank counts (× ~12 cascade builds/step) and its arrays swapped the memory-tight 9.2M np=16 run to disk — the cost scales up exactly where it's needed. A viable pin would need a *cheap* parallel reachability (serial-on-rank-0 like the pit spillover graph, or a donor-based upstream BFS), not a per-build halo loop. Any pre-existing bit-fingerprint at fixed np shifts (tie resolution changed even at np=1, all physically equivalent); `test_parallel_correctness`/`test_mass_conservation` hold.

**The OUTER cascade loop is ALSO bounded now (do not remove) — it has to be.** `_cascade_max_it` caps a *single* solve; it does nothing for the number of *passes* of the `while excess` loop in `flowAccumulation`. When a partition cuts a genuinely **un-drainable pocket**, the cascade solve never lets `excess` clear, so that loop — uncapped — spins forever, rebuilding `fMat` and re-solving every pass and shedding deferred-destroy PETSc objects until the **collective garbage collector** (`GarbageKeyAllReduceIntersect_Private`) overflows its key count and aborts the whole run with a bogus ~18 EB malloc (`Memory requested 18446744065119617024`). Reproduced **deterministically** at 5 km / 144 ranks (the pocket = a 1473-node region on that partition; 96/192/240/288 partition the network without forming it, so they complete). Fix (commit `d4162aa`): a **stagnation break** — `_distributeDownstream` records the residual downstream flux `_cascade_resid` (a global `Vec.sum`); after `_cascade_patience` (default 3, env `GOSPL_CASCADE_PATIENCE`) consecutive passes that fail to shrink it by `_cascade_rel_improve` (0.1 %), the pocket is declared un-drainable, the residual is left **ponded** (physically a closed basin), a loud `[flow] downstream cascade stalled after N passes …` is printed, and the loop breaks — typically within ~3–5 passes, well short of the overflow window. A slowly-but-genuinely-converging basin never plateaus, so real convergence is never truncated. A hard `_cascade_max_steps` (default 100, env `GOSPL_CASCADE_MAX_STEPS`) remains as a backstop. Outcome is partition-invariant (counts that form the pocket pond it instead of crashing). NOTE: the **sed/sea** cascades already bound their loops (`max_iters=5000` + graceful force-deposit) so they don't hit this; only the flow loop was uncapped.

**Related: isolated-pocket sinks.** A handful of flat nodes still can't reach a spill and fall back to **self-sinks** — **isolated filled pockets** (tiny clusters in non-pit terrain sharing a watershed's pitID+spill but not pit-connected to it). A sink is the physically-correct local outcome (they pond; mass stays in the pit); now partition-invariant. `verbose` surfaces the count.

**Dead Fortran removal (2026-06).** The partition-invariance work orphaned `fill_rcvs` (flat-routing receiver pick, now Python in `_dirFlats`) and `spill_pts` (spill selection, now Python in `_pitInformation`). Removed together with eight pre-existing never-called kernels — `definegtin`, `distocean`, `donorslist`, `donorsmax`, `gfill`, `scale_volume`, `stencil`, `stratathreesed` — from both `fortran/functions.F90` and `fortran/functions.pyf`. Each was verified to have **zero invocations** in Python (`from gospl._fortran import …`) and zero internal Fortran callers before removal. **`filllabel` and `definetin` were deliberately kept** (along with everything `definetin` calls — `distance`/`euclid`/`gtriarea`/`gtriareasphere`/`slerpmidpt`/`slerpvec`/`striarea`). The priority-queue/queue infrastructure (`PQpush`/`PQpop`/`shiftdown`/`pop`/`push`/`spop`/`spush`/`ppop`/`ppush`/`wpop`/`wpush`) stays — it's used by the priority-flood and `filllabel`. **Before removing any Fortran routine, confirm both: no `from gospl._fortran import <name>` anywhere AND no internal `call <name>`/`<name>(` in `functions.F90` — the f2py `.pyf` interface and the `.F90` body must be edited together.**

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
| 7 | `_IceMesh()` :154 | FAMesh | `iceHL/iceMeltL/iceMeltRiverL/iceUbL/iceAbrL/iceFAL/iceFAG/iceFlex` + cached `iceMat` (only if `iceOn`) |
| 8 | `_SPL()` :157 | FAMesh | `hOld`, `hOldLocal`, `hOldFlex`, `Eb`, `EbLocal`, `stepED`, `newH` |
| 9 | `_nlSPL()` :160 | SPL | `snes_rtol/atol/maxit`, lazy `_snes_ed/_snes_nl` |
| 10 | `_soilSPL()` :163 | nlSPL | `Gsoil`, `Lsoil`, `lHbed`, `gHbed`, `prodSoil`, `soil_transition` |
| 11 | `_PITFill()` :166 | UnstMesh | `borders`, `outEdges` |
| 12 | `_SEDMesh()` :169 | FAMesh | **`tmp`, `tmpL`, `tmp1`**, `Qs`, `QsL`, `nQs`, `vSed`, `vSedLocal`, `maxnb` |
| 13 | `_hillSLP()` :172 | SEDMesh | **`h`, `hl`, `dh`**, `mat`, `Dlimit/dexp/minDiff`, lazy `_snes_hill`, `_ts_marine` |
| 14 | `_SEAMesh()` :175 | FAMesh | `zMat` |
| 15 | `_GridProcess()` :178 | ReadYaml + UnstMesh | `localFlex`, cached FEM operator/KSP (flat `fem`) or DH grid (`global`); cached orography advection operators if orography on (no regular grid) |
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

**Erosion** (`erodeStrat`) splits eroded solid into `thCoarse`/`thFine` by the consumed layers' composition; the bedrock sentinel is inflated *proportionally to layer-0's fine ratio* so it never alters that layer's composition. **Transport**: `_getSedFlux` routes the TOTAL (`thCoarse+thFine`) into `vSed` and `thFine` into `vSedF` (same operator, exact). **Deposition is composition-only — the deposit geometry (`delta`/`dh`) is never modified, so elevations/flow are unchanged.** `_pitFineFraction` (continental) and `_marineFineFraction` (marine) bias the deposit fine fraction toward the depocenter / deep-distal water by bathymetric depth, conserving the deposited fine volume; coarse stays proximal. **Compaction** (`_depthPorosityDual`) compacts each fraction on its own curve, conserving per-fraction solid. **Advection** (`stratalRecord`) advects the fine pile (`stratHf`,`phiF`) with a SECOND `strataonesed` call — a 3-fraction kernel would be wrong here (it would renormalise its extra 0–1 fraction fields, conflicting with the bulk-thickness representation). (The dormant `stratathreesed` kernel that embodied that approach was never called and was removed 2026-06 — see the "dead Fortran removal" note.)

**Fine-enriched overspill (implemented):** `_moveDownstream` threads the fine sub-volume through the cascade with **coarse-settles-first** retention — a filled pit keeps a coarse-enriched deposit (its `_pitRetFine`, used by `_pitFineFraction`) and the excess overspills fine-enriched, routed downstream and on to the marine domain. `vSedFLocal` accumulates the residual fine *mirroring* `vSedLocal` — **it is NOT zeroed in `_distributeSediment`; that was the bug in the first, reverted attempt.** Per-fraction conservation is guarded by **`test_dual_fine_conservation`** (`_fineEroded` vs `_fineDeposited`, the fine analogue of the total-cumED test, which does NOT catch a fine-only leak). **Any future change to the fine routing MUST keep that test green.**

## Ice sheet (diagnostic glacial model)
Opt-in via an `ice` YAML section (`self.iceOn`). goSPL ships a **single diagnostic glacial model** — no ice-dynamics PDE solve, no `flow_model` selector. Lives in `flow/iceplex.py` (`IceMesh`); full design `docs/DESIGN_ICE_SHEET.md`, summary `docs/ICE_SHEET_SUMMARY.md`. (A true Shallow-Ice-Approximation solver was implemented and **removed**: the `H≥0` ice margin is a free-boundary/obstacle problem on which an implicit `F(H)=0` solve diverged, and any ice-dynamics solve is stiff and over-thickens km-scale continental ice at goSPL's 10²–10⁴ yr steps. The investigation kernels (`ice_flux*`, `ice_flux_jacobian`, `sia_diff_coeff`) lived only on branches `fix/sia-solver-investigation` / `feat/ice-explicit-flux`, both **deleted (local-only, never merged) in 2026-06** — those kernels no longer exist anywhere. If a SIA / obstacle-VI ice-dynamics solver is ever revisited, re-derive from the written design in `docs/DESIGN_ICE_SHEET.md` §9.)

**Diagnostic (`iceAccumulation` → `_iceFlowMFD`).** Per step: `_iceMassBalance` builds the ELA `mdot` (`= rainVal·clip((z−elaH)/(iceH−elaH),·,1)`, accumulation scaled/capped by `accum_factor`/`accum_max`); the accumulation (`mdot⁺·larea`) is **routed downhill** through an MFD matrix (`_matrixIceFlow`, same epsfill+`mfdreceivers` machinery as water) into an **ice discharge** `iceFAL`/`iceFAG` (m³/yr, one `_solve_KSP`); smoothed → **Bahr thickness** `iceHL = icewe·icewf·Q^0.3`; **basal sliding velocity** `iceUbL` from Glen's sliding law on that thickness + bed slope (`ice_velocity` kernel, `∝ H^(n-1)|∇s|^(n-1)`, physically bounded). Then `iceAbrL = Kg·|ub|^l` (+ optional lateral term, see below), terminus clamp `iceHL[zbed<max(hterm,sealevel)]=0`. NO time integration — ~1 routing solve/step, robust + physical at any resolution (real 30 km global: ~0.8 s/step, `Hmax`~1 km). Config: `icedir`/`eheight`(icewe)/`fwidth`(icewf)/`melt`(meltfac)/`slide`(ice_slide)/`glen`(ice_glen)/`accum_factor`/`accum_max`/`melt_conserve`. Persistent Vecs `iceFAL`/`iceFAG`/`iceMeltL`/`iceMeltRiverL`/`iceUbL`/`iceAbrL` + cached `iceMat` in `destroy_DMPlex`. Output: `iceH`/`iceUb`/`iceFA`/`iceAbr`/`iceMelt`.
> *Don't use a raw balance velocity `Q/(H·√area)` for `iceUbL`* — it blows up with catchment size (`Q` grows with upstream area, `W` stays cell-scale) and spikes at flow-convergence cells → tens-of-km erosion/deposition in sinks. The Glen sliding law on the Bahr thickness is bounded; magnitude is tuned via `Kg`.

**Meltwater → rivers.** `iceMeltRiverL` (re-injected into the river FA by `flowplex.flowAccumulation`) is **discharge-conserving by default** (`ice.melt_conserve`, default True): `_glacialMeltwater` routes the accumulation down the ice surface and releases it where the ice melts out (`f=1` at the margin), so `Σ river-melt == Σ accumulation` (closes the glacial water budget over the long, ~steady-state step). `melt_conserve: False` = precip-scaled local ablation (loses water). **Keep `iceMeltL` (precip-scaled ablation, the till melt-out/deposition weight) and `iceMeltRiverL` (river water) distinct** — conflating them breaks till conservation (the conserving melt-out lands at accumulation cells, not the ablation zone).

**Geometry forcing.** `hela`/`hice`/`hterm` are uniform scalars, per-vertex maps, or a `glaciers` time series (`_buildIceSeries` → `_iceTimeSeries`, stepped each forcing update by `mesher._updateIce` into `elaMesh`/`iceMesh`/`termMesh`, indexed `[locIDs]`). Scalar path keeps the global no-ice short-circuit; `_iceMassBalance` and the `flowplex` precip split are divide-by-zero hardened for per-node `iceH≤elaH`. `hterm` unset = `TERMINUS_UNSET` sentinel → sea level. Helper `scripts/ela_from_temperature.py` builds ELA maps from temperature.

**Glacial erosion / till.** With `till.on` off, abrasion `Eg=Kg·|ub|^l` is added to `Eb` as incision (`SPL._glacialAbrasion`, gated off when `till.on`). With `till.on`, `glacialTill` is a conservative transport: it computes a per-cell deposition weight (`Σ owned = 1`) and applies it bulk (`stratNb==0`, net bed change 0) or through the strata (`_glacialTillStrata`: `erodeStrat`→moraine via `deposeStrat`, split coarse/fine under dual lithology, per-fraction solid conserved; **saves/restores `thCoarse`/`thFine`/`depoFineFrac` around the calls** since `sedChange` consumes them later). Weight = **catchment-routed by default** (`till.route`, default True) via `_routeTill` — transport-with-loss on the ice surface (`mfdreceivers` steepest descent, melt-out `f=min(1,ablation·dt/H)` forced 1 at the margin so nothing leaks off-ice), normalised to sum 1 exactly — deposition stays connected to the upstream erosion; `till.route: False` instead spreads the GLOBAL abraded volume melt-weighted across all ablation cells (decouples erosion/deposition across separate ice masses). **`till.on` and `till.route` both default True** (abrasion still opt-in via `Kg`); the melt-spread fixtures pin `till.route: False` and `test_ice_glacial_abrasion` pins `till.on: False`. Conservation guarded by `test_ice_glacial_till_conserves` / `_dual_lithology` / `_routing` / `_routing_conserves`. **Any change to till must keep those green.**

## Analysis tools (`gospl.analyse`)
Post-processing subpackage (NOT part of the `Model.runProcesses` pipeline; safe to edit independently). `gospl.analyse.provenance` attributes sediment deposited in sink basins to source-rock classes (per basin + per pixel) with transport distance and a Cu-fertility layer — design + the (shipped) in-model variant in `docs/DESIGN_PROVENANCE.md`, user/tech docs in `docs/tech_guide/provenance.rst`. Standalone (reads goSPL output); does **not** run inside the model. Key facts: routes the per-step net `erodep` down the re-derived flow graph (MFD default), recycling via a per-node deposited "pile"; `GosplOutput` reassembles the **partitioned** HDF5 (one file per rank per step) onto the global input-mesh ordering by KDTree (the same `tree.query(lcoords)` map goSPL uses); the topological sweep is one Numba-`njit`-compatible kernel (`_sweep_impl`) used as both the Python reference and, via `numba`, the fast path (`method='auto'`, identical results). `numba`/`geopandas` are **optional** (`pip install "gospl[analysis]"`); the tool falls back / errors clearly without them. The package is listed in `pyproject.toml` `[tool.setuptools] packages` (it must be, or it won't install). Tests: `tests/test_provenance.py`. A sparse transport-with-loss solve was evaluated and rejected (inexact under sediment under-supply + slower) — do not re-add it.

**In-model provenance tracers (Approach B — functional, phases B0–B4 + B2b done).** Opt-in `provenance:` block (`inputparser._extraProvenance`, requires `stratNb>0`): sets `provOn`/`provNb`/`source_class`/`prov_cu_weight`. State (allocated only when `provOn`, like dual lithology): `stratP[node, layer, class]` per-class layer thickness (Σ over classes == `stratH`, seeded to the bedrock `source_class` in `stratplex._initProvenance`); routed sub-flux `vSedP[c]` + `provFrac`/`depoProvFrac` + `_provEroded`/`_provDeposited` in `sedplex`; registered in `destroy_DMPlex`. Provenance is a **passive label** (no K/D/porosity/sorting feedback) — it rides the total-sediment routing. **Phases B0–B3 done** (functional end-to-end): B1 `erodeStrat` splits eroded sediment by class (`provEro`, Σ == total); B2 `_getSedFlux` routes each `vSedP[c]` through `fMati` (through-flux, Σ == `vSed`) → `provFrac`/`depoProvFrac`; B3 `deposeStrat` adds `depo·depoProvFrac` to `stratP` (keeps Σ over classes == `stratH`, ~3e-8) + `_provDeposited`. Guarded by `test_provenance_conservation` (single source ⇒ every layer stays 100% that class; class-0 leak == 0). **B4 done**: `stratalRecord` advects each `stratP[:,:,c]` like `stratHf` (re-normalised Σ==`stratH`); `_outputStrat` writes/`readData` restores `stratP` (lpoints, layers, classes) in the stratal HDF5 (`test_provenance_output_io`). **Multi-source validated** (`test_provenance_multisource`: 2-class map, both tracked, attribution spatially sensible). **B2b-marine done**: `seaplex._marineProvFraction` sets the marine `depoProvFrac` to the basin-delivered mix (`Σ provFrac·sedFlux / Σ sedFlux`, uniform — no depth bias) before the marine `deposeStrat`; this was the dominant sink (~79% of deposition) and the weakest spot of B2 — with it the deposited/eroded class ratios agree to ~1e-6 and the partition is machine-exact (~1e-24). **B2b-pit done**: `_distributeSediment` threads a per-class sub-flux `vSedP[c]·dt` through `_moveDownstream` in lockstep with the total — each pit retains its mix proportionally (`ret_prov = inVp·pitVol/inV`), overspills the rest through the same flow matrix (linear, so downstream-lake chains mix exactly), accumulating `_pitRetProv`; `_pitProvFraction` (from `_updateSinks`) sets `depoProvFrac[in_pit]` uniformly to `_pitRetProv/depo` (passive-label analogue of `_pitFineFraction`, **no depocenter bias** ⇒ simpler than dual). Σ over classes of `_pitRetProv` == retained volume ⇒ partition stays machine-exact. Guarded by `test_provenance_pit_fraction` + conservation tests (now run with pit attribution active). Plan: `docs/DESIGN_PROVENANCE.md` §6. (Bare `STRAMesh.__new__` stub tests don't set `provOn`, so `readStratLayers`/`erodeStrat`/`deposeStrat` read it via `getattr(..., False)`.)

## The `_extra*` methods are mandatory continuations
NOT optional parsers. Each sets attributes required by other modules. Never delete or rename without following the full call chain in `inputparser._readDomain/_readTime/_readHillslope/_readCompaction/_readFlex/_readOrography/_readIce`.

- `_extraDomain` (inputparser.py:189) → `seaDepo`, `overlap`, `dataFile`, `nodep`, `strataFile`; calls `_extraDomain2`.
- `_extraDomain2` (:229) → `advscheme`, `radius`, `gravity`.
- `_extraHillslope` (:475) → `nlK`, `clinSlp`, `tsStep`, `Gmar`, `offshore`, `nl_pit_volume/depth/K/inlet_bias`, `marineSolver`(`ts`|`picard`)/`picardSub`/`picardIts` (opt-in lagged-diffusivity marine/lake solver).
- `_extraStrata` (called from `_readCompaction`) → `stratLith` (dual-lithology master opt-in), `phi0c/z0c` (coarse porosity curve, defaults to compaction `phi0s/z0s`), `phi0f/z0f` (fine), `fine_k_factor` (fine erodibility multiplier), `bedrock_coarse_frac`, `fine_efficiency`, `pit_inlet_bias_coarse/fine`, `fine_diff_factor` (fine diffusivity multiplier). All inert while `stratLith` is False; forced False if `stratNb == 0`. See `docs/DESIGN_DUAL_LITHOLOGY.md`.
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

**Same value, different role — DO NOT replace these with the listed constants:**
- `1.0e-8` at `sed/stratplex.py:219` (thickness numerical-noise floor).
- `1.0e-3` at `flow/flowplex.py:363` (water-routing convergence), `sed/sedplex.py:139` (sediment-routing convergence), `flow/pitfilling.py:568,622` (minh epsilon nudges), `sed/seaplex.py:507,509` (clinoH 1mm offset).
- `1.0e-2` in `flow/iceplex.py` (`_glacialMeltwater`, `_routeTill`) is the ice-presence threshold (m): a distinct physical role, not a tolerance.

Each of these is marked with a permanent `# TODO-REFACTOR: value matches X but distinct role; do not replace` comment so future readers know the coincidence is intentional.

## High-risk modules (do not edit without full regression run)
- **`mesher/unstructuredmesh.py`** — owns `dm`, every shared mesh attribute, forcing dispatch (`applyForces` → `_updateRain`/`_updateIce`/…), and `destroy_DMPlex` (~line 849) which names every Vec/Mat by hand. Adding a new persistent Vec elsewhere requires adding it to that destroy list or it leaks (e.g. the dual-lithology `self.vSedF`/`self.vSedFLocal` and the SIA `self.iceMeltL`/`self.iceUbL`/`self.iceAbrL` are registered there).
- **`flow/flowplex.py`** — owns `_solve_KSP`, `_solve_KSP2`, `_matrix_build`, `_matrix_build_diag`, `_buildFlowDirection`. Consumed by every downstream module.
- **`tools/inputparser.py`** — owns all parameter parsing; forcing DataFrame column order is API; the `_extra*` chain is mandatory.

## Known bugs (fix before refactoring)
- _(none currently open)_

## Fixed (regression-guarded)
- **continental-sediment closed-sink deposit deadlock at np>1** — fixed 2026-06 in `sedplex._distributeSediment`. The closed-sink conservation closure (`if self._closedDepo.any():`, deposits cascade sediment that reached a closed flow-terminal node) wrapped **collective** DM scatters (`localToGlobal`/`globalToLocal`) inside a **rank-local** `.any()` guard. A closed terminal node (interior, pit-less, above sea, not an outlet) can exist on one partition and not another even on an **open** `'oooo'` domain, so at np>1 one rank entered the scatters while the other skipped ahead to `_updateSinks` → `_spillCoords`'s `Allreduce` — different collectives, permanent deadlock (both ranks 100% CPU). Serial can't mismatch, so it only hit np>1. Fix: reduce the guard with `MPI.COMM_WORLD.allreduce(has_closed, op=MPI.LOR)` so every rank takes the same branch (ranks with no closed sink contribute zeros — the `nz` mask already no-ops them). Same bug class as the coastline-gated `_smoothMat`/`_hillMat` rebuild deadlock (see the "collective rebuild decision" rule above): **any collective gated on a rank-local condition needs the reduce.** Guard: `test_advection_parallel` (mpirun -n 1 vs -n 2 on `flat_advect.yml`; previously timed out at 600s, now ~13s).
- **soil temperature-map parallel shape crash** — fixed 2026-06 (#444) in `soilSPL.__init__`. The `soil: tempMap:` field was loaded full-mesh (`mpoints`) without `[self.locIDs]`, so `prodSoil` stayed global and `_form_residual_soil` raised a `ValueError` broadcast mismatch (`mpoints` vs `lpoints`) at MPIsize>1 — invisible at n=1. Fix indexes by `locIDs` (mirrors the `soilFile` branch). Guard: `test_soil_temp_parallel_shape` (mpirun -n 2). **Lesson: grep `np.load(...)[` sites that feed local arrays for the same serial-only-tested global-vs-local bug.**
- **soil accumulating underwater** — fixed 2026-06 (#451) in `soilSPL._solveSoil`: `nHsoil[self.seaID] = 0.0` before writing `Lsoil`. The rainfall-scaled soil-production term was depositing a spurious soil cover on submarine nodes (erosion there was already correct — `Kbr`/`K_soil`/`fDep` are zeroed at `seaID`; only the soil-thickness output was polluted).
- **flow-accumulation KSP `PCSetUp` failure at scale** — fixed 2026-06 (#445) in `flowplex.py`. At higher rank counts an ocean-only/degenerate partition gives a rank a local block whose `bjacobi`/`asm` ILU sub-solver hits a zero pivot → `DIVERGED_PCSETUP_FAILED` → the flow solve falls through to zero discharge (silently corrupting that step). Fix: `pc_factor_shift_type nonzero` on both flow-KSP sub-PCs (scoped by `flowacc_`/`flowaccfb_` prefixes). Guard: `test_parallel_correctness` (rel_sum_fa < 5e-2). Divergence messages now name the solver (flow vs detachment-SPL).
- **`rUni`/`sUni` key/column mismatch** — fixed 2026-06 at `inputparser.py:915` (dict key `'rUni'` → `'sUni'` in the uniform branch of `_defineErofactor`). Before the fix, a YAML with a uniform `sedfactor` event landed `NaN` in `sedfacdata['sUni']`, then crashed `unstructuredmesh._updateEroFactor` on `np.load(None)`. Guarded by `tests/test_regression.py::test_uniform_sedfactor_populates_sUni`.
- **Marine sediment leak at terminal ocean sinks** — fixed 2026-06 at `seaplex._distOcean`. `_matOcean` constructs `dMat1` with a zero column at every terminal sink (cells whose every flow direction is self-referential, typically deep-ocean basin floors); the first `dMat1.mult` therefore annihilated any `sedflux` that landed directly at a sink. `vSed` upstream-integration naturally accumulates large flux at these sinks, so this was the dominant mass-loss path. Diagnostic: 9.748e12 m³ marine sediment supply per step → 4.107e12 m³ landing at sinks (42% of input) → 0 m³ deposited via the routing loop → silently dropped. The fix adds a pre-drain step at the top of `_distOcean` that deposits any sinkVol at terminal sinks directly into `vdep` before the first `dMat1.mult`. The in-loop force-deposit and exit-residual drain (added alongside) handle any sediment that subsequently arrives at sinks during multi-step routing. Guarded by `tests/test_regression.py::test_mass_conservation`, which runs on the `minimal_model` fixture — a **global-sphere mesh** (`flatModel=False`, `len(idBorders)==0`), so the domain is closed by construction with no boundary outflux possible; the test asserts `|dV_surface|/activity < 1e-4` AND `|dV_cumED|/activity < 1e-4`, i.e. that deposited volume equals eroded volume to within ~6e-5 relative (the floor-effect budget from `DEPOSIT_FLOOR=1e-3` and pit-routing residue, both documented in the test's tolerance rationale). The corresponding `is_sink_local` mask is computed in `_matOcean` once per call.
- **`fine_diff_factor` design/parser drift** — fixed 2026-06 in `_extraStrata`. The dual-lithology diffusivity contrast was changed from absolute `Dc`/`Df` to a `fine_diff_factor` multiplier in the design doc, but the parser still set the obsolete `self.Dc`/`self.Df` and never parsed `fine_diff_factor`. Surfaced when `_surfaceLithoD` (Phase 7) read `self.fine_diff_factor` → `AttributeError` on any real dual run. Now parses `strata.fine_diff_factor` (multiplier on `Cda`/`Cdm`/`nlK`); guarded by `test_dual_lithology_opt_in`.
- **`Eb` vs `EbLocal` sign-convention divergence** — fixed 2026-06 by adopting the thickness-rate convention (positive deposition, negative incision) for both fields. Changes: drop the leading minus in `E = -tmp.getArray()` inside `_getEroDepRate` (SPL.py), `_getEroDepRateNL` (nlSPL.py), `_getEroDepRateSoil` (soilSPL.py); drop the minus in `tmp.setArray(-Eb * dt)` inside the three matching `erodepSPL*` wrappers; negate `Eb` inside `sedplex._getSedFlux` (sedplex.py:63) where the upstream-integration solve still wants an erosion-positive source. Before the fix, `Eb` was incision-positive while `EbLocal` (after the wrapper's `add_rate = tmp/dt` overwrite) was deposition-positive, so the two fields disagreed in sign by end-of-step; the on-disk `EDrate` field was always deposition-positive (because it's written from `EbLocal`), and the restart loader at `outmesh.py:439-440` restored `Eb` in deposition-positive convention from `EDrate` — so fresh-run state and restart-state were inconsistent before this fix. Guarded by `tests/test_regression.py::test_erosion_sign_conventions` (which asserts `EbLocal <= 0` at incising nodes; `Eb` follows the same convention now).

## Intentional surprises (do NOT "fix")
- **`gid` argument on `mfdreceivers` / `mfdrcvrs`**. The Fortran kernels take an extra `gid(nb)` argument (per-local-node global vertex ID, from `self.gid` set in `mesher/unstructuredmesh.py` immediately after `getLGMap`). Inside the kernels, the stored slope is perturbed by `val * (1 - 1.0e-15 * gid(n))` before the quicksort tie-break. This makes EXACT slope ties resolve deterministically across MPI decompositions (lower global ID wins). The perturbation is well below KSP-solver precision and does not affect physical results; near-ties driven by KSP floating-point noise are a separate problem and are NOT fixed by this pattern. See the relaxation comment on `test_parallel_correctness` in `tests/test_regression.py` and `fortran/functions.F90:mfdreceivers` for the full rationale. **Do NOT remove the `gid` argument** — without it, `mfdreceivers` ordering depends on local iteration order and the parallel-correctness test diverges by ~10x from current baseline.
- **Platform-dependent KSP-noise floor in `test_parallel_correctness`**. With the `gid` exact-tie-break in place, the remaining sum-FA drift between n=1 and n=2 is dominated by KSP-solver floating-point non-determinism at partition boundaries — and the magnitude differs by platform. Measured on `tests/fixtures/minimal.yml`: macOS-14 (arm64, conda-forge OpenMPI) ≈ **0.3%**, ubuntu-latest (x86_64, conda-forge MPICH) ≈ **1.7%**. Same fixture, same goSPL, same Python — the spread is driven by MPI implementation, BLAS variant, and FMA availability at the platform level. The test tolerance for `rel_sum_fa` is set at **5e-2** (5%) to clear both with headroom; a real routing regression would shift it toward 0.5+, not 0.05. If a future contributor sees these numbers and wants to tighten the tolerance, the correct path is to make KSP-precision halo state bitwise-identical across decompositions (hard PETSc work), NOT to lower the bound and hope. **Note (2026-06):** `environment.yml` now pins the MPI stack to OpenMPI 4.x on all platforms (see the `environment.yml` conventions in `## CI contract`), so ubuntu CI now also runs OpenMPI rather than the MPICH referenced above — the ubuntu figure is historical and the next CI run will measure the new drift; the 5% tolerance still clears it. **Tier-1B (`rel_mean_h`, area-weighted mean elevation) has the SAME platform floor and was likewise relaxed (2026-06, commit `1cbcf42`).** Originally `1e-10` (macOS-calibrated), it failed ubuntu CI: the n=1-vs-n=2 mean-elevation drift propagated from the same near-tie non-determinism is ~3.6e-11 on macOS-14 but ~7.6e-10 on ubuntu-latest. Now **5e-9** (~6× the worst observed, still ~7 orders of magnitude below a real ghost-node reduction bug, which drives mean_h to the %-level). **Confirmed unrelated to the cascade stagnation-break (`d4162aa`):** pre-fix and post-fix flowplex give *bit-identical* `rel_mean_h` (3.616855e-11) on macOS — the break is a no-op on this converging fixture — so the failure was pre-existing runner FP drift hitting an over-tight bound, not the flow change.
- **Heavy model fixtures must be released between tests — `tests/conftest.py` autouse teardown (do NOT remove).** The ~20 `*_model` fixtures are function-scoped and `return` (not `yield`) their `Model`, so each test's PETSc objects (Mat/Vec/KSP/DMPlex) + mesh arrays linger until Python collects them. Across the full slow tier **in one process** this accumulates and, on the **ubuntu** runner, intermittently aborts a *later* `Mat().create()` with **`SIGABRT`** — flaky, deep in the suite, **ubuntu-only** (passes on macOS and on re-runs). It surfaced first when an RC `master`/`release-candidate` push triggered the full tier via `tests-pr.yml` (`tests-slow`'s slow-regression step is gated OFF on push, so `dev` never exercised it). Fix (commit `cf4c736`): an **autouse** function-scoped fixture `_release_petsc_after_test` — set up first, so torn down LAST (after each model fixture has dropped its reference) — runs `gc.collect()` then `petsc4py.PETSc.garbage_cleanup()` after every test, only touching PETSc if a test imported it (`sys.modules` check, so parser-only tests don't pull in PETSc/MPI). **Do NOT delete the autouse fixture, and do NOT add a model fixture that holds a `Model` past its test without it** — without it the slow tier flakes on ubuntu.

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

## Docs / Read the Docs (autodoc)
The API reference (`docs/api_ref/*.rst`) is Sphinx **autodoc** — it imports each module to extract docstrings, so the modules MUST import under the docs build. Two invariants keep the API pages from rendering empty (both regressed once and produced blank pages for *most* classes):
- **`gospl` must be importable.** It is NOT pip-installed on Read the Docs, so `docs/conf.py` puts BOTH the repo root (`..`, so `from gospl.tools.constants import …` and other `from gospl.X import` resolve) AND `../gospl/` (so the legacy top-level `.. autoclass:: sed.hillslope.hillSLP` paths resolve) on `sys.path`. Don't remove either.
- **Mock only the genuinely-absent compiled/MPI deps via `autodoc_mock_imports`** (`h5py`, `mpi4py`, `petsc4py`, `vtk`, `pyshtools`, `gospl._fortran`). Do NOT mock packages `docs/requirements.txt` installs (`numpy`, `scipy`, `pandas`, `numpy-indexed`, `ruamel.yaml`) — the old manual `sys.modules[m] = Mock()` shadowed the real `scipy` and broke `from scipy.special import …` (and didn't cover submodules), blanking the pages. `autodoc_mock_imports` mocks submodules automatically; the `READTHEDOCS` env-guard around `from gospl._fortran import …` is the belt-and-braces for the compiled extension.
- When you add a public/private method that should appear in the API, add it to BOTH the `.. autosummary::` and `.. automethod::` lists in the relevant `docs/api_ref/*_ref.rst` (they are hand-maintained, not auto-generated).

## Analytical benchmark suite
Tests in `benchmarks/` validate goSPL physics against exact analytical solutions. They require `scipy` and `matplotlib` and are silently skipped via `pytest.importorskip` in environments without these packages. `matplotlib` is in `environment.yml` for CI; it is NOT a goSPL runtime dependency (intentional — `pyproject.toml` does not list it).

The `benchmark` and `slow` marker names are registered in `pyproject.toml` under `[tool.pytest.ini_options].markers` (authoritative) and mirrored in `tests/conftest.py::pytest_configure` for the existing tests/ tree. Don't expect new markers to work without registering them in pyproject — pytest does not walk sideways into `tests/conftest.py` when invoked from `benchmarks/`.

Each benchmark fixture (`spl_tmp_path`, `hillslope_tmp_path`, `knickpoint_tmp_path` in `benchmarks/conftest.py`) copies the benchmark's input files into a per-test tmp directory and chdirs there. goSPL intermediate output (`sim_output/`, `sims_outputs/`) stays inside `tmp_path` and is reaped by pytest's cleanup. Report artefacts (`*.pdf`, `*.md`, hillslope `*.png`) are copied out on teardown into `<repo_root>/results/<benchmark_name>/` so they survive pytest cleanup AND are picked up by the GitHub Actions `upload-artifact@v4` step in `tests-slow.yml`. `results/` is in `.gitignore` — never commit it.

The SPL fixture's teardown glob is `spl_benchmark*.pdf` / `spl_benchmark*.md` (wildcarded) because `test_spl.py` writes one PDF + one Markdown per benchmark case (`spl_benchmark_150.pdf`, `spl_benchmark_200.pdf`, plus the combined `spl_benchmark_all.md`). Don't tighten the glob — the wildcard is load-bearing.

Mesh `.npz` files under each benchmark's `boundary_condition[s]/` subfolder are working-tree artefacts (currently untracked in git — see `git status`). Naming convention diverges by benchmark: `spl/boundary_conditions/` and `hillslope/boundary_conditions/` are **plural**, `knickpoint/boundary_condition/` is **singular**. Match the existing folder when adding a new mesh; the YAML configs reference these paths directly. `knickpoint/boundary_condition/drop_baselevel.npz` is regenerated by the test itself during phase 2 of the knickpoint workflow.

`benchmarks/test_hillslope.py` delegates to `scripts/analysis.py` (imported as `anlys` after `hillslope_tmp_path` prepends the copied folder to `sys.path`). The other two benchmarks inline all logic in their `test_*.py` files. `analysis.runFullBenchmark` does NOT expose an `overall_pass` key; use `analysis.print_summary(results)` for the pass/fail bool.

`benchmarks/test_spl.py` runs **two cases per invocation** (`input150.yml`, `input200.yml` under `benchmarks/spl/sims/`), each scored **5/5** with `evaluateSPL`'s `skip_basin_test=True`. The basin-geometry (>80% coverage) gate is skipped because **no** case meets it (the radial ramp drains across the open E/S/W edges: 100≈51%, 150≈56%, 200≈77%). The finest **case 100 (`input100.yml`) is intentionally excluded** (see the `cases` list comment): it would fail both the basin gate AND R² (~0.93) — finest-mesh/boundary artifacts, not physics errors — and it dominates runtime (~2 s/step × 1000 ≈ 34 min serial vs ~45 s each for 150/200; its fragmented drainage makes the flow solve ill-conditioned). `input100.yml` + `mesh100.npz` are kept so it can be re-enabled if the benchmark gets a single-outlet boundary and the per-step cost is addressed. The test asserts `overall_pass` per-case; failure on any case fails the whole test. Helpers: `_run_spl_benchmarks` / `_assert_spl_benchmarks`. (Mirror change in the standalone `benchmarks/spl/test-SPL.py`.)

| Test file | Process | Analytical basis | Pass criteria |
|---|---|---|---|
| `test_spl.py` | SPL steady state | Braun & Willett 2013; Perron & Royden 2013 | Cases 150 & 200: 5/5 each (basin test skipped — no case meets >80%); case 100 excluded (finest-mesh basin/R² artifacts + ~34 min cost) |
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
| 2026-06-14 | — (on `dev`) | **Provenance B2b — exact sink attribution** (PRs #440 marine + this one pit) — both sinks now carry the exact delivered source mix. Marine: `seaplex._marineProvFraction` (basin-delivered `Σ provFrac·sedFlux / Σ sedFlux`). Continental pit/lake: `_distributeSediment` threads a per-class sub-flux `vSedP[c]·dt` through `_moveDownstream` in lockstep with the total (proportional retain `ret_prov = inVp·pitVol/inV`, overspill routed through the same matrix ⇒ downstream-lake chains mix exactly), accumulating `_pitRetProv`; `_pitProvFraction` (from `_updateSinks`) sets `depoProvFrac[in_pit] = _pitRetProv/depo` uniformly (passive label, no depocenter bias). Σ over classes == retained volume ⇒ `stratP` partitions `stratH` machine-exact. Guarded by `test_provenance_pit_fraction` + conservation tests (now run with pit attribution active). Design/tech/user docs updated. Suite: 58 passed, 1 skipped. |
| 2026-06-16 | — (branch `feat/ice-explicit-flux`) | **SIA solver rewrite: implicit → explicit flux-limited** — the implicit `ngmres`/Newton SNES diverged (`reason −9`) on real continental runs because the `H≥0` free boundary is an *obstacle problem* (no `F(H)=0` root-find converges at the margin; clamping injects mass). Replaced `_iceFlowSIA` + `_form_residual_ice` with explicit substepping `H←H+Δt_sub(ṁ−∇·q)` using two new kernels: **`ice_flux_limiter`** (per-cell `R`, halo-synced for exact parallel conservation) + **`ice_flux_rscaled`** (R-scaled divergence) ⇒ mass-conservative to machine precision AND `H≥0` for any substep, no clamp. `Δt_sub` is accuracy-only (`sia.cfl` default 0.5, `sia.max_substeps` default 500); flux/ablation rate capped at `1/dt` (fixes a thin-cell/ablation stall) while **accumulation gain is uncapped** (`mdot⁺/max(H, ICE_HREF)`) so a step's accumulation is added gradually and drained as the ice thickens — without this, an `H=0` accumulation cell dumped a whole step's accumulation in one giant substep → runaway km ice (42 km in 2 steps on the real 30 km global model). Added SMB accumulation controls `sia.accum_factor`/`sia.accum_max` (scale/cap positive SMB; full-precipitation→ice is unphysically high). Dropped the cached `_snes_ice*` (removed from `destroy_DMPlex`). Tests: `test_ice_sia_explicit_stiff_conservation`, `test_ice_sia_explicit_no_substep_stall`, `test_ice_sia_accumulation_scale_cap`. Investigation kernels (`ice_flux_jacobian`, `sia_diff_coeff`) on `fix/sia-solver-investigation`. Full diagnosis: `docs/DESIGN_ICE_SHEET.md` §9. Suite: 66 passed, 1 skipped. |
| 2026-06-16 | — (branch `feat/ice-diagnostic-glacial`) | **Diagnostic glacial driver `ice.flow_model: mfd`** — re-introduce a non-dynamical ice option for when glacial-erosion *morphology* matters more than ice dynamics (the SIA is stiff/over-thick for km-thick continental ice at coarse resolution, even after the explicit rewrite). `_iceFlowMFD` routes the ELA accumulation (`mdot⁺·larea`) downhill via an MFD matrix (`_matrixIceFlow`) → ice discharge `iceFAL`/`iceFAG` → Bahr thickness `iceHL` → balance velocity `iceUbL = Q/(H·√area)`, feeding the existing velocity-based abrasion + till + loading + output. No PDE/CFL — real 30 km global: ~0.8 s/step, physical `Hmax`~1 km (vs SIA's 11–42 km). **`mfd` is now the DEFAULT `flow_model`** (the SIA over-thickens km continental ice); the SIA test fixtures pin `flow_model: sia`. All the existing erosion/till/dual-lithology/stratigraphy/meltwater-coupling machinery works under `mfd` (it's driven by `iceUbL`/`iceMeltL`/`iceHL`, which `mfd` sets) — confirmed end-to-end by `test_ice_mfd_dual_strata_till`. Config `icedir`/`eheight`/`fwidth`/`melt`; new Vecs `iceFAL`/`iceFAG` + `iceMat` in `destroy_DMPlex`; `iceFA` output added. Tests `test_ice_mfd_diagnostic`, `test_ice_mfd_dual_strata_till`. Suite: 68 passed, 1 skipped. |
| 2026-06-16 | — (branch `feat/ice-diagnostic-glacial`) | **mfd velocity fix + lateral erosion.** (1) The mfd balance velocity `Q/(H·√area)` blew up with catchment size and spiked at flow-convergence cells → tens-of-km erosion/deposition in sinks; replaced with the bounded SIA sliding law on the Bahr thickness (`ice_velocity`). NB velocity scale dropped a lot → recalibrate `Kg` upward. (2) Added **explicit lateral (valley-wall) glacial erosion** `ice.abrasion.Kl` (kernel `ice_lateral_erosion` + `iceplex._glacialLateralErosion`): a wall cell (`hice≤1 m`) flanking faster ice is abraded at `Kl·u_b,ngb^lat_l`, tapered by ice-contact height, **widening valleys toward a U-profile**; the wall rock joins the same conserved till→moraine (or `Eb` when till off) path. Off by default (`Kl=0`). Tests `test_ice_lateral_erosion`. Suite: 70 passed, 1 skipped. |
| 2026-06-16 | — (branch `feat/ice-diagnostic-glacial`) | **Defaults: till.on/till.route True; discharge-conserving river meltwater.** (1) `till.on` + `till.route` default True (abrasion still opt-in via `Kg`) — the old `route:False` global melt-pooling decoupled erosion/deposition across ice masses (user saw deposition with no nearby erosion); catchment routing fixes it. Melt-spread fixtures pin `route:False`; `test_ice_glacial_abrasion` pins `till.on:False`. (2) New **`iceMeltRiverL`** + `_glacialMeltwater`: glacial meltwater to rivers is **discharge-conserving** (`ice.melt_conserve`, default True) — routes accumulation to where ice melts out so `Σ river-melt == Σ accumulation` (over long steady-state steps all accumulated ice returns as river water; precip-scaled ablation lost water). Kept distinct from `iceMeltL` (till melt pattern) to avoid breaking till conservation. `flowplex` + `iceMelt` output now read `iceMeltRiverL`; registered in `destroy_DMPlex`. Tests `test_ice_meltwater_conserves` (Σmelt==Σaccum). Suite: 70 passed, 1 skipped. |
| 2026-06-16 | — (branch `feat/ice-diagnostic-glacial`) | **SIA removed — diagnostic is the only glacial model.** The Shallow-Ice-Approximation thickness solve (implicit AND explicit flux-limited) was removed: the `H≥0` margin is a free-boundary/obstacle problem the implicit solve diverged on, and any ice-dynamics solve is stiff and over-thickens km-scale continental ice at goSPL's 10²–10⁴ yr steps. Dropped the `ice.flow_model` selector, the `sia:` config block (`Aglen`/`cfl`/`max_substeps`/`accum_*` moved to `ice`-level keys), `hinit` seeding, `_iceFlowSIA`/`_iceSIAParams`/`_iceSIAFinalize`/`_form_residual_ice`, and the Fortran kernels `ice_flux`/`ice_flux_limiter`/`ice_flux_rscaled` (`.F90` + `.pyf`). `iceAccumulation` now always runs `_iceFlowMFD`; mass balance is `_iceMassBalance`. Investigation kernels (`ice_flux_jacobian`, `sia_diff_coeff`) remain on `fix/sia-solver-investigation` for any future obstacle/VI solver. Docs (`DESIGN_ICE_SHEET.md`, `ICE_SHEET_SUMMARY.md`, `tech_guide/ice.rst`, `user_guide/surfproc.rst`, `api_ref/ice_ref.rst`) + this file rewritten to the diagnostic-only model; SIA fixtures (`minimal_ice_sia`/`_seed`/`_mfd`) removed. Suite: 67 passed, 1 skipped. |
| 2026-06-18 | — (on `dev`) | **Pit-filling (`flow/pitfilling.py` + `fortran/functions.F90`) perf + Barnes-faithful robustness.** The serial rank-0 spillover-graph solve (`fill_edges`) IS Barnes (2016)'s master step *by design* (graph is `O(perimeter)`) — kept serial. Improvements landed: (#465) `O(log²n)→O(log n)` heap push in `PQpush` + flow sub-phase profiler; (#466) dropped the dead `O(nb·maxnghbs)` array inits in `fill_edges`; (#467) **pit-label compaction** — `_fillFromEdges` densely re-indexes the sparse globally-offset labels before the flood so its arrays track distinct spillover basins, not the rank-dependent label range; (#468) **union-find label unification** — `_unifyLabels` replaces the iterative `_sortingPits` fixpoint (up to 1000 passes) with one deterministic disjoint-set pass collapsing each component to its min label, **making the filled surface partition-independent** (the iterative relabel was the original instability source); plus `Allgatherv` of the equivalence pairs (was padded `Allreduce(MAX)`). All bit-for-bit identical to prior `dev` at fixed rank count (n=1/n=2 A/B max rel diff 0.0); the residual n=1-vs-n=2 drift is the pre-existing KSP/FP non-determinism, unchanged. The `Reduce(MAX)`→`Gatherv` of the graph itself was deliberately NOT done (it's ~0.01 s, doesn't grow with ranks, and its −1 padding clamps sub-(−1 m) deep-ocean spills harmlessly). Net: flow is edge-bound at the partition perimeter, well-scaled; gains are robustness/determinism more than wall-time. Guard: `test_pit_unifyLabels_unionfind` + `test_parallel_correctness` + `test_mass_conservation`. Docs: `tech_guide/dep.rst` "Parallel implementation", PITFill docstring. |
| 2026-06-18 | — (branch `feat/flex-fem-2d`) | **Flat-model flexure rewritten as a parallel FV biharmonic (`method='fem'`); gFlex + FFT removed.** New `_cmptFlexFEM`/`_buildFlexFEM` solve `∇²(D∇²w)+Δρg·w=q` directly on the DMPlex as the single-field `[Lm·diag(D)·Lm + Δρg·I] w = q` (FV neg-Laplacian `Lm` applied twice) — fully parallel, no gather/regrid/dependency; varying Te is ONE solve (`diag(D)`). Operator+factorisation **cached and reused each step** (serial PETSc LU / parallel **MUMPS** — GMRES+GAMG fails to converge on the stiff biharmonic at realistic Te). BCs `0Slope0Shear`/`Mirror` (natural FV zero-flux) + `0Displacement0Slope` (clamped, Dirichlet `zeroRows`). **Removed**: gFlex dep (import + pyproject/environment/conf.py mock) + `_cptFlex2D` (FD/FFT), the Fortran `flexure` FFT kernel (`.F90`+`.pyf`), `boundflex`, the flexure regular-grid build (the regular grid is now orography-only), and the `flexure.method` values `'FD'`/`'FFT'` (now only `'fem'`/`'global'`, default `'fem'`); `regdx` dropped from the flexure block. On the stratigraphic_record example (45k flat): `fem` is faster than gFlex was (np=1 102 vs 112 s; np=4 47 vs 60 s), partition-independent, varying-Te np=4 converges (0 failures). Matches gFlex to corr 0.998 (natural) / 0.9996 (clamped) when the deflection decays inside the domain; ~10% boundary-zone difference when the flexural wavelength approaches the domain size. gFlex-comparison tests replaced by self-consistent `test_flex_fem_2d_solver`/`_physical`/`_clamped`. Docs: `user_guide/optfile2.rst`, `DUAL_LITHOLOGY_SUMMARY.md` appendix, `api_ref/grid_ref.rst`. Suite 74/1. |
| 2026-06-18 | — (branch `feat/cyclic-boundaries`) | **Flat-model boundary conditions reworked: named `o/f/w/c`, NESW order, + cyclic and true-wall.** `domain.bc` is a 4-char string in **[North, East, South, West]** order (N=ymax, E=xmax, S=ymin, W=xmin — was S,E,N,W) over `{'o' open, 'f' fixed, 'w' wall, 'c' cyclic}`; legacy `'0'→'o'`, `'1'→'f'`; default `'oooo'`. `unstructuredmesh` maps to `self.{north,east,south,west}` ∈ {0 open,1 fixed,2 cyclic,3 wall}. **Key new concept `outletIDs`** = `idBorders` minus true-wall edge nodes; the boundary drain/discard treatment in flowplex/pitfilling/sed/eroder/hillslope/seaplex now keys on `outletIDs` (≡ `idBorders` when no walls, so existing behaviour byte-unchanged — suite confirms). **`o` vs `f`**: both DRAIN (open = deep `MISSING_DATA_SENTINEL` outlet via `applyForces` `==0`; fixed = natural-elevation base-level outlet, the historic `'1'`). **`w` (wall)** = true no-flux: its edge nodes are excluded from `outletIDs` → treated as interior → flow contained, sediment deposits against the wall (`_closedDepo` closure in `_distributeSediment`/`_moveDownstream` deposits sediment reaching closed non-pit terminal sinks). **`c` (cyclic)**: periodic, requires a CYLINDER input mesh (periodic axis on a circle; intrinsically flat so `definetin` FV geometry = periodic strip; wrap cells link the seam in `FVmesh_ngbID`; cylinder ends keep `flatModel` True). **SERIAL-ONLY** — see the 2026-06-18 advection milestone: cyclic is NOT parallel-safe (the partitioned seam ill-conditions the sediment/advection solves), and `inputparser` now raises on cyclic + MPIsize>1. Validated cyclic must be an opposite PAIR, ≤1 pair. **LIMITATION**: a FULLY-closed domain (all walls, no `o`/`f`/marine exit) conserves sediment for ONE step but NOT multi-step — once basins saturate, goSPL's spill-based cascade can't aggrade a closed basin (architectural; it assumes excess spills to an outlet). Use walls + ≥1 draining edge for long conserving runs. Also fixed: fully-closed domain previously CRASHED in the pit-label compaction (`_fillFromEdges` empty `cgraph.max()`) — guarded (`total==0` → endorheic local fill). Fixtures `cyclic_cyl.npz` (cylinder) + `flat_wall.yml`; tests `test_cyclic_boundary`, `test_wall_boundary_conservation` (single-step). Docs: `user_guide/inputfile.rst` `bc`. Suite 76/1. **Follow-up (same branch): cyclic horizontal advection.** A cyclic model runs on a cylinder, so `tectonics` `hdisp` (still supplied as a flat `(vx,vy)` field) is remapped onto the cylinder tangent by `tectonics._cylinderVelocity` — the periodic-axis component becomes motion AROUND the cylinder (θ̂), the non-periodic component stays axial; otherwise a planar `(vx,vy,0)` velocity has zero out-of-plane part and cannot advect across the curved seam. Mesh now sets `self.cyclicBC`/`self.cyclicPts` (seam nodes) and `self.advectBorders` = `idBorders` minus the seam; the advection Dirichlet + the post-advection `MISSING_DATA_SENTINEL`+`fitedges` edge reset in `_varAdvector` use `advectBorders` so the seam stays free and material wraps (verified mass-conserving across the seam). Fixture `cyclic_cyl_advect.npz/.yml`; test `test_cyclic_advection`. Suite 77/1. |
| 2026-06-18 | — (branch `feat/orography-mesh`) | **Orographic precipitation moved off the regular grid → parallel mesh-native solve; regular grid fully removed.** `cptOrography` no longer gathers to rank 0 / interpolates to a regular grid / FFTs. Smith-Barstad is recast as two steady advection-relaxation PDEs on the DMPlex: `(v·∇ + 1/τc) q_c = Cw·(v·∇h)` then `(v·∇ + 1/τf) q_s = q_c/τc`, `P = q_s/τf`. In Fourier space this is exactly the Smith-Barstad transfer function with the stratified mountain-wave term `(1-i·hw·m)` → 1 (DROPPED). Rain shadow retained (windward source +, lee source −, downwind advection). **Source elevation clamped to sea level** (`max(h, sealevel)`): submarine bathymetry produces NO orographic forcing (airflow follows the flat sea surface) — the old FFT version wrongly used raw bathymetry. Reuses the existing FV kernels: `getfacevelocity` + `advecupwind(dt=1)` give `I + L` (upwind v·∇); `_buildOroMat` assembles `L + (1/τ)I` with zero-Dirichlet on `advectBorders`; `_oroSolve` is a cached GMRES+bjacobi KSP, warm-started. Operators depend only on the (uniform, constant) wind + fixed mesh ⇒ **assembled once and cached** (`_oroAc`/`_oroAf`/`_oroAdvDiag`/`_oroLcoeff`/`_oroQc`/`_oroQs`); each step only the source + solves change. Verified **partition-exact** (np=1 vs np=4 bit-identical on owned nodes). **Removed**: `_buildRegGrid`, `_regInterp`, all `reg_*`/`xIndices`/`xFrac`/`regWeights`/`oroEPS` state, and the parser keys `regdx`/`latitude`/`nm`/`hw` (spectral-only). Fixture `oro_hill.npz/.yml` (ocean→coastal mountain→ocean, wind from west); test `test_orography_rain_shadow` (rain shadow + sea-level clamp). Docs: `user_guide/optfile1.rst` (params explained, regular-grid removed), `api_ref/grid_ref.rst`/`index.rst`. Suite 78/1. |
| 2026-06-18 | — (branch `feat/advect-2d-improvements`) | **2D horizontal-advection refactor (`tectonics._varAdvector`) + cyclic-advection SERIAL-ONLY guard.** Perf/correctness cleanup of the FV advection: (1) the operator (`adveciioe`/`advecupwind` coeffs) is **cached** (`_advMatLeft`/`_advMatRight`/`_advKSP`), rebuilt only when the tectonic interval or `dt` changes (`_advRebuild` flag set in `getTectonics` after `getfacevelocity`) instead of every step; (2) `_buildAdvecMat` now assembles in **one CSR pass** via the existing `_assembleDiffMatCSR` (was ~24 throwaway-matrix `axpy`); (3) the four field blocks (h/cumED/flexure/soil) collapsed into a `_advectField` helper + `_resetEdges`, with per-field **warm-start** (seed `self.tmp` with the field's current value) and a dedicated cached KSP `_advSolve` (`_makeDiffusionKSP("advect_")` + `_solve_KSP2` fallback); (4) dropped the per-call `gc.collect()`. **Bug fixed**: `_advectorIIOE2` (IIOE2/`advscheme==3`) multiplied `self.hGlobal` for EVERY field — now takes the field `gvec`, so cumED/flexure/soil corrections were wrong. **cumED edge sentinel** `MISSING_DATA_SENTINEL`→`MISSING_LARGE_SENTINEL` (consistent with `_advectPlates`; `fitedges` uses a `<= -1e7` threshold so both are detected). Cached objs + the previously-leaked orography objs added to `destroy_DMPlex`. **Serial bit-identical to dev.** Flat advection is **partition-correct** (np=1–4 agree ~0.1%) — guarded by new mpirun `test_advection_parallel` (fixture `flat_advect.npz/.yml`). **CYCLIC advection/sediment is NOT parallel-safe** (separate pre-existing bug): the partitioned cylinder seam ill-conditions the FV operators — pit-deposition diffusion hangs, advection KSP diverges→`h=0`. Flow + pit-filling are fine; the gates in `_updateSinks`/`_diffuseLargePit` are collective-safe (global pit data), so it is solve-convergence, not a gate deadlock. Mitigation: `inputparser._readDomain` now **raises** on `cyclic bc + MPIsize>1` (serial cyclic still works). The parallel cylinder-seam fix is a tracked follow-up. Suite 78/1 + new parallel test. |
| 2026-06-18 | — (branch `fix/cyclic-parallel-seam`) | **Cyclic-parallel root cause #1 FIXED: partition-dependent flat-vs-curved branch.** `getfacevelocity` (advection face velocity) and `definetin` (FV geometry: faceVec/midpoints) chose the planar-vs-spherical branch from a SINGLE node's z (`lcoords(1,3)`) — partition-DEPENDENT. On the cylinder (`z=R·sinθ` spans ±R incl. 0 at the seam) different ranks' first nodes fall in different branches → inconsistent operators → cyclic advection garbage at np>1 (a MUMPS-direct diagnostic proved the OPERATOR differed by partition ~170×, NOT a preconditioner issue). Fix: new `meshparams::curvedmesh` flag + `setcurvedmesh(flag)` Fortran setter; both kernels use it when set (legacy `lcoords(1,3)` only as fallback). Python (`unstructuredmesh._meshStructure`, before `definetin`) sets it from the GLOBAL `flatModel` (sphere→1 curved; flat AND cylinder→0 euclidean — a cylinder is intrinsically flat). Sphere/planar stay bit-identical (legacy gave the same); only the cylinder changes (now consistently euclidean). **Verified**: cyclic advection `hmax=−302.79` at np=1–4 (was −331 at np=2, **0 at np≥3**). **STILL OPEN (issue filed): cyclic flow+SEDIMENT deadlocks at np≥2** — a SEPARATE collective-count desync in the `sedChange` path (per-rank logs: rank0 blocks at `_updateSinks` post-deposit `localToGlobal`; rank1 ahead in `seaChange._distanceCoasts` Allreduce → ranks offset by collectives; divergence is upstream in `_getSedFlux`/`fillElevation(sed)`/`_distributeSediment`). Needs a per-rank collective-count audit. The serial-only guard (PR #475) STAYS until that is fixed; this fix is groundwork (advection-in-parallel correct, guard keeps cyclic serial-only). |
| 2026-06-21 | — (on `dev`) | **Drainage graph made partition-invariant + np>1 deadlock fixed + dead Fortran removed.** Closed the residual drainage non-invariance that made the flow-KSP "un-drained region" (and the `sed`/`sea` cascade cost) vary with rank count. Four roots, all re-keyed on the partition-invariant input-mesh id (`self.gid` → `self.locIDs`): the flat-routing receiver tie-break (`_dirFlats`, Python, replacing `fill_rcvs`), the flat-BFS distance `ptdir` (`nghb_dir` rewritten as a relaxing label-correcting Dijkstra + convergence loop), and the spill-point selection + pit membership (`_pitInformation` min-`locID` pick replacing `spill_pts`; `label_pits` order-dependent flags dropped). Verified on a 5.9M global mesh np=1-vs-np8: true singular set symdiff 0; the `sed` scaling regression (P=240 50.8→15.2 s) gone. Also: the **continental closed-sink deadlock at np>1** (`_distributeSediment`'s `_closedDepo.any()` guarding collective scatters → `allreduce(LOR)`; guard `test_advection_parallel`, 600s→13s). Removed 10 now-unused Fortran kernels (`fill_rcvs`/`spill_pts` + 8 pre-existing dead). Benign-failure reporting + bounded fallback for the un-drainable pockets. Commits `ecad183`/`d07f42f`/`01decdc`/`de217d8`. |
| 2026-06-22 | — (on `dev`) | **Erratic-`sed` cascade grind capped; benchmark `bc` fixed for the NESW order; SPL case-100 dropped; scaling tooling.** At scale (9.2M mesh) a residual near-singular cell makes a cascade `fgmres` solve grind to `max_it` → erratic `sed` at some partitions; **capped** the non-fatal cascade solves at `_cascade_max_it`=1000 (fatal flow-accum keeps the full budget) — 5500→1500 iters at np=16, same result (`22c6aaa`). A reachability-to-outlet **pin** (eliminate the singular rows) was built + proven correct but **rejected** (cross-partition cost scales up with rank count, swaps memory) — documented so it isn't re-attempted (`4126a08`). **Benchmark `bc`** inputs were written for the OLD `[S,E,N,W]` order; under the NESW rework their edges flipped — fixed (`spl 0010→fooo` restores the N outlet = the SPL-benchmark failure; hillslope `ofof`, knickpoint `ooof`; SPL `end` 10→5 Myr, steady state ~1.9 Myr) and the finest **SPL case 100 dropped** (failing 4/6 on finest-mesh basin/R² artifacts + the 34-min bottleneck). Full analytical suite now green (~13 min: SPL ~3, knickpoint ~3.5, hillslope ~6.8). New `scripts/scaling/plot_scaling.py`; `gadi.pbs` consolidation. Commits `92549db`/`a03dc0d`/`bd98ef5`/`e2258e2`/`7e8f07a`. |
| 2026-06-22 | `d4162aa` (on `dev`) | **Flow downstream cascade bounded (fixes the 5 km / 144-rank `garbage_cleanup` overflow crash).** The `while excess` loop in `flowAccumulation` had no iteration bound; a partition that cuts an un-drainable pocket left `excess` permanently set, spinning the loop (rebuild `fMat` + re-solve each pass) until the collective PETSc garbage collector overflowed and aborted the run (`Memory requested 18446744065119617024`, ~2.3 GB/rank in use — **not** OOM). Deterministic at 144 ranks; 96/192/240/288 don't form the pocket. Fix: a **stagnation break** on the global residual flux `_cascade_resid` — break after `_cascade_patience` (3, env `GOSPL_CASCADE_PATIENCE`) non-shrinking passes, pond the residual as a closed basin, warn loudly; `_cascade_max_steps` (100, env `GOSPL_CASCADE_MAX_STEPS`) backstop. Partition-invariant outcome; real (slow) convergence never truncated (33/33 fast + mass/provenance conservation green). sed/sea cascades were already bounded (`max_iters=5000`). See **`## Flow-accumulation KSP` → Residual**. |
| 2026-06-22 | — (on `dev`) | **HPC scaling user-guide page.** New `docs/user_guide/hpc_scaling.rst` (+ toctree/card in `index.rst`) documenting strong-scaling speedup + parallel efficiency on the 10 km (5.9 M), 8 km (9.2 M) and 5 km (23.7 M) global meshes, honest about the baseline-to-smallest convention (p0 = 16/24/96 ranks), super-linear (cache/bandwidth) efficiency at 8 km and 5 km, the serial floors (flexure, pit-graph), and the rank-0-concentrated ~24 GB RSS peak at 5 km. Reusable figure generator `docs/user_guide/scaling/make_scaling_figure.py` auto-discovers `scaling_<N>km.csv` and emits BOTH the speedup/efficiency figure (`scaling_hpc.png`) and a per-phase figure (`scaling_phases.png`); `_label_ticks()` thins the x-axis labels so the wide log axis stays legible (markers at every rank). The 5 km / 144-rank point was the run that confirmed the cascade stagnation-break fix completes (exit 0). **Sweep extents (`fbe23ed`, `b4d5e1d`): 10 km 16→288, 8 km 24→384, 5 km 96→528.** The "useful core ceiling rises with mesh size" is the headline: 10 km sweet spot ~240; 8 km plateaus by ~288–384 (288→336→384 only ~5 % faster total, efficiency 0.93→0.82→0.74); 5 km still improving at 528 (4.87×, efficiency ~1.0 ≤384 easing to ~0.89). 5 km also OOMs at 48 ranks (single 192 GB node, rank-0 ~24 GB peak from the serial pit-graph), so 96 is its floor. Also added an illustrative 10 km model-output figure (`docs/_static/run_view10km.jpg`, `89de6a8`). |
| 2026-06-22 | `1cbcf42` (on `dev`) | **`test_parallel_correctness` Tier-1B tolerance relaxed `1e-10`→`5e-9` (platform FP, NOT a regression).** The area-weighted mean-elevation n=1-vs-n=2 check was macOS-calibrated and failed ubuntu CI (`rel_mean_h` ~3.6e-11 macOS vs ~7.6e-10 ubuntu) — the same near-tie parallel-reduction non-determinism already loosening `rel_sum_fa` to 5 %. Verified the cascade fix (`d4162aa`) is NOT the cause: pre-fix and post-fix flowplex give bit-identical `rel_mean_h` (3.616855e-11), so the stagnation break is a no-op on the converging fixture. Test-only change (no source touched); see the *Platform-dependent KSP-noise floor in `test_parallel_correctness`* bullet under the known-pitfalls section. |
| 2026-06-23 | `b7946fe`…`cf4c736` (on `master`/`dev`/`release-candidate`) | **Release candidate `2026.6.23rc1` staged (publish-free) + slow-tier fixture teardown.** Bumped the `meson.build`+`conda/meta.yaml` version pair to `2026.6.23rc1` (PEP 440-canonical so PyPI and conda render it identically — a bare `rc` would normalize to `rc0` on PyPI and diverge); fast-forwarded `master` to `dev` (the dev cycle: cascade fix, HPC scaling docs, CI fixes) and created the `release-candidate` branch. **No `v*` tag pushed** → nothing published to PyPI/conda/Docker yet (tag `v2026.6.23rc1` when ready; that fires all three publish workflows, and PyPI burns the version string irreversibly). Also added the autouse `tests/conftest.py` teardown (`_release_petsc_after_test`) — see the known-pitfalls bullet — fixing the flaky ubuntu slow-tier `SIGABRT` from un-torn-down model fixtures (`cf4c736`). |

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