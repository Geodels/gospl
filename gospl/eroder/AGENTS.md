# River incision / SPL eroders (gospl/eroder/)

Read before editing `SPL.py`, `nlSPL.py` or `soilSPL.py`. All three must use the `i`-suffix drainage arrays (`rcvIDi`, `wghtVali`, `fMati`, ...) and keep `Eb` in the thickness-rate convention (root `AGENTS.md`).

This file is part of the goSPL agent guide. Read the root `AGENTS.md` first; it holds
the rules that apply everywhere (MPI contract, solver lifecycles, scratch Vecs,
conventions, commit checklist). The text below was moved verbatim from the root
file in the 2026-10 split, so dates and `## ...` cross-references still refer to
the original section names (most now live in the files listed in the root
`AGENTS.md` > Map).

## Cached solvers in this package (full rows)

| File | Method | Cached attribute | Solver |
|---|---|---|---|
| `eroder/nlSPL.py` | `_solveNL_ed` | `self._snes_ed` (+ `self._snes_ed_fb`) | SNES, transport-limited (**`qn` default**, ngmres fallback) |
| `eroder/nlSPL.py` | `_solveNL` | `self._snes_nl` | SNES, detachment-limited (nrichardson + analytic Jacobian) |
| `eroder/soilSPL.py` | `_solveSoil` | `self._snes_soil` (+ `self._snes_soil_fb`) | SNES, soil-aware (**`qn` default**, ngmres fallback) |
| `eroder/soilSPL.py` | `diffuseSoil` | `self._ts_soil` | TS (rosw soil diffusion). Stage solve **gmres + gasm, rtol 1e-8** (was `preonly` + gasm until 2026-11: inexact and partition-dependent, 0.2 m np=1 vs np=3 on a small valley); error tolerance **atol 1e-3 m, rtol 0** (rtol 1e-3 on the absolute surface allowed ~0.2 m error per 200 m of elevation). Guard: `test_soil_diffusion_partition_invariant`. |

## Selectable solver + complementary fallback

**Selectable solver + complementary fallback** (soilSPL `_solveSoil`, nlSPL `_solveNL_ed`): built by a `_build_*_snes(primary=)` helper. The primary defaults to `qn` (L-BFGS + critical-point line search, set via the YAML `solver:` key). A *bare* `ngmres` ignores its configured KSP/PC and **stalls or diverges** on these stiff residuals (this was the historical default; it cost soil ~2.4× and nlSPL-transport ~10×, and diverged outright on the since-removed SIA ice solver; only the diagnostic glacial model remains, see `gospl/flow/AGENTS.md` > Ice sheet). On non-convergence the solve retries from the same initial guess with the **complementary** solver (`qn` ⇄ `ngmres`+`nrichardson`+hypre). **soilSPL only — adaptive sub-stepping (`_adaptiveSubstepSoil`):** if BOTH solvers still diverge, the fluvial solve is retried in `N = 4→8→16` increments of `Δt/N` (the demanded per-step incision `Kbr·Sⁿ ∝ Δt·Aᵐ` is what makes the stiff residual overshoot; `Δt/N` cuts it by `N` until it converges) and the increments sum to the full step's erosion — so erosion is *retained*, not discarded. Divergence is triggered by a **de-armoured node that just captured a large discharge** (`surfK` jumps ~10× as the duricrust/top-stratum erodes off, `PAᵐ` jumps as drainage reroutes) demanding a ~tens-of-km single-step incision; the old "continue with the best (diverged) iterate" then spiked elevation to ~−10⁶ m and cascaded (and its phantom sediment over-deposited downstream). Only a step *still* diverging after `N=16` is reverted to its prior elevation (no fluvial erosion that step) as a last resort. Helpers `_soilErodibility` (rebuild `Kbr`/`K_soil` per sub-step) and `_soilSNESsolve` (primary+fallback); the full-Δt converging path is unchanged. Each SNES carries its own options prefix (`soilspl_`/`soilsplfb_`/`nlspled_`/`nlspledfb_`) so line-search options don't leak between solvers. The `_fb` fallback SNES **and** its `_f` residual vec are cached and are in the `destroy_DMPlex` list. YAML knobs (`soil:` / `spl:` blocks): `maxIter`, `rtol`, `atol`, `pcType`, `solver`. The two **coupled two-field (h/q) deposition** solves — **linear SPL deposition** (`SPL._coupledEDSystem`, `spl_ed_` prefix) and **marine deposition** (`seaplex._depMarineSystem`, `marine_dep_` prefix) — both use a **Schur-complement** fieldsplit (defaults injected into the options DB guarded by `hasName` so `PETSC_OPTIONS`/`-spl_ed_*`/`-marine_dep_*` override). The two blocks are near-triangular M-matrices their per-block ILU solves almost exactly, so the Schur factorisation captures the full coupling and converges in ~2 outer iterations vs the old additive split (SPL ~60×/122→2 fewer Krylov iters; marine 31→2 on the minimal mesh, more at scale), solution unchanged. (Both remain AD-HOC per-call create/destroy — the Schur change is only the preconditioner, not the lifecycle.) **Block solvers must be set through the options database** (`spl_ed_fieldsplit_{h,q}_*`, `marine_dep_fieldsplit_{h,q}_*`, hasName-guarded): `subksps[i].setType(...)` before `setFromOptions` has no effect once the split becomes Schur, and until 2026-11-12 both blocks silently ran GMRES (SPL solve ~5,000 s/step at 5 km). Current: `preonly` everywhere; SPL q exact per-rank LU; marine q block-Jacobi ILU (exact LU and GASM are slow on 12-direction marine routing). Check with `-spl_ed_ksp_view`. Guard: `test_coupled_deposition_block_solves_are_direct`.

## SNES stopping criteria (2026-11-24)

**The soil SNES primary is `ngmres` + hypre by default (2026-11-25).** `qn` (L-BFGS, `cp` line search) was the default; with `stol=0` it failed at nearly every step on the 47k/190k regional margin and on the global continental_geochem example (reason -5 at 500 its, or -9), and the `ngmres` fallback finished at `rtol*100` = 1e-4, the tolerance the user guide warns gives needles. `ngmres` primary: no fallback, full tolerance, same or lower wall time (109 % vs 187 % of the base step at 47k). `solver: qn` remains available.

Every outer SNES (`soilspl_`, `soilsplfb_`, `nlspled_`, `nlspledfb_`, the no-deposition `_snes_nl`, `hillnl_`) sets `stol=0.0`: PETSc's default `stol = 1e-8` let the quasi-Newton solves report `CONVERGED_SNORM_RELATIVE` when the line search stalled, with the residual ~1e3 times `rtol`. Check with `-<prefix>snes_converged_reason`; a solve must end on `CONVERGED_FNORM_ABS`/`FNORM_RELATIVE` (or diverge to the fallbacks). Guard `test_snes_converge_on_the_residual`. Any new SNES: set `stol=0.0` too.

## The `fDep` cap (transport-limited deposition)

See `gospl/tools/AGENTS.md` > Debugging an elevation spike for the full text. Rule of
thumb: keep `spl: G ≲ rain/2` (rain in m/yr); at the 0.99 cap the coupled `(1−fDep)`
block is near-singular and trunk-river loads land on a handful of cells.

## Detachment-limited solve in increment form (`SPL._getEroDepRate`, 2026-11-24)

With `fDepa == 0` the implicit SPL system `E h = hOld` is solved for the increment `dh = h - hOld` (`E dh = hOld - E hOld`, zero initial guess), then `h = hOld + dh`. Solving for the absolute elevation let the KSP `rtol` of 1e-8 act on elevations of thousands of metres: an error of ~1e-4 m that changed with the block-Jacobi partition. The increment form puts the tolerance on the erosion itself. Guard: `test_base_model_partition_exact`.
