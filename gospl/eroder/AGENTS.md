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
| `eroder/soilSPL.py` | `diffuseSoil` | `self._ts_soil` | TS (rosw soil diffusion) |

## Selectable solver + complementary fallback

**Selectable solver + complementary fallback** (soilSPL `_solveSoil`, nlSPL `_solveNL_ed`): built by a `_build_*_snes(primary=)` helper. The primary defaults to `qn` (L-BFGS + critical-point line search, set via the YAML `solver:` key). A *bare* `ngmres` ignores its configured KSP/PC and **stalls or diverges** on these stiff residuals (this was the historical default; it cost soil ~2.4× and nlSPL-transport ~10×, and diverged outright on the SIA — which is why SIA is now solved explicitly, see "Ice sheet"). On non-convergence the solve retries from the same initial guess with the **complementary** solver (`qn` ⇄ `ngmres`+`nrichardson`+hypre). **soilSPL only — adaptive sub-stepping (`_adaptiveSubstepSoil`):** if BOTH solvers still diverge, the fluvial solve is retried in `N = 4→8→16` increments of `Δt/N` (the demanded per-step incision `Kbr·Sⁿ ∝ Δt·Aᵐ` is what makes the stiff residual overshoot; `Δt/N` cuts it by `N` until it converges) and the increments sum to the full step's erosion — so erosion is *retained*, not discarded. Divergence is triggered by a **de-armoured node that just captured a large discharge** (`surfK` jumps ~10× as the duricrust/top-stratum erodes off, `PAᵐ` jumps as drainage reroutes) demanding a ~tens-of-km single-step incision; the old "continue with the best (diverged) iterate" then spiked elevation to ~−10⁶ m and cascaded (and its phantom sediment over-deposited downstream). Only a step *still* diverging after `N=16` is reverted to its prior elevation (no fluvial erosion that step) as a last resort. Helpers `_soilErodibility` (rebuild `Kbr`/`K_soil` per sub-step) and `_soilSNESsolve` (primary+fallback); the full-Δt converging path is unchanged. Each SNES carries its own options prefix (`soilspl_`/`soilsplfb_`/`nlspled_`/`nlspledfb_`) so line-search options don't leak between solvers. The `_fb` fallback SNES **and** its `_f` residual vec are cached and are in the `destroy_DMPlex` list. YAML knobs (`soil:` / `spl:` blocks): `maxIter`, `rtol`, `atol`, `pcType`, `solver`. The two **coupled two-field (h/q) deposition** solves — **linear SPL deposition** (`SPL._coupledEDSystem`, `spl_ed_` prefix) and **marine deposition** (`seaplex._depMarineSystem`, `marine_dep_` prefix) — both use a **Schur-complement** fieldsplit (defaults injected into the options DB guarded by `hasName` so `PETSC_OPTIONS`/`-spl_ed_*`/`-marine_dep_*` override). The two blocks are near-triangular M-matrices their per-block ILU solves almost exactly, so the Schur factorisation captures the full coupling and converges in ~2 outer iterations vs the old additive split (SPL ~60×/122→2 fewer Krylov iters; marine 31→2 on the minimal mesh, more at scale), solution unchanged. (Both remain AD-HOC per-call create/destroy — the Schur change is only the preconditioner, not the lifecycle.)

## The `fDep` cap (transport-limited deposition)

See `gospl/tools/AGENTS.md` > Debugging an elevation spike for the full text. Rule of
thumb: keep `spl: G ≲ rain/2` (rain in m/yr); at the 0.99 cap the coupled `(1−fDep)`
block is near-singular and trunk-river loads land on a handful of cells.
