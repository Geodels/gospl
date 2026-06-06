"""
Hillslope diffusion benchmark for goSPL.

Adapted from benchmarks/hillslope/test-hillslope.py. The orchestration —
configure/tolerances, the 7-model run, the 7 extractOutputs calls, and
runFullBenchmark — is preserved verbatim from the original; only the
entry point is wrapped in a pytest function and each Model() lifecycle
is enforced with try/finally.

Imports `scripts.analysis as anlys` from the COPIED hillslope folder
(prepended to sys.path by hillslope_tmp_path).

Analytical basis: z(x) = (U / 2κ) · x · (L − x)
Reference: Roering, Kirchner & Dietrich (1999); Perron (2017).
See AGENTS.md > Analytical benchmark suite.
"""

import os
import time

import pytest

# Optional dependency guards — silent skip in environments without these.
scipy = pytest.importorskip("scipy")
matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")   # non-interactive backend

from petsc4py import PETSc
PETSc.Options().setValue("-options_left", "0")
from gospl.model import Model as sim


# ══════════════════════════════════════════════════════════════════════════════
# Model parameters — edit these to match your goSPL run
# ══════════════════════════════════════════════════════════════════════════════
U     = 5e-4   # m/yr  uplift rate
KAPPA = 0.1    # m²/yr hillslope diffusivity
L     = 1000.0 # m     domain length (x direction)

# ══════════════════════════════════════════════════════════════════════════════
# Test tolerances
# ══════════════════════════════════════════════════════════════════════════════
TOL_TAU_PCT      = 10.0   # % — fitted τ vs theoretical τ
TOL_PEAK_PCT     = 5.0    # % — asymptotic z_ss vs analytical z_max
TOL_PROFILE_R2   = 0.990  # R² of spatial profile vs parabola
TOL_CONVERGENCE  = 99.0   # % — model must reach this % of its own asymptote
TOL_MONO         = True   # z_max must increase monotonically
TOL_EARLY_PCT    = 5.0    # % — early-time linear growth rate vs analytical


@pytest.mark.benchmark
@pytest.mark.slow
def test_hillslope_diffusion(hillslope_tmp_path):
    """
    Hillslope diffusion benchmark.

    Validates: diffusion solver, steady-state convergence, parabolic
    profile shape.
    Analytical basis: z(x) = (U/2κ)·x·(L-x)
    Pass criteria: all tolerances in TOL_* constants must be met.
    See AGENTS.md: Analytical benchmark suite.
    """
    # Late import — only reachable after hillslope_tmp_path has prepended
    # the copied benchmark folder onto sys.path. Module-level import would
    # resolve to the SOURCE-tree scripts/analysis.py and fight chdir.
    from scripts import analysis as anlys

    anlys.configure(U, KAPPA, L)
    anlys.tolerances(TOL_TAU_PCT, TOL_PEAK_PCT, TOL_PROFILE_R2,
                     TOL_CONVERGENCE, TOL_MONO, TOL_EARLY_PCT)

    # ------------------------------------------------------------------
    # 1. Run 7 goSPL models — each wrapped in try/finally per AGENTS.md
    #    KSP lifecycle contract.
    # ------------------------------------------------------------------
    outfolder = "sims"
    yml_files = [
        "hillslope_base.yml",
        "hillslope_dx25.yml", "hillslope_dx10.yml",
        "hillslope_dt_5000.yml", "hillslope_dt_500.yml",
        "hillslope_kappa2.yml", "hillslope_U2.yml",
    ]

    os.makedirs("sims_outputs", exist_ok=True)

    time_start = time.time()
    for fname in yml_files:
        model = sim(f"{outfolder}/{fname}", False, False)
        try:
            model.runProcesses()
        finally:
            model.destroy()
    print(f"Total run time : {time.time() - time_start:.1f} s")

    # ------------------------------------------------------------------
    # 2. Extract outputs for each run.
    #    NOTE: configure() must be called BEFORE each extractOutputs so
    #    the analysis module's _set dict carries the right U/KAPPA for
    #    that run. This is exactly the order the original script uses.
    # ------------------------------------------------------------------
    anlys.configure(U=5e-4, KAPPA=0.1, L=1000.0)
    r_base   = anlys.extractOutputs(80, tol=20, model_path="sims_outputs/hillslope_base",    plotdir="results")
    r_dx25   = anlys.extractOutputs(80, tol=20, model_path="sims_outputs/hillslope_dx25",    plotdir="results")
    r_dx10   = anlys.extractOutputs(80, tol=20, model_path="sims_outputs/hillslope_dx10",    plotdir="results")
    r_dt5000 = anlys.extractOutputs(80, tol=20, model_path="sims_outputs/hillslope_dt_5000", plotdir="results")
    r_dt500  = anlys.extractOutputs(80, tol=20, model_path="sims_outputs/hillslope_dt_500",  plotdir="results")

    anlys.configure(U=1e-3, KAPPA=0.1, L=1000.0)
    r_2U     = anlys.extractOutputs(80, tol=20, model_path="sims_outputs/hillslope_U2",      plotdir="results")

    anlys.configure(U=5e-4, KAPPA=0.2, L=1000.0)
    r_kappa2 = anlys.extractOutputs(80, tol=20, model_path="sims_outputs/hillslope_kappa2",  plotdir="results")

    # ------------------------------------------------------------------
    # 3. Combined summary across all groups (G1..G5).
    #    runFullBenchmark internally calls write_ci_report and
    #    write_ci_annotations — both go to results/.
    # ------------------------------------------------------------------
    anlys.configure(U=5e-4, KAPPA=0.1, L=1000.0)   # restore base
    results = anlys.runFullBenchmark(
        results_base   = r_base,
        results_2U     = r_2U,
        results_kappa2 = r_kappa2,
        results_dx25   = r_dx25,
        results_dx10   = r_dx10,
        results_dt5000 = r_dt5000,
        results_dt500  = r_dt500,
        plotdir        = "results",
    )

    # ------------------------------------------------------------------
    # 4. Pipe the markdown report into $GITHUB_STEP_SUMMARY if running
    #    under GitHub Actions.
    # ------------------------------------------------------------------
    md_path = hillslope_tmp_path / "results" / "benchmark_report.md"
    github_summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if github_summary and md_path.exists():
        with open(md_path) as f:
            content = f.read()
        with open(github_summary, "a") as f:
            f.write(content)

    # ------------------------------------------------------------------
    # 5. Assert pass/fail.
    #    anlys.runFullBenchmark does NOT expose an overall_pass key on
    #    its return dict, but anlys.print_summary(results) returns the
    #    bool we need (True iff every boolean flag in results is True).
    # ------------------------------------------------------------------
    overall_pass = anlys.print_summary(results)
    assert overall_pass, (
        "Hillslope benchmark failed — see "
        f"{md_path} and {hillslope_tmp_path / 'results'}/ for details"
    )
