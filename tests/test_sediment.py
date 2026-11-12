"""
Sediment routing and marine deposition: mass conservation, marine diffusion, oFill.

Protects: AGENTS.md > Fixed (marine leak, steep-bathymetry clamp, oFill).

Split out of the former tests/test_regression.py. Run this group alone with
`pytest -m sediment`; see tests/README.md for the marker list.
"""

from __future__ import annotations

import numpy as np
import pytest

# Skip the whole module (rather than erroring at collection) when the
# goSPL runtime stack is not installed.
pytest.importorskip(
    "gospl.tools.inputparser",
    reason="goSPL runtime deps (petsc4py / ruamel.yaml / scipy) not installed",
)

pytestmark = [pytest.mark.sediment]


# ---------------------------------------------------------------------------
# TEST 6 - True mass conservation on a closed sphere
# ---------------------------------------------------------------------------


# @pytest.mark.xfail(
#     strict=True,
#     reason=(
#         "Discovered 2026-06 on a global sphere fixture: ~25% of "
#         "redistributed sediment vanishes per run, with dV_surface == "
#         "dV_cumED (h<->cumED sync is correct). Suspects: (1) seaplex."
#         "_distOcean discards residual `sinkVol` when the convergence "
#         "loop exits (seaplex.py:340-417); (2) seaplex._matOcean builds "
#         "a dMat1 with zero columns at true sinks, so sediment that "
#         "lands at a saturated ocean-basin floor is annihilated by the "
#         "next dMat1.mult. Fix in a separate sprint; strict=True means "
#         "the suite WILL fail once the leak is closed, forcing this "
#         "marker to be removed and turning the test into a permanent "
#         "regression guard. See diagnostic numbers in the assertion "
#         "messages below for the original baseline."
#     ),
# )
@pytest.mark.slow
def test_mass_conservation(minimal_model):
    """
    Protects: TRUE sediment conservation on a closed domain.

    Silent failure prevented: any kernel that creates or destroys
    sediment mass — or scales an erosion/deposition rate by the wrong
    factor — will drift the global cumED integral away from zero on a
    closed sphere. The previous skeleton only verified that hGlobal
    and cumED moved in lockstep; it would pass green if BOTH were
    wrong in the same way (e.g. a doubled scaling applied to both
    `hGlobal.axpy` and `cumED.axpy` calls).

    Strong invariants on a closed sphere with no tectonics, compaction,
    flexure, or paleo-Z resetting:

        |dV_cumED|   / total_activity < TOLERANCE
        |dV_surface| / total_activity < TOLERANCE

    where `total_activity = sum(|cumED_change| * larea)` is the volume
    of sediment redistributed during the run. Using activity (not 1.0)
    as the scale makes the bound meaningful: we want 0.01 % imbalance
    against the mass that actually moved, not against 1 m^3.

    Non-applicability gate: this test SKIPS cleanly when the fixture
    has any kernel that writes hGlobal WITHOUT a paired cumED update
    (tectonic uplift via `upsub`, compaction, flexure, paleo-Z reset).
    Those processes are not bugs; they are real h-modifiers outside
    the sediment budget. To exercise this test, the fixture must be a
    closed sphere with only sediment-conserving kernels active.
    """
    model = minimal_model

    # ---- Gate: this strong test only applies on a closed domain ----
    reasons = []
    if getattr(model, "flatModel", True):
        reasons.append("flatModel=True (2D plane, has boundary outflux)")
    if getattr(model, "tecdata", None) is not None:
        reasons.append("tectonics is active (upsub adds/removes mass without cumED)")
    if getattr(model, "flexOn", False):
        reasons.append("flexure is active (hGlobal moves without cumED)")
    if getattr(model, "stratNb", 0) > 0:
        reasons.append("stratigraphy is active (compaction shrinks h without cumED)")
    if getattr(model, "paleoZ", None) is not None:
        reasons.append("paleoZ reset is active (overwrites h below sea level)")
    if reasons:
        pytest.skip(
            "Strong mass conservation requires a closed sphere with no "
            "non-sediment h-writers. This fixture has: "
            + "; ".join(reasons)
            + ". NEEDS_HUMAN_REVIEW: either commit a second clean-sphere "
            "fixture (tests/fixtures/closed_sphere.yml) and bind it to "
            "this test via its own fixture, or extend the budget below "
            "to subtract dV_tectonic / dV_compaction / dV_flexure."
        )

    # ---- Run and capture state ----
    h_before = model.hLocal.getArray().copy()
    cumED_before = model.cumEDLocal.getArray().copy()
    larea = model.larea
    owned = model.inIDs == 1

    model.runProcesses()

    h_after = model.hLocal.getArray().copy()
    cumED_after = model.cumEDLocal.getArray().copy()

    dh = h_after - h_before
    dED = cumED_after - cumED_before

    # ---- Reduce over owned cells across ranks ----
    from mpi4py import MPI
    dV_surface_local = float(np.sum((dh * larea)[owned]))
    dV_cumED_local = float(np.sum((dED * larea)[owned]))
    activity_local = float(np.sum((np.abs(dED) * larea)[owned]))

    dV_surface = MPI.COMM_WORLD.allreduce(dV_surface_local, op=MPI.SUM)
    dV_cumED = MPI.COMM_WORLD.allreduce(dV_cumED_local, op=MPI.SUM)
    total_activity = MPI.COMM_WORLD.allreduce(activity_local, op=MPI.SUM)

    if total_activity < 1.0:
        pytest.skip(
            f"Total sediment activity ({total_activity:.3e} m^3) is below "
            f"1 m^3 over the run. The fixture did not move enough sediment "
            f"for a conservation check to be meaningful. "
            f"NEEDS_HUMAN_REVIEW: lengthen the run or steepen the gradient."
        )

    # ---- Tolerance pinned 2026-06 against tests/fixtures/minimal.yml ----
    # 1e-4 relative against total_activity. The underlying KSP runs at
    # rtol=1e-10 so true closure is tighter; the gap is dominated by
    # legitimate floor effects:
    #   - DEPOSIT_FLOOR=1e-3 (seaplex.py:465, sedplex.py:156) drops
    #     sub-mm deposits as numerical noise; over many cells this
    #     accumulates to a small but non-zero mass loss.
    #   - Pit-routing residue (sedplex.py:184 threshold ~1e-3 m^3) leaves
    #     a small unrouted excess inside the pit volumes.
    # Before raising this bound: demonstrate that a 1% scaling error in
    # SPL.py:352 (multiply `-Eb * self.dt` by 1.01) still trips the
    # assertion. If it doesn't, 1e-4 is already too loose — tighten,
    # don't loosen.
    TOLERANCE = 1.0e-4

    rel_cumED = abs(dV_cumED) / total_activity
    rel_surface = abs(dV_surface) / total_activity

    diagnostic = (
        f"\n  total_activity = {total_activity:.3e} m^3 "
        f"(sediment redistributed)\n"
        f"  dV_cumED       = {dV_cumED:+.3e} m^3  "
        f"(relative {rel_cumED:.3e})\n"
        f"  dV_surface     = {dV_surface:+.3e} m^3  "
        f"(relative {rel_surface:.3e})\n"
        f"  tolerance      = {TOLERANCE:.0e}"
    )

    # ---- Strong invariant A: sediment closure ----
    assert rel_cumED < TOLERANCE, (
        f"Sediment is being created or destroyed: cumED integral does "
        f"not close to zero on a closed sphere. Either a kernel is "
        f"scaling its `cumED.axpy` term incorrectly, or the SPL / "
        f"sed / hillslope kernels are not internally mass-conserving."
        + diagnostic
    )

    # ---- Strong invariant B: surface volume closure ----
    assert rel_surface < TOLERANCE, (
        f"Surface volume is drifting on a closed sphere with no "
        f"tectonics/compaction/flexure: a kernel is writing to "
        f"hGlobal without a matching cumED update, OR the sediment "
        f"kernels are not internally mass-conserving."
        + diagnostic
    )


# ---------------------------------------------------------------------------
# TEST 7 - Sediment conservation must hold under evaporation
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_sediment_balance_with_evap(minimal_model_with_evap):
    """
    Protects: DESIGN_EVAPORATION.md §2.5 — the evap feature couples to
    sediment ONLY one-way (reduced FA → reduced erosion → reduced
    sediment flux). The sediment pit-fill machinery in sedplex.py uses
    its own pitVol state read fresh at sedplex.py:648, which is
    invariant under whatever the water-side `_distributeDownstream`
    did with its local `pitVol`. This test verifies that decoupling
    still holds when evap is active.

    Silent failures prevented:
      1. A future refactor that promotes `pitVol` from a local var to
         `self.pitVol` inside flowplex.py would corrupt sediment's view
         of the pit (sedplex reads `self.pitParams[:, 0]` at line 648,
         but if the water side has already modified it, sediment over-
         or under-deposits).
      2. A future "improvement" that subtracts `self.evapLoss` from the
         sediment budget (mistakenly trying to balance water leaks
         against sediment) would create or destroy sediment mass.
      3. A future refactor that makes lake-removal write to `cumED`
         (treating "removed lake water" as if it deposited sediment)
         would create mass from nothing.

    Strategy: identical to test_mass_conservation, but with evap injected
    via the `minimal_model_with_evap` fixture (uniform evap = 0.5 × rUni).
    If sediment closure remains within 1e-4 relative, the decoupling is
    intact.
    """
    model = minimal_model_with_evap

    # ---- Same gate as test_mass_conservation: closed-sphere only ----
    reasons = []
    if getattr(model, "flatModel", True):
        reasons.append("flatModel=True (2D plane, has boundary outflux)")
    if getattr(model, "tecdata", None) is not None:
        reasons.append("tectonics is active (upsub adds/removes mass without cumED)")
    if getattr(model, "flexOn", False):
        reasons.append("flexure is active (hGlobal moves without cumED)")
    if getattr(model, "stratNb", 0) > 0:
        reasons.append("stratigraphy is active (compaction shrinks h without cumED)")
    if getattr(model, "paleoZ", None) is not None:
        reasons.append("paleoZ reset is active (overwrites h below sea level)")
    if reasons:
        pytest.skip(
            "Strong mass conservation requires a closed sphere with no "
            "non-sediment h-writers. This fixture has: "
            + "; ".join(reasons)
        )

    # ---- Run and capture state ----
    h_before = model.hLocal.getArray().copy()
    cumED_before = model.cumEDLocal.getArray().copy()
    larea = model.larea
    owned = model.inIDs == 1

    model.runProcesses()

    h_after = model.hLocal.getArray().copy()
    cumED_after = model.cumEDLocal.getArray().copy()

    dh = h_after - h_before
    dED = cumED_after - cumED_before

    from mpi4py import MPI
    dV_surface_local = float(np.sum((dh * larea)[owned]))
    dV_cumED_local = float(np.sum((dED * larea)[owned]))
    activity_local = float(np.sum((np.abs(dED) * larea)[owned]))

    dV_surface = MPI.COMM_WORLD.allreduce(dV_surface_local, op=MPI.SUM)
    dV_cumED = MPI.COMM_WORLD.allreduce(dV_cumED_local, op=MPI.SUM)
    total_activity = MPI.COMM_WORLD.allreduce(activity_local, op=MPI.SUM)

    if total_activity < 1.0:
        pytest.skip(
            f"Sediment activity {total_activity:.3e} m^3 too small for a "
            f"meaningful check. With 50% evap reducing erosion, this can "
            f"happen on fixtures where the no-evap version barely cleared "
            f"the threshold. NEEDS_HUMAN_REVIEW: lower the fixture evap "
            f"rate in conftest.py::minimal_model_with_evap, or use a "
            f"steeper-gradient fixture."
        )

    # Same tolerance as test_mass_conservation — see that test for the
    # rationale (KSP rtol=1e-10 but DEPOSIT_FLOOR and pit-routing residue
    # dominate the gap).
    TOLERANCE = 1.0e-4

    rel_cumED = abs(dV_cumED) / total_activity
    rel_surface = abs(dV_surface) / total_activity

    # Sanity: with evap enabled, evapLoss should be > 0.
    evap_total = MPI.COMM_WORLD.allreduce(float(model.evapLoss), op=MPI.SUM)

    diagnostic = (
        f"\n  evap_total     = {evap_total:.3e} m^3 (sanity > 0)"
        f"\n  total_activity = {total_activity:.3e} m^3 "
        f"(sediment redistributed)"
        f"\n  dV_cumED       = {dV_cumED:+.3e} m^3  "
        f"(relative {rel_cumED:.3e})"
        f"\n  dV_surface     = {dV_surface:+.3e} m^3  "
        f"(relative {rel_surface:.3e})"
        f"\n  tolerance      = {TOLERANCE:.0e}"
    )

    # Sanity check: the fixture is supposed to have evap on. If evapLoss
    # is zero, the fixture is broken and the test below is meaningless.
    assert evap_total > 0.0, (
        "Fixture minimal_model_with_evap did not accumulate any evap "
        "this run. The fixture or applyForces wiring may be broken; "
        "test_sediment_balance_with_evap cannot prove the decoupling "
        "if evap never fires." + diagnostic
    )

    assert rel_cumED < TOLERANCE, (
        "Sediment mass is being created or destroyed when evap is "
        "active. The evap feature has leaked into the sediment side. "
        "Check that:\n"
        "  (a) flowplex.py:_distributeDownstream still modifies the "
        "LOCAL `pitVol` var (not `self.pitParams[:, 0]`).\n"
        "  (b) sedplex.py:648 still re-initialises sediment's pitVol "
        "from `self.pitParams[:, 0]` (water-side hasn't touched it).\n"
        "  (c) No new caller has added `self.evapLoss` into a cumED "
        "or hGlobal axpy."
        + diagnostic
    )

    assert rel_surface < TOLERANCE, (
        "Surface volume drifting on closed sphere with evap. Same "
        "diagnosis as the cumED case." + diagnostic
    )


# ---------------------------------------------------------------------------
# TEST 8d - Cached marine TS: per-call step-counter reset
# ---------------------------------------------------------------------------
#
# hillslope._diffuseImplicit (marine diffusion) and soilSPL.diffuseSoil reuse a
# cached PETSc TS. `ts.setTime(0.0)` resets the clock each call but NOT the step
# counter (`getStepNumber()`), so without an explicit `ts.setStepNumber(0)` the
# counter accumulates across calls. `ts.setMaxSteps(self.tsStep)` then acts as a
# CUMULATIVE cap: after ~tsStep total substeps (a few hundred model steps at the
# default tsStep=2000) the cap is already exceeded on entry, TSSolve returns
# immediately, and the marine deposit is left un-diffused — silently. (It also
# made the verbose substep/iteration counts grow without bound.)
#
# This guard runs the marine fixture for several steps and asserts the per-call
# TS step count stays bounded (independent per call) instead of accumulating.
# ---------------------------------------------------------------------------


def test_marine_ts_step_counter_resets(minimal_model):
    """
    Protects: hillslope._diffuseImplicit must `ts.setStepNumber(0)` each call so
    the cached TS's `setMaxSteps(tsStep)` is a PER-CALL budget, not a cumulative
    cap that eventually stops marine diffusion silently. Without the reset,
    `getStepNumber()` grows monotonically across the reused TS (7, 14, 21, ...).
    """
    m = minimal_model
    if getattr(m, "flatModel", False):
        pytest.skip("needs a marine domain (global sphere) to run _diffuseImplicit")

    orig = m._diffuseImplicit
    counts = []

    def wrap(*args, **kwargs):
        out = orig(*args, **kwargs)
        if getattr(m, "_ts_marine", None) is not None:
            counts.append(m._ts_marine.getStepNumber())
        return out

    m._diffuseImplicit = wrap
    try:
        m.runProcesses()
    finally:
        m.destroy()

    if len(counts) < 3:
        pytest.skip(
            f"marine diffusion ran only {len(counts)} time(s); too few to "
            f"distinguish per-call from cumulative step counts."
        )

    # Per-call (fixed): counts stay ~flat -> max ~ min. Cumulative (bug): counts
    # grow ~linearly with the call index -> max == N*min. The 2x bound clears the
    # mild physical drift in substep count while failing hard on accumulation.
    assert max(counts) <= 2 * min(counts), (
        "Cached marine TS step counter accumulates across calls "
        f"(getStepNumber per call = {counts}). hillslope._diffuseImplicit must "
        "call ts.setStepNumber(0) each invocation; otherwise setMaxSteps(tsStep) "
        "is a cumulative cap that silently stops marine diffusion on long runs."
    )


# ---------------------------------------------------------------------------
# TEST 8e - Opt-in lagged-diffusivity (Picard) marine solver
# ---------------------------------------------------------------------------
#
# `diffusion: marineSolver: picard` selects hillslope._diffuseImplicitPicard
# (lagged-diffusivity backward-Euler with linear solves) instead of the default
# adaptive non-linear TS. It is an opt-in approximation: on the minimal fixture
# it matches the TS deposit to ~1e-5; on large stiff marine inputs it is much
# faster (no kink rejections) at a small deposit difference. This test checks
# (a) the opt-in is parsed, (b) the Picard run conserves mass, and (c) its
# deposit matches the default TS on the minimal fixture (closed sphere).
# ---------------------------------------------------------------------------


def test_marine_picard_solver(minimal_model, minimal_picard_model):
    """
    Protects: the opt-in `diffusion: marineSolver: picard` path
    (inputparser._extraHillslope -> hillslope._diffuseImplicitPicard). Asserts
    the flag is parsed, the Picard marine/lake diffusion conserves mass, and its
    cumulative erosion/deposition matches the default TS solver on the minimal
    closed-sphere fixture (the approximation is exact there).
    """
    mp = minimal_picard_model
    assert getattr(mp, "marineSolver", "ts") == "picard", (
        "diffusion.marineSolver: picard was not parsed into self.marineSolver"
    )

    def _run(m):
        try:
            m.runProcesses()
            owned = m.inIDs == 1
            ed = m.cumEDLocal.getArray()[owned]
            la = m.larea[owned]
            ednorm = float(np.sum(ed ** 2)) ** 0.5
            vol = float(np.sum(ed * la))
            act = float(np.sum(np.abs(ed) * la))
            return ednorm, abs(vol) / max(act, 1.0e-30)
        finally:
            m.destroy()

    if getattr(minimal_model, "flatModel", False):
        pytest.skip("needs a closed sphere (flatModel=False) for mass conservation")

    ed_pic, mass_pic = _run(mp)
    ed_ts, mass_ts = _run(minimal_model)

    # Mass conservation on the closed sphere (same gate as test_mass_conservation).
    assert mass_pic < 1.0e-4, (
        f"Picard marine solver breaks mass conservation (rel {mass_pic:.2e})."
    )
    # Deposit must match the default TS solver to a tight tolerance on the
    # minimal fixture (the lagged-diffusivity approximation is exact here).
    rel = abs(ed_pic - ed_ts) / max(ed_ts, 1.0e-30)
    assert rel < 1.0e-2, (
        f"Picard deposit diverges from the TS solver on the minimal fixture "
        f"(cumED rel diff {rel:.2e} > 1e-2): ed_pic={ed_pic:.6e} "
        f"ed_ts={ed_ts:.6e}. The opt-in approximation should be near-exact here."
    )


def test_marine_diffusion_conserves_on_steep_bathymetry(minimal_model):
    """
    Protects: mass conservation of the marine sediment diffusion
    (`seaplex.seaChange` -> `hillslope._diffuseOcean`) on steep seafloor relief.

    Silent failure prevented: `_diffuseImplicit` / `_diffuseImplicitPicard`
    solve the marine diffusion on the ABSOLUTE seafloor (`bed + deposit`) and
    then clamp the result to >= 0 (a fresh deposit cannot erode the pre-existing
    bed). Over steep bathymetry that clamp is asymmetric -- it keeps the
    basin-fill half of the otherwise volume-conserving smoothing and discards
    the paired highs-erosion half -- so the diffused deposit can carry several
    times its input volume, i.e. marine sediment appears from nowhere and breaks
    the global mass budget (observed ~5-8x on a real-Earth run). `_diffuseOcean`
    rescales the (>=0) diffused deposit back to the incoming marine volume,
    mirroring the per-pit rescale the continental large-pit path already uses
    (`sedplex._diffuseLargePit`; dep.rst "A final per-pit mass rescale absorbs
    any boundary drift introduced by the diffusion").

    Why `test_mass_conservation` did not catch this: its backing fixture
    (minimal.yml) uses a tiny `nlK`/`dt` so its marine deposit barely diffuses,
    and it deposits little marine sediment overall, so the clamped volume stays
    under the 1e-4 budget tolerance. This test drives `_diffuseOcean` directly
    with a THIN deposit on the mesh's natural STEEP bathymetry and a strong
    diffusivity, where the pre-fix clamp inflates the volume by ~2x.
    """
    from mpi4py import MPI

    model = minimal_model
    if getattr(model, "flatModel", False):
        pytest.skip("needs a sphere fixture with submarine relief")

    owned = model.inIDs == 1
    larea = model.larea

    # Natural (steep) mesh bathymetry as the bed; raise sea level so a large
    # steep-relief region is submarine. A THIN deposit on that steep relief is
    # what drives the absolute-surface diffusion below the pre-deposition bed.
    bed = model.hLocal.getArray().copy()
    model.sealevel = 2000.0
    model.seaID = bed < model.sealevel

    n_sea = MPI.COMM_WORLD.allreduce(int(np.sum(model.seaID & owned)), op=MPI.SUM)
    if n_sea < 10:
        pytest.skip(
            f"fixture has too few submarine cells ({n_sea}) at sealevel=2000 "
            f"for a meaningful steep-bathymetry marine-diffusion test"
        )

    # Strong marine diffusion so the deposit spreads far over the steep bed and
    # the clamp actually fires (with the fixture's tiny nlK/dt it would not).
    model.dt = 5.0e4
    model.nlK = 5.0e4

    dh = np.zeros(model.lpoints)
    dh[model.seaID] = 5.0                      # thin deposit (m)
    vin = MPI.COMM_WORLD.allreduce(
        float(np.sum((dh * larea)[owned])), op=MPI.SUM
    )
    assert vin > 0.0

    # Run ONLY the marine diffusion; it writes the applied deposit into self.tmp.
    model._diffuseOcean(dh)
    model.dm.globalToLocal(model.tmp, model.tmpL)
    applied = model.tmpL.getArray().copy()
    vout = MPI.COMM_WORLD.allreduce(
        float(np.sum((applied * larea)[owned])), op=MPI.SUM
    )

    rel = abs(vout - vin) / vin
    # Pre-fix this is ~1.0 (deposit volume ~doubles); the rescale drives it to 0.
    assert rel < 1.0e-6, (
        f"marine diffusion does not conserve volume on steep bathymetry: "
        f"Vin={vin:.4e} Vout={vout:.4e} ratio={vout / vin:.3f} rel={rel:.3e}. "
        f"The _diffuseOcean mass-conservation rescale is missing or wrong."
    )

    # Physical: a fresh marine deposit must never erode the pre-existing bed.
    min_applied = MPI.COMM_WORLD.allreduce(
        float(np.min(applied[owned])) if owned.any() else 0.0, op=MPI.MIN
    )
    assert min_applied >= -1.0e-9, (
        f"marine diffusion eroded the pre-existing bed (min deposit "
        f"{min_applied:.3e} m < 0)."
    )


@pytest.mark.slow
def _steep_marine_setup(model, dt=5.0e4, nlK=5.0e4):
    """Raise sea level over the fixture's natural relief (steep bathymetry)."""
    bed = model.hLocal.getArray().copy()
    model.sealevel = 2000.0
    model.seaID = bed < model.sealevel
    model.dt = dt
    model.nlK = nlK
    return bed


@pytest.mark.parametrize("solver", ["ts", "picard"])
def test_marine_diffusion_bare_bed_does_not_shed(minimal_model, solver):
    """
    Protects: the upslope (donor) face diffusivity of the marine/lake deposit
    diffusion (`hillslope._donorDiffCoeffs`, used by the TS residual/Jacobian
    and the Picard operator).

    Silent failure prevented: with the averaged face diffusivity
    `0.5*(Cd_k + Cd_n)` a bare cell above a deposit lost BED material into it;
    the >= 0 clip dropped that loss and kept the gain, so a 5 m deposit on
    steep bathymetry came out at 2.17x its volume with a 456 m peak, under BOTH
    solvers. `_diffuseOcean`'s mass rescale then scaled every deposit by ~0.5
    (seen as `Marine diffusion mass rescale: 0.51` on a 10 km global run).
    The raw solver output must conserve volume BEFORE any rescale.
    """
    from mpi4py import MPI

    model = minimal_model
    if getattr(model, "flatModel", False):
        pytest.skip("needs a sphere fixture with submarine relief")
    try:
        _steep_marine_setup(model)
        owned = model.inIDs == 1
        dh = np.zeros(model.lpoints)
        dh[model.seaID] = 5.0
        model.marineSolver = solver
        model.picardSub, model.picardIts = 40, 4
        vin = MPI.COMM_WORLD.allreduce(float(np.sum((dh * model.larea)[owned])))
        nd = model._diffuseImplicit(dh, model.seaID, model.nlK, label="marine")
        vout = MPI.COMM_WORLD.allreduce(float(np.sum((nd * model.larea)[owned])))
        peak = MPI.COMM_WORLD.allreduce(float(nd[owned].max()), op=MPI.MAX)
    finally:
        model.destroy()

    # Pre-fix: vout/vin = 2.17 (both solvers), peak 452-457 m.
    assert abs(vout / vin - 1.0) < 1.0e-6, (
        f"{solver}: raw marine diffusion volume ratio {vout / vin:.4f} != 1; "
        f"bare bed is shedding material into the deposit (face diffusivity)."
    )
    assert peak < 100.0, f"{solver}: 5 m deposit diffused to a {peak:.1f} m peak"


def test_marine_picard_supply_limited(minimal_model):
    """
    Protects: the supply limiter of the Picard deposit diffusion
    (`_donorDiffCoeffs(limit=(avail, dt_sub))` in `_diffuseImplicitPicard`).

    Silent failure prevented: with a large diffusion number `Cd*dt_sub/dx^2`
    (fine meshes: 10 km cells, `nonlinKm 2e4`) a deposit cell on a steep slope
    exported more than its deposit in one sub-step, cut below its bed, and the
    `>= 0` clip created volume: `Marine diffusion mass rescale` read 0.36-0.53
    on a 10 km global run even with the donor face diffusivity. Here a 5 m
    deposit with `nlK = 5e6` came out at 3.27x its volume without the limiter.
    """
    from mpi4py import MPI

    model = minimal_model
    if getattr(model, "flatModel", False):
        pytest.skip("needs a sphere fixture with submarine relief")
    try:
        _steep_marine_setup(model, nlK=5.0e6)
        owned = model.inIDs == 1
        dh = np.zeros(model.lpoints)
        dh[model.seaID] = 5.0
        model.marineSolver = "picard"
        model.picardSub, model.picardIts = 40, 4
        vin = MPI.COMM_WORLD.allreduce(float(np.sum((dh * model.larea)[owned])))
        nd = model._diffuseImplicit(dh, model.seaID, model.nlK, label="marine")
        vout = MPI.COMM_WORLD.allreduce(float(np.sum((nd * model.larea)[owned])))
    finally:
        model.destroy()

    assert abs(vout / vin - 1.0) < 1.0e-3, (
        f"Picard diffusion volume ratio {vout / vin:.4f} != 1 at a large "
        f"diffusion number: deposit cells drain below their bed (supply limiter)."
    )


def test_marine_ts_failure_falls_back_to_picard(minimal_model):
    """
    Protects: `_diffuseImplicit` falls back to the Picard solver when the
    adaptive TS gives up (`DIVERGED_STEP_REJECTED`), instead of raising
    PETSc error 91 and aborting the run (seen on a 10 km global run with
    `marineSolver: ts`). The fallback is sticky (`marineSolver` -> picard).
    """
    model = minimal_model
    if getattr(model, "flatModel", False):
        pytest.skip("needs a sphere fixture with submarine relief")
    try:
        bed = _steep_marine_setup(model, dt=5.0e5, nlK=5.0e6)
        dh = np.zeros(model.lpoints)
        dh[model.seaID & (bed > 1000.0)] = 500.0
        model._diffuseImplicit(dh, model.seaID, model.nlK, label="marine")
        assert model.marineSolver == "ts"
        # Any rejected step now exhausts the TS.
        model._ts_marine.setMaxStepRejections(0)
        nd = model._diffuseImplicit(dh, model.seaID, model.nlK, label="marine")
        solver = model.marineSolver
    finally:
        model.destroy()

    assert solver == "picard", "TS failure did not switch to the Picard solver"
    assert np.all(np.isfinite(nd)) and nd.min() >= 0.0


def test_ofill_is_relative_to_sea_level(minimal_model):
    """
    Protects: `oFill` means the SAME thing everywhere — a depth below sea level.
    `pitfilling.fillElevation` / `fillIceElevation` compute the depression-fill
    cut-off as `max(minh, sealevel + oFill)`, and `seaplex._matOcean` must hand
    `epsfill` the same `sealevel + oFill` (recorded as `_oceanFillLevel`).

    Silent failure prevented: `_matOcean` used to treat `self.oFill` as an
    ABSOLUTE elevation. That is the same number only when the sea level is 0 —
    which is true of every fixture, benchmark and Earth example, so nothing here
    would ever have caught it — and it diverges by `sealevel` otherwise.

    Why it is not cosmetic: `epsfill` only eps-fills the cells AT OR ABOVE its
    cut-off (everything below is pre-flagged in its seeding loop and never
    revisited). A cut-off placed above sea level therefore skips the ENTIRE
    marine domain, so every closed bathymetric pocket on the shelf survives into
    the marine flow-direction surface and river sediment is trapped at the coast
    instead of routing basinward. A model with a deep datum hits this with a
    perfectly reasonable-looking YAML: `sea: position: -2200` with
    `oFill: -1500` gave a cut-off of -1500 m, i.e. 700 m ABOVE sea level.

    The `minimal` fixture is a global sphere with `sea: position: -100.`; the
    test sets `oFill: -5000.` so it decides the cut-off (above the smoothed
    minimum floor, -5653 m): the correct cut-off is -5100 m, the old absolute
    reading gave -5000 m.
    """
    model = minimal_model
    # Since 2026-11-12 the global cut-off floor is the SMOOTHED surface's
    # minimum + 0.1 (-5653 m on this fixture, whose raw -16 km needle smoothing
    # raises), so the default oFill (-6000 -> -6100) no longer decides it. Use
    # an oFill that does: sealevel + oFill = -5100 m (relative, correct) vs
    # -5000 m (absolute reading, the old bug).
    model.oFill = -5000.0
    # `_matOcean` runs inside `seaChange`, so the marine path has to execute
    # once. The minimal fixture is a short global-sphere run with deposition on.
    model.runProcesses()

    level = getattr(model, "_oceanFillLevel", None)
    assert level is not None, (
        "`seaplex._matOcean` no longer records `_oceanFillLevel`; the oFill "
        "contract is then untestable from here — re-add it or update this test."
    )

    expected = model.sealevel + model.oFill
    assert expected > -5653.0, "oFill must dominate the smoothed-minimum floor here"
    assert level == pytest.approx(expected), (
        f"marine eps-fill cut-off is {level} m but `sealevel + oFill` gives "
        f"{expected} m (sealevel={model.sealevel}, oFill={model.oFill})."
    )

    # Guard the specific regression: the raw (absolute) oFill must NOT be used.
    assert model.sealevel != 0.0, (
        "this fixture must keep a non-zero sea level or the test cannot "
        "distinguish the relative reading from the absolute one."
    )
    assert level != pytest.approx(model.oFill), (
        "the marine eps-fill cut-off equals the raw `oFill`, i.e. it is being "
        "read as an absolute elevation again (the pre-fix behaviour)."
    )


def _icosphere(level, radius):
    """Regular icosphere (no degenerate cells), unit triangles scaled to radius."""
    t = (1 + 5 ** 0.5) / 2
    v = [np.array(p, float) for p in (
        [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t],
        [0, -1, -t], [0, 1, -t], [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1])]
    v = [p / np.linalg.norm(p) for p in v]
    f = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4],
         [11, 10, 2], [10, 7, 6], [7, 1, 8], [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8],
         [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]]
    for _ in range(level):
        cache, nf = {}, []

        def mid(a, b):
            key = (min(a, b), max(a, b))
            if key not in cache:
                p = v[a] + v[b]
                v.append(p / np.linalg.norm(p))
                cache[key] = len(v) - 1
            return cache[key]
        for a, b, c in f:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            nf += [[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]
        f = nf
    return np.array(v) * radius, np.array(f, dtype=np.int32)


def test_marine_epsfill_seeded_on_global_mesh(tmp_path):
    """
    Protects: on a global mesh `seaplex._matOcean` takes the eps-fill cut-off
    from the SMOOTHED surface it fills (`_hillSlope(smooth=2)`).

    Silent failure prevented: the cut-off was `hGlobal.min() + 0.1` (raw
    elevation). Smoothing raises the deepest cell by more than 0.1 m whenever
    it is a sharp low, so no cell seeded `epsfill`, which then returned the
    surface unchanged: the marine routing ran on unfilled bathymetry with
    hundreds of terminal sinks, and river loads were force-deposited km-thick
    in the first shallow pit (a 23.6M-node Earth run; fixed 2026-11-12,
    docs/dev/FIXED_BUGS.md).

    Invariant: a sphere whose deepest point is a one-cell low, under shallow
    pits: epsfill is seeded (some cell is below the cut-off), it changes the
    surface, and at most a handful of terminal sinks remain.
    """
    import os

    from gospl.model import Model
    from gospl.sed import seaplex

    R = 6378137.0
    v, c = _icosphere(4, R)                      # 2562 nodes
    lat = np.arcsin(v[:, 2] / R)
    lon = np.arctan2(v[:, 1], v[:, 0])
    z = -3000.0 + 2500.0 * np.cos(lat) * np.cos(lon)          # a continent around lon 0
    rng = np.random.default_rng(3)
    ocean = np.where(z < -500.0)[0]
    pits = rng.choice(ocean, 40, replace=False)
    z[pits] -= 300.0                                          # shallow closed pockets
    deep = ocean[np.argmin(z[ocean])]
    z[deep] -= 1500.0                                         # sharp one-cell deepest point
    np.savez(tmp_path / "sph.npz", v=v, c=c, z=z)
    (tmp_path / "sph.yml").write_text("""name: marine epsfill seeding
domain:
    npdata: ['sph','v','c','z']
    flowdir: 5
    seadepo: True
time:
    start: 0.
    end: 1000.
    tout: 1000.
    dt: 1000.
spl:
    K: 1.e-6
    d: 0.
    m: 0.5
diffusion:
    hillslopeKa: 0.1
    hillslopeKm: 1.0
sea:
    position: 0.
climate:
  - start: 0.
    uniform: 1.
output:
    dir: 'out'
""")
    seen = {}
    orig = seaplex.epsfill

    def spy(minh, elev):
        out = orig(minh, elev)
        seen["seeds"] = int((np.asarray(elev) < minh).sum())
        seen["changed"] = int((np.asarray(out) != np.asarray(elev)).sum())
        return out

    cwd = os.getcwd()
    os.chdir(tmp_path)
    seaplex.epsfill = spy
    model = None
    try:
        model = Model("sph.yml", verbose=False)
        model.runProcesses()
        nsink = int((model.is_sink_local & (model.inIDs == 1)).sum())
    finally:
        seaplex.epsfill = orig
        if model is not None:
            model.destroy()
        os.chdir(cwd)
    assert seen.get("seeds", 0) > 0, "epsfill got no seed below its cut-off (no fill at all)"
    assert seen["changed"] > 0, "epsfill left the marine routing surface unchanged"
    assert nsink <= 5, f"{nsink} terminal sinks left on the marine routing surface"


def test_coupled_deposition_block_solves_are_direct(minimal_model):
    """
    Protects: the coupled two-field deposition solves (`SPL._coupledEDSystem`,
    `seaplex._depMarineSystem`) apply their block preconditioners once per
    outer iteration (`fieldsplit_{h,q}_ksp_type preonly`), set through the
    options database.

    Silent failure prevented: the `subksps[i].setType("preonly")` calls ran
    before `setFromOptions` switched the split to Schur, so they never took
    effect: both blocks ran GMRES (rtol 1e-5, up to 10,000 its) inside every
    outer iteration. Results were unchanged, so only the cost showed: the SPL
    solve grew from ~300 s to ~5,000 s per step on a 23.6M-node mesh (fixed
    2026-11-12, docs/dev/FIXED_BUGS.md).

    Invariant: with G > 0 and Gmar > 0 on the minimal sphere, no GMRES
    orthogonalisation happens inside either coupled solve (outer TFQMR,
    preonly blocks), and both solves run.
    """
    import sys

    import petsc4py

    PETSc = petsc4py.PETSc
    import gospl.eroder.SPL  # noqa: F401
    from gospl.sed import seaplex
    SPLmod = sys.modules["gospl.eroder.SPL"]

    model = minimal_model
    model.fDepa = 0.1
    model.Gmar = 0.5
    PETSc.Log.begin()
    orth = PETSc.Log.Event("KSPGMRESOrthog")
    counts = {}

    def counted(name, f):
        def g(self, *a, **k):
            before = orth.getPerfInfo().get("count", 0)
            out = f(self, *a, **k)
            counts.setdefault(name, []).append(orth.getPerfInfo().get("count", 0) - before)
            return out
        return g

    spl_orig = SPLmod.SPL.__dict__["_coupledEDSystem"]
    sea_orig = seaplex.SEAMesh.__dict__["_depMarineSystem"]
    SPLmod.SPL._coupledEDSystem = counted("spl", spl_orig)
    seaplex.SEAMesh._depMarineSystem = counted("marine", sea_orig)
    try:
        model.runProcesses()
    finally:
        SPLmod.SPL._coupledEDSystem = spl_orig
        seaplex.SEAMesh._depMarineSystem = sea_orig
    assert counts.get("spl"), "the coupled SPL deposition solve never ran (G > 0)"
    assert counts.get("marine"), "the coupled marine deposition solve never ran (Gmar > 0)"
    assert max(counts["spl"]) == 0, f"GMRES iterations inside the SPL coupled solve: {counts['spl']}"
    assert max(counts["marine"]) == 0, f"GMRES iterations inside the marine solve: {counts['marine']}"


def test_routing_pond_is_deposited_in_place(minimal_model):
    """
    Protects: sedplex._keepRoutingPond. When a cascade routing solve is
    repaired (flowplex._pondRoutingRegion), the volume its local sinks trap
    must join the closed-sink deposit (minus cells already counted as pit
    sinks), not vanish.
    """
    from mpi4py import MPI

    m = minimal_model
    try:
        m.flowAccumulation()
        n_owned = m.hGlobal.getLocalSize()
        pond = np.zeros(n_owned)
        pond[: min(5, n_owned)] = 2.0              # owned rows of a global Vec
        m._routingPond = pond
        m._closedDepo = np.zeros(m.lpoints)
        m._keepRoutingPond(1.0)
        got = MPI.COMM_WORLD.allreduce(float(m._closedDepo.sum()))
        # expected: the ponded rows mapped to local, owned, non-sink nodes
        m._pondG.setArray(pond)
        m.dm.globalToLocal(m._pondG, m._pondL)
        loc = m._pondL.getArray() * m.inIDs
        want = MPI.COMM_WORLD.allreduce(float(loc[~m.lsink].sum()))
        assert want > 0.0
        assert abs(got - want) < 1e-12 * want
        assert m._routingPond is None
    finally:
        m.destroy()
