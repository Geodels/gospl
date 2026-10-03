"""
Flow routing: flow-accumulation KSP, pit graph, downstream cascade, evaporation.

Protects: AGENTS.md > Flow-accumulation KSP, Input-mesh contract.

Split out of the former tests/test_regression.py. Run this group alone with
`pytest -m flow`; see tests/README.md for the marker list.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest

# Skip the whole module (rather than erroring at collection) when the
# goSPL runtime stack is not installed.
pytest.importorskip(
    "gospl.tools.inputparser",
    reason="goSPL runtime deps (petsc4py / ruamel.yaml / scipy) not installed",
)

pytestmark = [pytest.mark.flow]


def test_pit_unifyLabels_unionfind():
    """
    Protects: PITFill._unifyLabels — the union-find that collapses cross-rank
    depression-label equivalence pairs to one canonical id per connected
    component (replacing the iterative sort_ids fixpoint). Guarantees:
      - canonical id is the component MINIMUM (deterministic, partition-free),
      - multi-hop chains and merges-via-shared-node fully collapse,
      - labels not in any pair (and the -1 border sentinel) pass through.
    Pure function (self unused), so call it unbound with no Model/mesh.
    """
    pitfilling = pytest.importorskip("gospl.flow.pitfilling")
    PITFill = pitfilling.PITFill

    def unify(pairs, label):
        df = (pd.DataFrame(pairs, columns=["p1", "p2"]) if pairs
              else pd.DataFrame({"p1": [], "p2": []}))
        return PITFill._unifyLabels(None, df, np.asarray(label, dtype=int))

    # multi-hop chain 1-4-7 collapses to min (1); separate comp {2,6}->2;
    # merge via shared node 8: {3,5,8}->3; lone label 9 and -1 unchanged.
    label = np.array([1, 4, 7, 2, 6, 3, 5, 8, 9, -1])
    expected = np.array([1, 1, 1, 2, 2, 3, 3, 3, 9, -1])
    out = unify([(1, 4), (4, 7), (2, 6), (3, 8), (5, 8)], label)
    assert np.array_equal(out, expected), f"{out} != {expected}"

    # order/direction independence: shuffled pairs (some p1>p2) give the same
    # canonical mapping — union-find is symmetric and roots to the component min.
    out2 = unify([(5, 8), (8, 3), (7, 4), (4, 1), (6, 2)], label)
    assert np.array_equal(out2, expected), f"{out2} != {expected}"

    # empty pair set is a no-op
    assert np.array_equal(unify([], label), label)


# ---------------------------------------------------------------------------
# TEST 2c - channel-evap hook reduces FA (DESIGN_EVAPORATION.md §4 T1)
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_evap_reduces_FA(minimal_model, minimal_model_with_evap):
    """
    Protects: DESIGN_EVAPORATION.md §2.1 — channel-evap hook subtracts
    evap from rainA at flowplex.py:429 before the IDA solve. Because FA
    is linear in the IDA right-hand side, halving rainfall via evap
    should roughly halve FA.

    Silent failure prevented: a future refactor that drops the
    flowplex.py hook, or reverses its sign, would leave FA unchanged
    when evap is enabled. evapLoss would also stay at zero, no error.

    The two model fixtures load the same `minimal.yml`; only the second
    has uniform evap injected = 50% of the uniform rain rate (so the
    expected reduction is ~50%, allowing tolerance for the IDA solver
    and any non-linearities introduced by pit filling).
    """
    # Baseline run — no evap.
    minimal_model.runProcesses()
    fa_baseline = float(minimal_model.FAG.sum())
    assert minimal_model.evapLoss == 0.0, (
        "evapLoss should be zero when evapdata is None"
    )
    assert fa_baseline > 0, "Baseline FA should be positive"

    # With-evap run — same YAML + injected uniform evap.
    minimal_model_with_evap.runProcesses()
    fa_with_evap = float(minimal_model_with_evap.FAG.sum())
    assert minimal_model_with_evap.evapLoss > 0, (
        "evapLoss should accumulate when channel-evap hook fires; got "
        f"{minimal_model_with_evap.evapLoss}"
    )
    # Generous bound: 50% evap should drop FA below 70% of baseline.
    assert fa_with_evap < 0.7 * fa_baseline, (
        f"channel-evap did not meaningfully reduce FA. "
        f"baseline={fa_baseline:.3e}, with_evap={fa_with_evap:.3e}, "
        f"ratio={fa_with_evap / fa_baseline:.3f}"
    )


# ---------------------------------------------------------------------------
# TEST 2d - lake-evap dominates → no lake forms (DESIGN_EVAPORATION.md §4 T2)
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_lake_not_formed_when_evap_dominates(minimal_model):
    """
    Protects: DESIGN_EVAPORATION.md §2.2 — lake-evap hook in
    `_distributeDownstream` AND the partial-fill mask refinement at
    flowplex.py:346 (the `(inV > 0)` clause).

    With evap = 100x rain, the combined effect of channel-evap and
    lake-evap hooks must prevent every lake from forming. `waterFilled`
    is "depth above hl" after `flowAccumulation` (line 486), so a missing
    lake means waterFilled == 0 everywhere.

    Two silent failures this guards against:
      1. The lake-evap hook is missing — water would reach pits, eV < 0,
         partial-fill branch would raise waterFilled to a positive value.
      2. The mask refinement at line 346 is missing — pits with
         inV=0-after-evap would erroneously go through the partial-fill
         branch and end up with waterFilled bumped by the epsilon nudge
         from pitfilling.py:567 (~1e-3 m).
    """
    import pandas as pd

    if minimal_model.raindata is None:
        pytest.skip("minimal.yml has no rainfall")
    rUni = minimal_model.raindata.at[0, "rUni"]
    if pd.isnull(rUni):
        pytest.skip("minimal.yml rain is not uniform")

    # Massive evap: 100x rain. Channel-evap consumes most cells'
    # rainfall; whatever survives meets lake-evap's max-fill budget.
    minimal_model.evapdata = pd.DataFrame(
        [{"start": 0.0, "eUni": 100.0 * float(rUni),
          "eMap": None, "eKey": None}]
    )
    minimal_model.evapNb = -1
    minimal_model.evapVal = None
    minimal_model.evapMesh = None
    minimal_model.evapLoss = 0.0

    minimal_model.runProcesses()

    # The hooks fired (evapLoss accumulated).
    assert minimal_model.evapLoss > 0, (
        "evapLoss should be > 0 with massive evap; got "
        f"{minimal_model.evapLoss}"
    )

    # No lake anywhere. Tolerance 1e-6 m is well below the 1e-3 m
    # epsilon nudge that the partial-fill bug would produce.
    max_depth = float(np.max(minimal_model.waterFilled))
    assert max_depth < 1e-6, (
        "Lakes should not form with massive evap. Max waterFilled "
        f"depth = {max_depth:.3e} m. If close to 1e-3, the "
        "(inV > 0) mask refinement at flowplex.py:346 may be missing."
    )


# ---------------------------------------------------------------------------
# TEST 2e - water balance with evap (DESIGN_EVAPORATION.md §4 T3)
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_water_balance_with_evap(minimal_model):
    """
    Protects: water-mass conservation through the evap accumulator.

    Strategy: with evap = 100x rain, the channel-evap clamp at
    `min(rainA, evapVal * larea)` consumes every cell's runoff exactly
    once (lake-evap then sees inV == 0 and contributes nothing). So
    `self.evapLoss` after the run must equal the total water that
    entered the system, accurate to floating-point.

    Total input = (rainVal × larea, summed over owned land cells) ×
    (tEnd − tStart). The fixture must use uniform rain (constant in
    space) and a steady time-series (constant in time over the run) for
    this prediction to hold; the test skips otherwise.

    Silent failure prevented: an off-by-one in the `self.dt` multiplier
    inside the channel-evap hook (e.g. multiplying twice or forgetting
    altogether) would skew `evapLoss` by a factor of `dt` or `1/dt`
    relative to input. A wrong sign on `channelLoss` would underflow
    to zero. Both would fail this assertion immediately.

    Edge cases NOT covered:
      - Mixed channel + lake evap under moderate rates — requires a
        fixture with a known closed depression and tuned rates to
        guarantee the lake-evap hook fires. Out of scope for v1.
      - Time-varying evap (multiple climate rows). Requires per-step
        bookkeeping the current accumulator does not expose.
    """
    import pandas as pd
    from mpi4py import MPI

    if minimal_model.raindata is None:
        pytest.skip("minimal.yml has no rainfall")
    rUni = minimal_model.raindata.at[0, "rUni"]
    if pd.isnull(rUni):
        pytest.skip("minimal.yml rain is not uniform; need scalar rain for budget")
    if len(minimal_model.raindata) > 1:
        pytest.skip("minimal.yml has time-varying rain; T3 needs steady forcing")

    # Massive evap: channel-evap clamp consumes all rain at every cell
    # every step. Lake-evap never fires (inV = 0 after subtraction).
    minimal_model.evapdata = pd.DataFrame(
        [{"start": 0.0, "eUni": 100.0 * float(rUni),
          "eMap": None, "eKey": None}]
    )
    minimal_model.evapNb = -1
    minimal_model.evapVal = None
    minimal_model.evapMesh = None
    minimal_model.evapLoss = 0.0

    # Drive a single flowAccumulation call directly. We avoid runProcesses
    # because (a) model.py:240 and 255 call flowAccumulation TWICE per
    # step (pre-SPL and pre-sedChange), each firing the evap hook with a
    # different mid-step seaID — the integrated budget over both calls
    # isn't predictable from a single final-state seaID snapshot; (b)
    # applyForces fires at the END of each iteration (model.py:273-274),
    # so the FIRST step's flowAcc would skip the hook entirely on
    # un-primed forcing. Calling these two methods explicitly accumulates
    # exactly one step's worth of evap against a stable seaID, which is
    # the smallest case that proves the unit-conversion and accumulator
    # math is sound.
    minimal_model.applyForces()       # populate rainVal, evapVal, sealevel
    minimal_model.flowAccumulation()  # accumulate evap once; set seaID

    # Total water input (m^3) over the run = ∫∫ rainVal × dA × dt.
    # Constant uniform rain → factor out: rainVal × land_area × duration.
    rainVal = minimal_model.rainVal
    larea = minimal_model.larea
    owned = minimal_model.inIDs == 1
    is_land = np.ones(len(rainVal), dtype=bool)
    is_land[minimal_model.seaID] = False

    # One flowAcc call accumulates evap over exactly one dt.
    duration = float(minimal_model.dt)
    rain_rate_local = float(np.sum((rainVal * larea)[owned & is_land]))
    input_local = rain_rate_local * duration
    input_total = MPI.COMM_WORLD.allreduce(input_local, op=MPI.SUM)
    rain_rate_total = MPI.COMM_WORLD.allreduce(rain_rate_local, op=MPI.SUM)

    if input_total < 1.0:
        pytest.skip(
            f"Total water input {input_total:.3e} m^3 below 1 m^3 — "
            f"fixture too small for a meaningful balance check"
        )

    # evapLoss is rank-local (each rank accumulates its own cells).
    evap_total = MPI.COMM_WORLD.allreduce(
        float(minimal_model.evapLoss), op=MPI.SUM
    )

    rel_error = abs(evap_total - input_total) / input_total
    diagnostic = (
        f"\n  dt          = {minimal_model.dt}"
        f"\n  duration    = {duration} (one flowAcc call)"
        f"\n  rain_rate   = {rain_rate_total:.6e} m^3/yr (global)"
        f"\n  input_total = {input_total:.6e} m^3"
        f"\n  evap_total  = {evap_total:.6e} m^3"
        f"\n  rel_error   = {rel_error:.3e}"
        f"\n  evap/input  = {evap_total / input_total:.6f}"
    )
    assert rel_error < 1.0e-4, (
        "Water mass not conserved through evap accumulator." + diagnostic
        + "\nCheck the channel-evap hook's `* self.dt` factor at "
        + "flowplex.py and the lakeLoss accumulation in "
        + "_distributeDownstream."
    )


# ---------------------------------------------------------------------------
# TEST 2f - losing-stream evaporation (debited from accumulated discharge)
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_evap_losing_stream(minimal_model):
    """
    Protects: the opt-in "losing stream" mode (`flowplex._losingStreamSolve`) —
    evaporation debited from the **accumulated discharge**, not just the local
    runoff. With high evaporation downstream (lower-elevation half) and none
    upstream, the source-side hook can only remove the small LOCAL runoff there
    (the through-flowing river keeps its discharge), whereas the losing-stream
    hook evaporates the through-flow and strongly reduces downstream FA. Both
    keep FA >= 0 and conserve water (evaporated <= total rain that fell).
    """
    import numpy as np
    from mpi4py import MPI

    m = minimal_model
    m.applyForces()
    rainVal = m.rainVal.copy()
    larea = m.larea
    owned = m.inIDs == 1
    rain_total = MPI.COMM_WORLD.allreduce(
        float(np.sum((rainVal * larea)[owned])), op=MPI.SUM)
    if rain_total < 1.0:
        pytest.skip("fixture too small for a meaningful losing-stream check")

    # Baseline (no evap) to locate the actual channel cells (high discharge);
    # evaporating LAKE-interior cells would do nothing (their FA is 0), so the
    # losing-stream effect must be tested on the through-flowing channels.
    m.evapVal = None
    m.evapStream = False
    m.evapLoss = 0.0
    m.flowAccumulation()
    fa0 = m.FAL.getArray().copy()
    fa_base = float(m.FAG.sum())
    pos = fa0[owned & (fa0 > 0.0)]
    if pos.size < 5:
        pytest.skip("no channels in the fixture")
    thr = np.percentile(pos, 90)            # top-decile discharge = main channels
    evap = np.where(fa0 > thr, 50.0 * rainVal, 0.0)   # strong channel evaporation
    evap[m.seaID] = 0.0

    # Source-side (current default): evaporates only the LOCAL runoff of those
    # cells, so the through-flowing river keeps most of its discharge.
    m.applyForces()
    m.evapVal = evap.copy()
    m.evapStream = False
    m.evapLoss = 0.0
    m.flowAccumulation()
    fa_src = float(m.FAG.sum())
    loss_src = MPI.COMM_WORLD.allreduce(float(m.evapLoss), op=MPI.SUM)

    # Losing-stream: evaporates the ACCUMULATED discharge passing through.
    m.applyForces()                       # reset the runoff source (bL)
    m.evapVal = evap.copy()
    m.evapStream = True
    m.evapLoss = 0.0
    m.flowAccumulation()
    fa_los = float(m.FAG.sum())
    loss_los = MPI.COMM_WORLD.allreduce(float(m.evapLoss), op=MPI.SUM)
    famin = MPI.COMM_WORLD.allreduce(float(m.FAL.getArray().min()), op=MPI.MIN)

    assert famin >= -1.0e-9, "losing-stream FA went negative"
    assert loss_los > 0.0, "no evaporation accumulated in losing-stream mode"
    # Can't evaporate more water than fell as rain (conservation bound).
    assert loss_los <= rain_total * m.dt * (1.0 + 1.0e-6)
    # Losing-stream removes MORE than source-side and cuts downstream FA further.
    assert loss_los > loss_src
    assert fa_los < fa_src < fa_base


def _singular_flow_system(ncycles, rhs_all=False):
    """Build a standalone singular `(I - W^T)`-like system for the fatal-solve
    recovery test. N=1000 identity matrix with `ncycles` disjoint 2-cycles
    (rows 2k, 2k+1): the 2x2 block ``[[1,-1],[-1,1]]`` is singular (null vector
    [1,1]), so Richardson cannot converge on a block whose RHS has a component
    along that null space -- exactly the knife-edge un-drainable configuration
    the fatal flow solve hits. Returns (A, b, x) with x seeded to b.
    """
    import petsc4py
    PETSc = petsc4py.PETSc
    N = 1000
    A = PETSc.Mat().create(PETSc.COMM_WORLD)
    A.setSizes(((None, N), (None, N)))
    A.setType("aij")
    A.setPreallocationNNZ((2, 2))
    A.setUp()
    istart, iend = A.getOwnershipRange()
    for gi in range(istart, iend):
        A.setValue(gi, gi, 1.0)
    for k in range(ncycles):
        i, j = 2 * k, 2 * k + 1
        if istart <= i < iend:
            A.setValue(i, j, -1.0)
        if istart <= j < iend:
            A.setValue(j, i, -1.0)
    A.assemble()
    b = A.createVecLeft()
    b.set(0.0)
    if rhs_all:
        for k in range(ncycles):
            i = 2 * k
            if istart <= i < iend:
                b.setValue(i, 1.0)
    elif istart <= 0 < iend:
        b.setValue(0, 1.0)
    b.assemble()
    x = b.duplicate()
    b.copy(result=x)
    return A, b, x


def test_fatal_flow_solve_ponds_small_undrained_region(minimal_model):
    """
    Protects: the fatal flow-accumulation discharge solve must RECOVER from a
    small, knife-edge un-drained region (pond it and continue) but still ABORT
    on a genuinely broken (large / non-finite) system
    (flowplex._solve_KSP2, fatal=True).

    Silent failure prevented: a near-flat region whose MFD tie-break forms a
    near-cycle makes `(I - W^T)` locally singular over O(100) cells. That is a
    floating-point knife-edge -- a restart's float32 elevation truncation, or a
    different partition, perturbs it away -- yet the fatal main solve used to
    abort the whole (potentially many-hour) run over it. It now ponds a small
    finite region (discharge -> local runoff `b`, plus a hard clamp on any cell
    whose discharge exceeds the total domain runoff, which no physical cell can)
    and continues, exactly as the non-fatal benign path does; it still aborts a
    large region or a non-finite RHS/matrix (a real broken state).

    The mesh cannot reproduce the knife-edge deterministically, so this drives
    `_solve_KSP2(fatal=True)` directly with a standalone singular matrix.
    """
    from mpi4py import MPI

    model = minimal_model
    cap = int(model._undrained_benign_cap)

    # --- SMALL region (1 singular cell): must POND, not abort ---
    A, b, x = _singular_flow_system(ncycles=1)
    try:
        raised = False
        try:
            model._solve_KSP2(A, b, x, fatal=True)
        except RuntimeError:
            raised = True
        total_runoff = float(b.sum())
        xmax = x.max()[1]
        xfinite = bool(np.isfinite(x.norm()))
    finally:
        A.destroy(); b.destroy(); x.destroy()

    assert not raised, (
        "fatal flow solve aborted on a SMALL un-drained region; it should pond "
        "those cells and continue (knife-edge routing degeneracy recovery)."
    )
    assert xfinite, "ponded discharge is not finite."
    # Every cell ponded/clamped to a physical value: no cell drains more than
    # the total domain runoff (the null-space blow-up is removed).
    assert xmax <= total_runoff + 1.0e-6, (
        f"ponded discharge {xmax:.3e} exceeds total runoff {total_runoff:.3e}: "
        f"the null-space blow-up was not clamped."
    )

    # --- LARGE region (> benign cap singular cells): must ABORT ---
    A, b, x = _singular_flow_system(ncycles=cap + 50, rhs_all=True)
    try:
        raised = False
        try:
            model._solve_KSP2(A, b, x, fatal=True)
        except RuntimeError:
            raised = True
    finally:
        A.destroy(); b.destroy(); x.destroy()

    assert raised, (
        f"fatal flow solve did NOT abort on a LARGE un-drained region "
        f"(> {cap} cells); a genuinely broken discharge must still abort."
    )


def test_undrained_cap_env_override(minimal_model):
    """
    Protects: `GOSPL_UNDRAINED_CAP` resolves the un-drained-region cap
    (`flowplex._undrainedCap`) as a fraction (< 1) or an absolute node count
    (>= 1), and rejects a malformed value instead of silently reverting.

    Why the knob exists: the default cap is `max(256, 0.5% of mpoints)`, tuned
    so a knife-edge micro-cycle ponds and a broken partition still aborts. A
    wide, genuinely closed near-flat basin (an endorheic interior at fine
    resolution) can make `(I - W^T)` singular over far more cells than that,
    and the FATAL main discharge solve then aborts a long run over a region
    that physically just ponds. This lets the cap be raised deliberately. It
    can NOT mask a NaN source or a broken operator: those trip the separate
    non-finite RHS/matrix checks and abort at any cap.

    The value gates COLLECTIVE branches (`benign` / `pond_fatal`), so it is
    broadcast from rank 0 and must be identical everywhere.
    """
    model = minimal_model
    mpoints = int(model.mpoints)
    default = max(256, int(0.005 * mpoints))

    assert model._undrainedCap() == default, (
        "with the variable unset the cap must be the 0.5%-of-mesh default."
    )

    previous = os.environ.get("GOSPL_UNDRAINED_CAP")
    try:
        # Absolute node count (>= 1).
        os.environ["GOSPL_UNDRAINED_CAP"] = "100000"
        assert model._undrainedCap() == 100000

        # Fraction of the mesh (< 1).
        os.environ["GOSPL_UNDRAINED_CAP"] = "0.02"
        assert model._undrainedCap() == int(round(0.02 * mpoints))

        # A fraction small enough to round below one cell still leaves a
        # usable cap: a zero cap would make EVERY solve failure fatal.
        os.environ["GOSPL_UNDRAINED_CAP"] = "1.0e-12"
        assert model._undrainedCap() == 1

        # Malformed / non-positive values must fail loudly, not revert to the
        # default -- silently ignoring the knob would read as "I raised the cap
        # and it still aborted".
        for bad in ("not-a-number", "-5", "0"):
            os.environ["GOSPL_UNDRAINED_CAP"] = bad
            with pytest.raises(ValueError):
                model._undrainedCap()
    finally:
        if previous is None:
            os.environ.pop("GOSPL_UNDRAINED_CAP", None)
        else:
            os.environ["GOSPL_UNDRAINED_CAP"] = previous

    # The model's own cap is untouched by the probing above.
    assert model._undrained_benign_cap == default


# ---------------------------------------------------------------------------
# TEST 12 - Downstream cascade: relative residual floor vs stagnation break
# ---------------------------------------------------------------------------
#
# Pure-logic test (no Model, no PETSc) — instantiates FAMesh via `__new__`
# to bypass the heavy init and exercises the cascade progress tracker
# directly. Runs in microseconds; belongs to the fast tier.
# ---------------------------------------------------------------------------


def _cascade_tracker(rel_floor=1.0e-3, patience=3, rel_improve=1.0e-3):
    """
    Bare FAMesh carrying only the cascade-tracker state (no mesh, no PETSc).
    """
    from gospl.flow.flowplex import FAMesh

    fa = FAMesh.__new__(FAMesh)
    fa._cascade_rel_floor = rel_floor
    fa._cascade_patience = patience
    fa._cascade_rel_improve = rel_improve
    fa._cascadeResetProgress()
    return fa


def test_cascade_relative_residual_floor():
    """
    Protects: `flowplex._cascadeStopReason` — the OUTER downstream-routing loop
    must stop once the residual flux is a negligible fraction of its initial
    value ("floor", calm), and must keep reporting a genuine plateau loudly
    ("stall"), and must be able to fall back to the pre-floor behaviour.

    Silent failure prevented: the only other cheap exit from that loop is
    ABSOLUTE (`_distributeDownstream` skips the solve below `maxarea`, i.e. ~1 m
    of water over the largest cell). On a coarse mesh that threshold is reached
    only after the cascade has chased a residual already 5-6 orders of magnitude
    below where it started, each extra pass paying a full `fMat` rebuild + KSP
    solve to route a physically irrelevant trickle (measured on a 10 km / np=8
    run: 5 of 18 passes, ~17% of the flow-accumulation wall-time, for 0.005% of
    the cascade flux). Regressing the floor silently restores that cost; getting
    its precedence wrong instead re-labels a normal completion as the loud
    "cascade stalled ... un-drainable" warning (or vice versa, hiding a real
    un-drainable pocket behind a calm message).

    Every input is a global `Vec.sum` in production, so the decision is a pure
    function of globally-identical scalars — that is what keeps the caller's
    `break` collective-consistent, and it is why this can be tested without MPI.
    """
    # --- 1. A steadily-shrinking residual runs until it crosses the floor ---
    fa = _cascade_tracker(rel_floor=1.0e-3)
    assert fa._cascadeStopReason(1.0e13) is None, (
        "the first pass only seeds the baseline; it can never stop the loop."
    )
    assert fa._cascadeStopReason(1.0e12) is None, "1e-1 of initial: keep routing."
    assert fa._cascadeStopReason(1.0e11) is None, "1e-2 of initial: keep routing."
    assert fa._cascadeStopReason(5.0e9) == "floor", (
        "residual fell to 5e-4 of the initial flux (below the 1e-3 floor) but "
        "the cascade did not take the calm converged-enough exit."
    )

    # --- 2. A plateau at HIGH residual is still the loud stagnation break ---
    fa = _cascade_tracker(rel_floor=1.0e-3, patience=3)
    assert fa._cascadeStopReason(3.79e12) is None
    reasons = [fa._cascadeStopReason(3.79e12) for _ in range(3)]
    assert reasons == [None, None, "stall"], (
        f"a frozen residual must trip the stagnation break exactly on the "
        f"patience-th non-improving pass, got {reasons}."
    )

    # --- 3. The floor OUTRANKS the stagnation break (calm beats loud) ---
    # A residual that is both negligible AND no longer shrinking is a completed
    # cascade, not an un-drainable pocket: it must NOT raise the loud warning.
    fa = _cascade_tracker(rel_floor=1.0e-3, patience=1)
    assert fa._cascadeStopReason(1.0e13) is None
    assert fa._cascadeStopReason(1.0e6) == "floor"
    assert fa._cascadeStopReason(1.0e6) == "floor", (
        "a negligible, plateaued residual was reported as a stalled "
        "un-drainable pocket; the floor must take precedence."
    )

    # --- 4. rel_floor = 0 restores the pre-floor (stagnation-only) behaviour ---
    fa = _cascade_tracker(rel_floor=0.0)
    seq = [1.0e13, 1.0e12, 1.0e11, 1.0e10, 1.0e-6]
    assert [fa._cascadeStopReason(r) for r in seq] == [None] * len(seq), (
        "with GOSPL_CASCADE_REL_FLOOR=0 the loop must behave exactly as before "
        "the floor existed (only the stagnation break and the absolute "
        "maxarea/step-cap exits stop it)."
    )

    # --- 5. The baseline is per-cascade: reset must re-arm the floor ---
    # flowAccumulation runs the cascade twice per timestep (water, then again
    # after erosion) and sedplex runs it for sediment. A second cascade legiti-
    # mately STARTS with a small residual; if the baseline leaked across calls
    # it would "floor" on its very first comparison and skip real routing.
    fa = _cascade_tracker(rel_floor=1.0e-3)
    fa._cascadeStopReason(1.0e13)
    fa._cascadeStopReason(1.0e9)          # would floor
    fa._cascadeResetProgress()
    assert fa._cascade_resid0 is None, "reset must clear the residual baseline."
    assert fa._cascadeStopReason(1.0e9) is None, (
        "after a reset the next cascade must re-baseline on its own first pass, "
        "not inherit the previous cascade's (much larger) initial residual."
    )
    assert fa._cascadeStopReason(1.0e8) is None, (
        "1e-1 of the NEW baseline must keep routing."
    )


def _fan_mesh(degree):
    """A single central vertex of the requested degree, ringed by triangles."""
    angle = np.linspace(0.0, 2.0 * np.pi, degree, endpoint=False)
    coords = np.zeros((degree + 1, 3), dtype=np.float64)
    coords[1:, 0] = np.cos(angle)
    coords[1:, 1] = np.sin(angle)
    cells = np.array(
        [[0, 1 + i, 1 + (i + 1) % degree] for i in range(degree)], dtype=np.int32
    )
    return coords, cells


def test_globalngbhs_reports_over_degree_vertices():
    """
    Protects: `globalngbhs` must NOT write past the fixed 12-slot global
    neighbour table `FVgnID(nt, 12)`, and must tell the caller when a vertex
    could not be stored in full.

    Silent failure prevented: the routine used to append with
    `nb = FVgnNb(n) + 1; FVgnID(n, nb) = ...` and no bound check, so a vertex
    of degree > 12 wrote outside the array. That is the `faceVel` bug class
    (see AGENTS.md > Intentional surprises): a heap overflow that macOS
    tolerates silently and glibc turns into a wandering SIGABRT at the next
    malloc, nowhere near the real culprit. goSPL's own meshes are degree <= 8,
    so only an externally generated or locally refined input mesh triggers it.

    The contract now: the excess connections are dropped rather than written,
    `nover` counts the vertices that overflowed, `FVgnNb` is clamped to 12 so
    a caller that ignores `nover` still cannot walk the table out of bounds,
    and `UnstMesh._buildMesh` raises on `nover > 0`.
    """
    fortran = pytest.importorskip(
        "gospl._fortran", reason="goSPL Fortran extension not built"
    )

    for degree in (6, 12):
        coords, cells = _fan_mesh(degree)
        assert fortran.globalngbhs(len(coords), cells) == 0, (
            f"a degree-{degree} vertex fits in the 12-slot table and must not "
            f"be reported as an overflow."
        )

    for degree in (13, 20):
        coords, cells = _fan_mesh(degree)
        nover = fortran.globalngbhs(len(coords), cells)
        assert nover == 1, (
            f"a degree-{degree} vertex overflows the 12-slot table; expected "
            f"exactly 1 over-degree vertex, got {nover}."
        )
        # The table must still be walkable: epsfill loops `k = 1, FVgnNb(i)`
        # and indexes FVgnID(i, k), so an unclamped count would read out of
        # bounds here. Elevations are strictly increasing away from the centre,
        # so the fill is a no-op and only the traversal is under test.
        elev = np.zeros(len(coords), dtype=np.float64)
        elev[0] = -1.0
        elev[1:] = 1.0 + 0.1 * np.arange(degree)
        filled = fortran.epsfill(elev.min() + 0.01, elev)
        assert np.isfinite(filled).all()
        assert (filled >= elev).all()


def _chain_flow_system(n, seed=7):
    """A well-posed but long SFD routing system: one drainage chain of `n`
    cells to an outlet, with the cell order randomly permuted (as mesh order is
    relative to drainage order). (I - W^T) is non-singular and its solution is
    the accumulated runoff along the chain, but ILU(0) in this ordering is
    inexact and fgmres(30) needs ~n iterations: the implicit_timestepping
    plateau failure (733-cell chains)."""
    import petsc4py
    PETSc = petsc4py.PETSc
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)            # perm[k] = matrix index of chain cell k
    A = PETSc.Mat().createAIJ([n, n], nnz=2, comm=PETSc.COMM_WORLD)
    for k in range(n):
        A.setValue(perm[k], perm[k], 1.0)
        if k > 0:                        # cell k-1 drains into cell k
            A.setValue(perm[k], perm[k - 1], -1.0)
    A.assemble()
    b = A.createVecLeft()
    b.set(1.0)                           # unit runoff everywhere
    x = b.duplicate()
    b.copy(result=x)
    exact = np.empty(n)
    exact[perm] = np.arange(1, n + 1, dtype=float)
    return A, b, x, exact


def _fresh_flow_ksps(m):
    """The cached flow KSPs were set up on the model's mesh-sized matrix; a
    standalone system of another size needs fresh ones (destroy_DMPlex skips
    the None attributes at teardown)."""
    for name in ("_ksp_main", "_ksp_fallback", "_ksp_exact"):
        k = getattr(m, name, None)
        if k is not None:
            k.destroy()
        setattr(m, name, None)


def test_long_chain_routing_solve_rescued_exactly(minimal_model):
    """
    Protects: flowplex._solveIDAExact. A long, well-posed routing chain that
    the primary fgmres+ILU cannot converge within the cascade cap used to be
    ZEROED by the bounded fallback (implicit_timestepping: 1949 cells' routed
    water dropped each step). The exact block-factor rescue must solve it, set
    the sticky switch, and solve the next one directly.
    """
    from mpi4py import MPI

    if MPI.COMM_WORLD.Get_size() > 1:
        pytest.skip("serial standalone matrix")
    m = minimal_model
    assert m._ida_exact is False, "fresh model must start on the primary solver"
    _fresh_flow_ksps(m)
    A, b, x, exact = _chain_flow_system(3000)
    try:
        m._solve_KSP(True, A, b, x, fatal=False, seed=True)
        got = x.getArray()
        assert np.allclose(got, exact, rtol=1e-6), (
            "routing solve not recovered: max |err| %.3e" % np.abs(got - exact).max())
        assert m._ida_exact is True, "a successful rescue must switch to exact"
        # Next routing solve goes straight to the exact solver.
        x.set(0.0)
        b.copy(result=x)
        m._solve_KSP(True, A, b, x, fatal=False, seed=True)
        assert np.allclose(x.getArray(), exact, rtol=1e-6)
        assert m._ksp_exact.getIterationNumber() <= 2
    finally:
        A.destroy(); b.destroy(); x.destroy()


def test_singular_routing_solve_still_falls_back(minimal_model):
    """
    A genuinely singular region (2-cycles) must NOT be "rescued": the exact
    solver fails too, the sticky switch stays off, and the bounded fallback
    zeroes the non-fatal solve as before.
    """
    m = minimal_model
    _fresh_flow_ksps(m)
    A, b, x = _singular_flow_system(ncycles=20, rhs_all=True)
    try:
        m._solve_KSP(True, A, b, x, fatal=False, seed=True)
        assert m._ida_exact is False
        assert np.all(np.isfinite(x.getArray()))
    finally:
        A.destroy(); b.destroy(); x.destroy()


def test_stale_warm_start_guess_is_reset(minimal_model):
    """
    Protects: flowplex._warmStartGuard. Callers that pass a scratch Vec as
    the solution start from whatever the last kernel left in it. On the first
    step of goSPL-examples glacial_erosion, _glacialMeltwater's tmp1 held a
    guess so far off that the KSP declared DIVERGED_DTOL before its first
    iteration (the system was the identity: no ice yet). A guess worse than
    zero (||b - A x0|| > ||b||) must be reset, and the solve must converge.
    """
    import petsc4py
    PETSc = petsc4py.PETSc
    m = minimal_model
    _fresh_flow_ksps(m)
    n = 2000
    A = PETSc.Mat().createAIJ([n, n], nnz=1, comm=PETSc.COMM_WORLD)
    lo, hi = A.getOwnershipRange()
    for i in range(lo, hi):
        A.setValue(i, i, 1.0)
    A.assemble()
    b = A.createVecLeft()
    b.set(3.0)
    for seed in (False, True):
        x = b.duplicate()
        x.set(1.0e15)                     # stale scratch content
        try:
            m._solve_KSP(True, A, b, x, fatal=False, seed=seed)
            assert m._ksp_main.getConvergedReason() > 0, (
                "seed=%s: stale guess not reset (reason %d)"
                % (seed, m._ksp_main.getConvergedReason()))
            assert np.allclose(x.getArray(), 3.0)
        finally:
            x.destroy()
    A.destroy(); b.destroy()
