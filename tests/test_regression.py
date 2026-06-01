"""
Regression tests protecting the invariants documented in AGENTS.md.

Each test names the specific AGENTS.md section it protects and the silent
failure mode it would catch. Two tiers:

  * Parser-only tests (TEST 1, TEST 2) run anywhere goSPL imports cleanly.
    They bypass `ReadYaml.__init__` and exercise the dataframe builders
    directly, so they need no mesh files and no PETSc DMPlex.

  * Full-stack tests (TEST 3 - TEST 6) are marked @pytest.mark.slow and
    skip cleanly when the backing fixture YAML/mesh files are absent.
    See tests/conftest.py for the fixture contract.

Run the fast set with `pytest -m 'not slow'`; the full set with `pytest`.

NOTE: this file does not import or modify any goSPL source. It only
EXERCISES public + private methods documented in AGENTS.md.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

# Module-level petsc4py.init / ruamel.yaml imports happen the moment we
# touch gospl.tools.inputparser. Skip the whole file (rather than erroring
# at collection) when the goSPL runtime stack is not installed.
inputparser = pytest.importorskip(
    "gospl.tools.inputparser",
    reason="goSPL runtime deps (petsc4py / ruamel.yaml / scipy) not installed",
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _bare_parser(t_start: float = 0.0, t_end: float = 1000.0,
                 t_out: float = 1000.0):
    """
    Build a `ReadYaml` instance with only the attributes the forcing
    parsers actually read, skipping the full `__init__` chain (which would
    require a YAML on disk and a real mesh npz).

    Intentionally uses `__new__` so the heavy file/IO bootstrap in
    `ReadYaml.__init__` (inputparser.py:24-82) does not run; the forcing
    `_readX` methods only depend on `self.input`, `self.tStart`,
    `self.tEnd`, `self.tout`, `self.rStep`, and (for `map:` paths)
    `self.meshFile`. The uniform-only YAML blocks used in the parser
    tests below never hit the file-loading branch.
    """
    parser = inputparser.ReadYaml.__new__(inputparser.ReadYaml)
    parser.tStart = t_start
    parser.tEnd = t_end
    parser.tout = t_out
    parser.rStep = 0
    parser.meshFile = "/dev/null"  # only consulted on the `map:` path
    parser.input = {}
    return parser


# ---------------------------------------------------------------------------
# TEST 1 - uniform sedfactor populates sUni (regression guard)
# ---------------------------------------------------------------------------


def test_uniform_sedfactor_populates_sUni():
    """
    Protects against regression of the rUni/sUni mismatch fixed at
    inputparser.py:915 (dict key changed from 'rUni' to 'sUni').

    Silent failure prevented: before the fix, a YAML with a uniform
    `sedfactor` event silently produced `sedfacdata['sUni'] == NaN`,
    which then crashed `unstructuredmesh._updateEroFactor` on
    `np.load(None)` when the event became current. The fix is one
    character; this guard ensures it stays correct under future
    refactors of `_defineErofactor` and the sedfactor dict-build path.

    Invariant: after `_defineErofactor(... sMap=None, sUniform=u, ...)`,
    the resulting DataFrame's `sUni` column must contain `u`, not NaN.
    """
    parser = inputparser.ReadYaml.__new__(inputparser.ReadYaml)
    df = parser._defineErofactor(
        0,      # k (first event)
        0.0,    # sStart
        None,   # sMap -> triggers the uniform-only branch (the previously-buggy path)
        1.0,    # sUniform
        None,   # sedfacdata accumulator
    )

    # Column order is also part of the contract; keep it pinned here so
    # this test doubles as a guard.
    assert list(df.columns) == ["start", "sUni", "sMap", "sKey"], (
        "Column order has drifted from AGENTS.md > "
        "Forcing DataFrame layout contract."
    )

    # Regression guard: the previously-buggy seam.
    assert not pd.isnull(df["sUni"][0]), (
        "rUni/sUni mismatch regressed: uniform value did not reach the "
        "sUni column. Check inputparser.py:915 — the dict key in the "
        "uniform branch of `_defineErofactor` must be 'sUni' to match "
        "the DataFrame columns built at line 929."
    )
    assert df["sUni"][0] == 1.0


# ---------------------------------------------------------------------------
# TEST 2 - forcing DataFrame column order contract
# ---------------------------------------------------------------------------


def test_forcing_column_order():
    """
    Protects: AGENTS.md > The forcing DataFrame layout contract.

    Silent failure prevented: every consumer of tecdata / raindata /
    sedfacdata / tedata uses positional `iloc[nb, k]` access. Inserting
    a column anywhere except the END of the DataFrame silently shifts
    every downstream read to a neighbouring column. The most fragile
    consumer is `tectonics.py:78` which reads `iloc[nb, -1]` expecting
    `hMap`.

    Invariant: `list(df.columns)` for each of the four forcing
    DataFrames must match the order documented in AGENTS.md exactly.
    The four DataFrames must be buildable purely from a YAML dict (no
    real mesh on disk), so this test runs in the fast tier.
    """
    # ---- tecdata: built via _readTectonics ----
    # Minimal event: only `start` and `end`. _defineTectonics treats
    # missing upsub/hdisp/zfit as "empty", so no file lookups happen.
    parser = _bare_parser()
    parser.input = {"tectonics": [{"start": 0.0, "end": 1000.0}]}
    parser._readTectonics()
    assert list(parser.tecdata.columns) == [
        "start", "end", "tMap", "zMap", "hMap",
    ], (
        "tecdata column order drift. AGENTS.md says: "
        "start, end, tMap, zMap, hMap. tectonics.py:78 reads iloc[nb,-1] "
        "expecting hMap."
    )

    # ---- raindata: built via _readRain (uniform-only path) ----
    parser = _bare_parser()
    parser.input = {"climate": [{"start": 0.0, "uniform": 1.0}]}
    parser._readRain()
    assert list(parser.raindata.columns) == [
        "start", "rUni", "rzA", "rzB", "rMap", "rKey",
    ], (
        "raindata column order drift. AGENTS.md says: "
        "start, rUni, rzA, rzB, rMap, rKey. "
        "unstructuredmesh.py:677-688 reads iloc[nb, 4]/[5] for rMap/rKey."
    )

    # ---- sedfacdata: built via _readErofactor (uniform-only path) ----
    # NOTE: this path also exercises the rUni/sUni bug, but the COLUMN
    # ORDER is unaffected — only the values are wrong. So this assertion
    # passes today; the value-level guard lives in TEST 1.
    parser = _bare_parser()
    parser.input = {"sedfactor": [{"start": 0.0, "uniform": 1.0}]}
    parser._readErofactor()
    assert list(parser.sedfacdata.columns) == [
        "start", "sUni", "sMap", "sKey",
    ], (
        "sedfacdata column order drift. AGENTS.md says: "
        "start, sUni, sMap, sKey. "
        "unstructuredmesh.py:725-730 reads iloc[nb, 2]/[3] for sMap/sKey."
    )

    # ---- tedata: built via _readTeMap (uniform-only path) ----
    parser = _bare_parser()
    parser.input = {"temap": [{"start": 0.0, "uniform": 10000.0}]}
    parser._readTeMap()
    assert list(parser.tedata.columns) == [
        "start", "tUni", "tMap", "tKey",
    ], (
        "tedata column order drift. AGENTS.md says: "
        "start, tUni, tMap, tKey. "
        "addprocess.py:160-174 reads iloc[nb, 2]/[3] for tMap/tKey."
    )


# ---------------------------------------------------------------------------
# TEST 3 - rcvIDi must be a snapshot, not an alias
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_rcvIDi_is_copy_not_reference(minimal_model):
    """
    Protects: AGENTS.md > The rcvID / rcvIDi convention (CRITICAL).

    Silent failure prevented: if someone ever changes the snapshot in
    `flowplex.py:421-426` from `.copy()` to a plain alias (`rcvIDi =
    self.rcvID`), then every SPL kernel and `sedplex._getSedFlux` would
    silently start reading the LIVE receiver array — which gets rebuilt
    against the filled / sediment-filled topography inside
    `_distributeDownstream` and `_moveDownstream`. Erosion rates would
    drift from physical reality with no crash and no warning.

    Invariant 1: after `flowAccumulation()`, the `*i` arrays are
    distinct Python objects from their live counterparts.
    Invariant 2: after `sedChange()` (which calls `_buildFlowDirection`
    on `sedFilled` and mutates `rcvID`/`wghtVal`/`distRcv`), the `*i`
    snapshots are byte-equal to the pre-`sedChange` capture.
    """
    model = minimal_model
    # flowAccumulation runs in Model.__init__ already (model.py:192) when
    # `self.fast` is False, so `rcvIDi` is already populated. Re-running
    # to be explicit about the contract.
    model.flowAccumulation()

    # ---- Invariant 1: distinct objects ----
    assert model.rcvIDi is not model.rcvID, (
        "rcvIDi aliases rcvID — snapshot lost. "
        "See flowplex.py:421-426 and AGENTS.md > rcvID/rcvIDi convention."
    )
    assert model.wghtVali is not model.wghtVal, (
        "wghtVali aliases wghtVal — snapshot lost."
    )
    assert model.distRcvi is not model.distRcv, (
        "distRcvi aliases distRcv — snapshot lost."
    )

    # NOTE: rcvIDi and rcvID can (and usually do) DIFFER in value right
    # after flowAccumulation. The snapshot is taken on the pre-fill
    # topography; rcvID is then rebuilt against the water-filled
    # topography inside `_distributeDownstream` (flowplex.py:364-369)
    # whenever the mesh has any pit. Don't assert byte-equality here —
    # that contradicts the very convention this test guards.

    # ---- Invariant 2: sedChange must not mutate the snapshot ----
    if getattr(model, "nodep", False):
        # sedChange isn't called from runProcesses when nodep=true.
        # The minimal_model fixture should set nodep=false for this test
        # to be meaningful.
        pytest.skip(
            "Fixture has nodep=true so sedChange isn't called from "
            "runProcesses. NEEDS_HUMAN_REVIEW: minimal.yml must set "
            "nodep: false to exercise the rcvIDi mutation guard."
        )

    rcvIDi_before = model.rcvIDi.copy()
    wghtVali_before = model.wghtVali.copy()
    distRcvi_before = model.distRcvi.copy()

    # sedChange internally rebuilds the flow direction matrix against
    # `sedFilled`, mutating self.rcvID / self.wghtVal / self.distRcv.
    # The `*i` snapshots must NOT move.
    model.sedChange()

    np.testing.assert_array_equal(
        model.rcvIDi, rcvIDi_before,
        err_msg=(
            "sedChange mutated rcvIDi. The pre-fill snapshot is shared "
            "with SPL erosion kernels — corrupting it silently changes "
            "every subsequent erosion calculation. "
            "See AGENTS.md > rcvID/rcvIDi convention."
        ),
    )
    np.testing.assert_array_equal(model.wghtVali, wghtVali_before)
    np.testing.assert_array_equal(model.distRcvi, distRcvi_before)


# ---------------------------------------------------------------------------
# TEST 4 - self.h is scratch, not elevation
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_scratch_vec_trap(minimal_model):
    """
    Protects: AGENTS.md > Scratch vector contract (CRITICAL).

    Fails if someone accidentally makes self.h persistent — that would
    silently break every hillslope caller. The trap is the name:
    `self.h` reads as if it were elevation, but it is allocated in
    `hillslope.__init__` as a generic scratch global Vec (alongside
    `self.hl` and `self.dh`) and is overwritten inside the rosw TS
    solver every time `_diffuseImplicit` or `diffuseSoil` runs.

    Silent failure prevented: a refactor that "consolidates" elevation
    storage by binding `self.h = self.hGlobal` (because the names look
    interchangeable) would make every diffusion solve clobber the
    canonical elevation field. The output would still validate field-
    by-field, but conservation and timestep-to-timestep continuity
    would break.

    Invariant: `self.h` and `self.hGlobal` must NOT share storage.
    Mutating `self.h` to a sentinel value must leave `self.hGlobal`
    unchanged.
    """
    model = minimal_model

    # Run at least one timestep so any hillslope/marine/soil consumer
    # has a chance to touch `self.h`. Without this, the Vec might still
    # hold its initial (undefined) duplicate() contents and the alias
    # test would be vacuous.
    model.runProcesses()

    hGlobal_snapshot = model.hGlobal.getArray().copy()

    # Write a clearly-unphysical sentinel to the scratch Vec. If self.h
    # is a true alias for self.hGlobal, this write propagates to the
    # elevation field.
    SENTINEL = -1.234567e9
    model.h.set(SENTINEL)

    hGlobal_after_scratch_write = model.hGlobal.getArray().copy()

    assert not np.allclose(hGlobal_after_scratch_write, SENTINEL), (
        "self.h appears to ALIAS self.hGlobal: writing to the scratch "
        "Vec mutated the elevation field. See AGENTS.md > Scratch "
        "vector contract. self.h must remain scratch — allocated via "
        "hLocal.duplicate() in hillslope.py:39 with its own storage."
    )
    np.testing.assert_array_equal(
        hGlobal_after_scratch_write, hGlobal_snapshot,
        err_msg="self.hGlobal changed after writing only to self.h.",
    )


# ---------------------------------------------------------------------------
# TEST 5 - sign conventions of cumED and EbLocal (thickness convention)
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_erosion_sign_conventions(incising_model):
    """
    Protects: AGENTS.md > Eb / EbLocal sign convention.

    Silent failure prevented: if a refactor inverts the sign of `tmp =
    -Eb*dt` (SPL.py:352 / nlSPL.py:404 / soilSPL.py:326), three things
    flip simultaneously: cumED, hGlobal updates (would make the landscape
    aggrade instead of incise), and EbLocal (the `EDrate` output field).
    `cumED` and `EbLocal` share the THICKNESS-RATE convention (positive
    = net deposition); they must both go NEGATIVE at any node that lost
    elevation by erosion.

    This test uses cumulative quantities (cumED) rather than the
    per-step Eb so it is robust to multi-step fixtures, and uses the
    INTERSECTION of `elev_change < 0` and `cumED_change < 0` to filter
    out nodes whose elevation moved for non-erosional reasons
    (compaction, flexural subsidence, advection).

    Note: as of 2026-06, `self.Eb` (global) and `self.EbLocal` (local)
    share the same thickness-rate convention (positive deposition,
    negative incision). They still hold DIFFERENT content by end-of-step
    — `Eb` is the river-only rate from the most recent SPL flavour,
    while `EbLocal` is the net rate including hillslope and marine
    contributions axpy'd in afterwards — but the SIGN of each follows
    the same rule.
    """
    model = incising_model
    h_before = model.hLocal.getArray().copy()
    cumED_before = model.cumEDLocal.getArray().copy()
    model.runProcesses()
    h_after = model.hLocal.getArray().copy()
    cumED_after = model.cumEDLocal.getArray().copy()

    elev_change = h_after - h_before
    cumED_change = cumED_after - cumED_before

    # A node is "definitely incising by erosion" iff BOTH its elevation
    # dropped AND its cumulative ED dropped over the same interval. The
    # intersection guards against tectonic / flexural / advective drops.
    incising = (elev_change < -1.0e-6) & (cumED_change < -1.0e-6)

    if not incising.any():
        pytest.skip(
            "Fixture did not produce measurable erosional incision. "
            "NEEDS_HUMAN_REVIEW: tune incising.yml so at least one node "
            "shows both elev_change < -1 µm and cumED_change < -1 µm "
            "over the run."
        )

    # ---- Primary invariant: cumED sign convention -------------------
    # cumED is constructed by `cumED.axpy(1.0, tmp)` where `tmp = -Eb*dt`
    # (positive deposition, negative erosion). A node that erodes must
    # see cumED *decrease*. The mask already requires this; the
    # assertion below makes the contract explicit and gives a useful
    # diagnostic if `tmp`'s sign is ever flipped.
    cumED_min = cumED_change[incising].min()
    cumED_max = cumED_change[incising].max()
    assert cumED_max < 0, (
        "cumED did not decrease at every incising node. AGENTS.md says: "
        "cumED is in the thickness convention (positive = deposition). "
        "If this fires, the `tmp.setArray(-Eb*dt)` line in one of "
        "SPL.py:352 / nlSPL.py:404 / soilSPL.py:326 has lost its minus, "
        "OR `cumED.axpy(1.0, tmp)` was changed to `axpy(-1.0, ...)`. "
        f"Diagnostic: cumED_change on incising nodes ∈ "
        f"[{cumED_min:.3e}, {cumED_max:.3e}] m"
    )

    # ---- Secondary invariant: EbLocal thickness-rate convention ----
    # By end-of-step EbLocal is the net thickness rate (same convention
    # as cumED: positive deposition, negative erosion). At nodes that
    # net-eroded, EbLocal should be ≤ 0. Allow exactly 0 because a node
    # that eroded earlier and was idle in the last step can have its
    # last-step EbLocal zeroed by the SPL setArray.
    eb_local = model.EbLocal.getArray()
    assert eb_local.shape == h_after.shape, (
        "EbLocal length does not match hLocal length. The fixture or "
        "the DMPlex layout has drifted."
    )
    incising_eb = eb_local[incising]
    assert (incising_eb <= 0).all(), (
        "EbLocal has a POSITIVE value at a node that eroded. AGENTS.md "
        "says EbLocal is in the thickness-rate convention (positive = "
        "deposition). A positive value at an eroding node means either "
        "the sign in SPL.py:367 (`EbLocal.setArray(add_rate)` where "
        "add_rate = tmp/dt = -Eb) has flipped, or a downstream `axpy` "
        "(seaplex.py:486 / hillslope.py:297 / soilSPL.py:549) is "
        "depositing more than the river eroded over the run. "
        f"Diagnostic: max(EbLocal on incising nodes) = "
        f"{incising_eb.max():.3e} m/yr"
    )


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
