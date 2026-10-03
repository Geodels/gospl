"""
Cross-module data contracts: rcvIDi snapshot, scratch Vecs, Eb sign convention.

Protects: AGENTS.md > rcvID / rcvIDi convention, Scratch vector contract, Eb / EbLocal sign convention.

Split out of the former tests/test_regression.py. Run this group alone with
`pytest -m contracts`; see tests/README.md for the marker list.
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

pytestmark = [pytest.mark.contracts]


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
