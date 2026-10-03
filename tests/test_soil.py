"""
Soil / regolith production, gating and hillslope coupling.

Protects: docs/DESIGN_SOIL_REGOLITH.md.

Split out of the former tests/test_regression.py. Run this group alone with
`pytest -m soil`; see tests/README.md for the marker list.
"""

from __future__ import annotations

import numpy as np
import pytest

from _helpers import FIXTURES_DIR

# Skip the whole module (rather than erroring at collection) when the
# goSPL runtime stack is not installed.
pytest.importorskip(
    "gospl.tools.inputparser",
    reason="goSPL runtime deps (petsc4py / ruamel.yaml / scipy) not installed",
)

pytestmark = [pytest.mark.soil]


def test_soil_subaerial_gate(minimal_ice_soil_model):
    """
    Protects: `Lsoil` is a **subaerial** regolith cover. `soilSPL._subaqueousMask`
    flags marine nodes (`seaID`) AND ponded continental lakes (`pitIDs > -1` with
    `lFill > hl`), and every `Lsoil` write-back (the soil solve + `updateSoilThickness`)
    holds soil at 0 there. This is the coherent subaerial gate that replaced the
    former marine-only add-then-zero (marine deposition added soil that only the
    next fluvial solve wiped, at `seaID` only — continental lakes kept a spurious
    cover). See `DESIGN_SOIL_REGOLITH.md` §3/§5.

    Three-part guard: (a) end-to-end, no soil survives on marine (`seaID`) nodes
    after a full run; (b) a deterministic unit check that a *ponded continental*
    cell is flagged subaqueous (independent of whether the fixture has natural
    lakes) while a non-pit / non-ponded cell is not; (c) ice-covered land is
    **frozen inert** — a soil increment does not change it, and it is preserved
    (not zeroed, unlike the subaqueous case).
    """
    from mpi4py import MPI

    m = minimal_ice_soil_model
    assert m.cptSoil
    m.runProcesses()

    soil = m.Lsoil.getArray()
    assert np.isfinite(soil).all() and (soil >= -1.0e-9).all()

    # (a) marine gate — seaID is stable across the step, so soil must be 0 there.
    assert (soil[m.seaID] == 0.0).all(), "soil left on submarine (seaID) nodes"
    assert MPI.COMM_WORLD.allreduce(len(m.seaID), op=MPI.SUM) > 0, (
        "no marine nodes on the fixture — marine gate not exercised"
    )

    # (b) lake gate — force a synthetic ponded continental cell (a non-marine node
    # made an in-pit cell whose fill/spill level sits above the bed) and confirm
    # the mask flags it; a non-pit node at the same elevation must NOT be flagged.
    hl = m.hLocal.getArray().copy()
    pit_save, fill_save = m.pitIDs.copy(), m.lFill.copy()
    try:
        land = np.ones(m.lpoints, dtype=bool)
        land[m.seaID] = False
        idx = np.where(land)[0]
        if len(idx) > 0:
            j = int(idx[0])
            m.pitIDs[j], m.lFill[j] = 0, hl[j] + 5.0   # ponded 5 m below spill
            assert m._subaqueousMask(hl)[j], (
                "ponded continental lake node not flagged subaqueous"
            )
            m.pitIDs[j] = -1                            # not in a pit anymore
            assert not m._subaqueousMask(hl)[j], (
                "non-pit continental node wrongly flagged subaqueous"
            )
    finally:
        m.pitIDs[:], m.lFill[:] = pit_save, fill_save

    # (c) ice-freeze gate — cover one land cell with ice, apply a uniform +0.5 m
    # soil increment, and confirm the ice cell is preserved unchanged (frozen
    # inert) while a non-ice land cell takes the increment. Distinguishes freeze
    # (preserve) from the subaqueous zero.
    assert getattr(m, "iceOn", False), "fixture is not ice-enabled"
    hl = m.hLocal.getArray()
    land_idx = np.where(~m._subaqueousMask(hl))[0]
    assert len(land_idx) >= 2, "need >=2 subaerial land cells to exercise (c)"
    L0 = m.Lsoil.getArray().copy()
    ice_save = m.iceHL.getArray().copy()
    try:
        j = int(land_idx[0])                            # will be ice-covered
        L0[j] = max(L0[j], 0.3)                         # ensure a nonzero column to preserve
        m.Lsoil.setArray(L0)
        m.iceHL.setArray(np.where(np.arange(m.lpoints) == j, 1.0, 0.0))
        m.tmp.set(0.5)                                  # +0.5 m soil increment everywhere
        m.updateSoilThickness()
        L1 = m.Lsoil.getArray()
        assert np.isclose(L1[j], L0[j]), "soil under ice not frozen (changed)"
        assert L1[j] > 0.0, "frozen ice soil was zeroed (should be preserved)"
        k = int(land_idx[1])                            # non-ice land cell (if not subaqueous)
        if not m._subaqueousMask(hl)[k]:
            assert L1[k] >= L0[k], "non-ice land soil should not lose the increment"
    finally:
        m.iceHL.setArray(ice_save)
        m.Lsoil.setArray(L0)
        m.dm.localToGlobal(m.Lsoil, m.Gsoil)


def test_soil_mode_regolith():
    """
    Protects: soil `mode: regolith` (DESIGN_SOIL_REGOLITH.md Option 2.5, step 2)
    on the intended config — soil + stratigraphy (`time: strat:` → `stratNb>0`).

    In regolith mode `Lsoil` is the WEATHERING-produced regolith only; deposited
    sediment (fluvial transport-limited, lake/pit, marine) is NOT routed into
    `Lsoil` — it lives in the stratigraphy and carries its own SOFT erodibility
    there: a freshly deposited layer gets `stratK = Ksoil/K` (so the SPL bedrock
    term `Kbr·stratK = Ksoil`, i.e. fresh sediment erodes like soil). In lumped
    mode the deposit becomes soil instead, so its `stratK` stays `1.0`.

    Default `mode: lumped` is byte-identical (`regolithSoil=False` → every gate
    takes the legacy branch; the full suite staying green confirms it). Regolith
    mode changes soil BOOKKEEPING (and the *next* step's deposit erodibility), not
    the current step's elevation SOLVE — the depositional growth is removed
    post-solve, the deposition→soil calls are skipped, and the fresh `stratK` is
    written after the erosion solve. So one step from an identical initial state
    gives **identical elevations**, while `Lsoil ≤` the lumped `Lsoil`.
    """
    import os
    from gospl.model import Model

    fx = str(FIXTURES_DIR)
    if not os.path.exists(os.path.join(fx, "minimal_soil_strata.yml")):
        pytest.skip("minimal_soil_strata.yml fixture not present")

    cwd = os.getcwd()
    os.chdir(fx)
    try:
        A = Model("minimal_soil_strata.yml", verbose=False, showlog=False)  # lumped
        B = Model("minimal_soil_strata.yml", verbose=False, showlog=False)  # regolith
    finally:
        os.chdir(cwd)

    try:
        assert A.cptSoil and A.stratNb > 0, "fixture must be soil + stratigraphy"
        assert A.regolithSoil is False, "default soil mode must be lumped"
        B.regolithSoil = True                     # regolith mode

        A.tEnd = A.tNow + 0.5 * A.dt              # exactly one step each
        B.tEnd = B.tNow + 0.5 * B.dt
        A.runProcesses()
        B.runProcesses()

        La, Lb = A.Lsoil.getArray(), B.Lsoil.getArray()
        ha, hb = A.hGlobal.getArray(), B.hGlobal.getArray()

        # (1) regolith mode does not perturb the FIRST-step elevation solve.
        assert np.allclose(ha, hb, rtol=1.0e-9, atol=1.0e-6), (
            "regolith mode changed the elevation solve (should be bookkeeping-only)"
        )
        # (2) regolith soil never exceeds lumped soil (deposition not added).
        assert (Lb <= La + 1.0e-9).all(), "regolith Lsoil exceeds lumped somewhere"
        # (3) finite, non-negative, subaerial gate intact.
        assert np.isfinite(Lb).all() and (Lb >= -1.0e-9).all()
        assert (Lb[B.seaID] == 0.0).all()
        # (4) fresh deposits carry the soft erodibility Ksoil/K in regolith mode,
        #     and stay at 1.0 (never Ksoil/K) in lumped mode.
        ratio = B.Ksoil / B.K
        assert not np.isclose(ratio, 1.0), "fixture must have Ksoil != K to distinguish"
        assert np.isclose(B.stratK, ratio).any(), (
            "no freshly deposited layer carries the soft stratK = Ksoil/K"
        )
        assert not np.isclose(A.stratK, ratio).any(), (
            "lumped-mode deposit stratK should be 1.0, never Ksoil/K"
        )
    finally:
        A.destroy()
        B.destroy()


def test_soil_hillslope_conservation(minimal_ice_soil_model):
    """
    Protects (soil step 4 — hillslope↔soil coupling audit, DESIGN_SOIL_REGOLITH.md
    §5): in a soil run the hillslope step delegates to soil creep
    (`getHillslope` → `diffuseSoil`, so there is no double-count with the plain
    `_hillSlope`/`_hillSlopeNL` path), and that creep is a well-behaved regolith
    transport:

      (1) **conserves volume** — divergence-form diffusion, so the net elevation
          change integrates to ~0 on a closed sphere (no borders);
      (2) **moves the soil with the surface** — ΔLsoil == Δelevation (the creep
          increment is added to both `Lsoil` and `hGlobal`);
      (3) **leaves the bedrock unchanged** — `lHbed = h − Lsoil` is preserved
          (creep moves soil, not rock).

    Guards against a future regression that decouples soil from elevation, moves
    bedrock by creep, or reinstates a double hillslope+soil diffusion.
    """
    m = minimal_ice_soil_model
    assert m.cptSoil
    assert len(m.idBorders) == 0, "expected a closed sphere (no boundary flux)"

    m.tEnd = m.tNow + 0.5 * m.dt
    m.runProcesses()                          # populate realistic soil / elevation

    area = m.larea
    L0 = m.Lsoil.getArray().copy()
    h0 = m.hLocal.getArray().copy()
    m.diffuseSoil()                           # isolated soil-creep step
    dh = m.hLocal.getArray() - h0             # elevation change from creep
    dL = m.Lsoil.getArray() - L0              # soil change
    dbed = (m.hLocal.getArray() - m.Lsoil.getArray()) - (h0 - L0)  # lHbed change

    activity = float(np.sum(np.abs(dh) * area))
    assert activity > 0.0, "soil creep did nothing — test not exercised"
    # (1) volume-conserving creep on the closed sphere (net ≈ 0).
    net = float(np.sum(dh * area))
    assert abs(net) / activity < 1.0e-3, "soil creep is not volume-conserving"
    # (2) soil follows the surface.
    assert np.max(np.abs(dL - dh)) < 1.0e-6, "soil change != elevation change"
    # (3) bedrock unchanged by creep.
    assert np.max(np.abs(dbed)) < 1.0e-6, "soil creep moved the bedrock (lHbed changed)"
