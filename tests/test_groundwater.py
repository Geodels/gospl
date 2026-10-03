"""
Water table, recharge, baseflow and duricrust.

Protects: docs/DESIGN_WATERTABLE_DURICRUST.md.

Split out of the former tests/test_regression.py. Run this group alone with
`pytest -m groundwater`; see tests/README.md for the marker list.
"""

from __future__ import annotations

import numpy as np
import pytest

from _helpers import FIXTURES_DIR, _gw_model

# Skip the whole module (rather than erroring at collection) when the
# goSPL runtime stack is not installed.
pytest.importorskip(
    "gospl.tools.inputparser",
    reason="goSPL runtime deps (petsc4py / ruamel.yaml / scipy) not installed",
)

pytestmark = [pytest.mark.groundwater]


def test_groundwater_opt_in(minimal_model):
    """
    Protects (water-table + duricrust, **Phase 0** — DESIGN_WATERTABLE_DURICRUST.md):
    the feature is OPT-IN and inert when off. Without a ``groundwater:`` block,
    ``gwOn`` is False and no groundwater/duricrust state is allocated (the full
    suite staying green confirms the default path is byte-identical). With the
    block, ``gwOn`` is True, the persistent state Vecs exist, the head is seeded
    to the surface, and a run completes — Phase 0 allocates state but does not yet
    solve, so it remains byte-identical.
    """
    import os
    from gospl.model import Model

    # (a) default OFF — no groundwater block ⇒ inert, nothing allocated.
    m = minimal_model
    assert m.gwOn is False, "groundwater must default off"
    assert not hasattr(m, "headG"), "no groundwater state should exist when off"

    # (b) opt-in ON.
    fx = str(FIXTURES_DIR)
    if not os.path.exists(os.path.join(fx, "minimal_gw.yml")):
        pytest.skip("minimal_gw.yml fixture not present")
    cwd = os.getcwd()
    os.chdir(fx)
    try:
        B = Model("minimal_gw.yml", verbose=False, showlog=False)
    finally:
        os.chdir(cwd)

    try:
        assert B.gwOn is True and B.duriOn is True
        for attr in ("headG", "headL", "duriHG", "duriHL", "rechargeL", "baseflowL"):
            assert hasattr(B, attr), f"missing groundwater state: {attr}"
        assert B.duriF.shape == (B.lpoints,)
        # head seeded at the aquifer BASE (a dry start, z_bed = z − aquifer_base),
        # so recharge fills it UP to the steady water table (seeding at the
        # surface would let the seepage clip pin it there — spurious saturation).
        assert np.allclose(
            B.headG.getArray(), B.hGlobal.getArray() - B.gwAquiferBase
        )
        B.runProcesses()          # Phase 0: state allocated but unused ⇒ completes
    finally:
        B.destroy()


def test_groundwater_recharge():
    """
    Protects (water-table + duricrust, **Phase 1** — DESIGN_WATERTABLE_DURICRUST.md
    §3): net recharge `R = f_infil · max(0, rain − evap)` (m/yr), held at 0 where
    the surface is not subaerial land — marine `seaID`, ponded continental lake,
    **or ice-covered** (rain doesn't infiltrate under ice). `f_infil` may be a
    scalar OR a per-vertex field. Verified by calling `updateGroundwater` directly
    so the mask matches the computation exactly.
    """
    import os
    from gospl.model import Model
    from gospl.tools.constants import ICE_COVER_MIN

    fx = str(FIXTURES_DIR)
    if not os.path.exists(os.path.join(fx, "minimal_gw.yml")):
        pytest.skip("minimal_gw.yml fixture not present")
    cwd = os.getcwd()
    os.chdir(fx)
    try:
        m = Model("minimal_gw.yml", verbose=False, showlog=False)
    finally:
        os.chdir(cwd)

    try:
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()                      # populate rain / seaID / drainage
        hl = m.hLocal.getArray()
        rain = m.rainVal

        def water_sub():
            """marine + ponded-lake mask on the current state (no ice)."""
            s = np.zeros(m.lpoints, dtype=bool)
            s[m.seaID] = True
            if getattr(m, "pitIDs", None) is not None and getattr(m, "lFill", None) is not None:
                s |= (m.pitIDs > -1) & (m.lFill > hl)
            return s

        # (humid, scalar f) R = f·max(0, rain−evap) with the under-water zeroing.
        m.updateGroundwater()
        R = m.rechargeL.getArray()
        exp = np.where(water_sub(), 0.0, m.gwInfiltration * np.maximum(0.0, rain))
        assert np.allclose(R, exp), "recharge != f·max(0, rain−evap) (gated)"
        assert (R[m.seaID] == 0.0).all(), "recharge left on submarine nodes"
        assert (R >= 0.0).all() and np.isfinite(R).all()
        assert R.max() > 0.0, "no recharge anywhere — test not exercised"

        # (spatial infiltration) a per-vertex f field must be honoured node-wise.
        f0 = m.gwInfiltration
        f_arr = np.linspace(0.1, 0.6, m.lpoints)
        m.gwInfiltration = f_arr
        m.updateGroundwater()
        Rs = m.rechargeL.getArray()
        exps = np.where(water_sub(), 0.0, f_arr * np.maximum(0.0, rain))
        assert np.allclose(Rs, exps), "per-vertex infiltration field not honoured"
        m.gwInfiltration = f0

        # (ice gate) ice-covered land gets no rain-recharge. Synthesise an ice
        # cover on a few subaerial cells (minimal_gw has no ice block).
        land = np.where(~water_sub())[0]
        assert len(land) >= 3, "need subaerial land to exercise the ice gate"
        icecells = land[:3]
        m.iceOn = True
        m.iceHL = m.hLocal.duplicate()
        icearr = np.zeros(m.lpoints, dtype=np.float64)
        icearr[icecells] = 10.0 * ICE_COVER_MIN
        m.iceHL.setArray(icearr)
        try:
            m.updateGroundwater()
            assert (m.rechargeL.getArray()[icecells] == 0.0).all(), (
                "ice-covered land received rain-recharge"
            )
        finally:
            m.iceHL.destroy()
            m.iceOn = False

        # (arid) evap > rain everywhere ⇒ net < 0 ⇒ no recharge at all.
        m.evapVal = np.full(m.lpoints, float(rain.max()) + 1.0)
        m.updateGroundwater()
        assert (m.rechargeL.getArray() == 0.0).all(), "recharge should be 0 when evap > rain"
    finally:
        m.destroy()


def test_groundwater_recharge_refinements():
    """
    Protects the opt-in recharge refinements (DESIGN_WATERTABLE_DURICRUST.md §3;
    all default off ⇒ unchanged): **slope-modulated** infiltration
    (`f/(1+slope/infil_slope_ref)` — less on steep terrain), **subglacial-meltwater
    recharge** (a fraction of `iceMeltRiverL` infiltrates under ice, the one
    recharge path the ice gate allows), and **lithology-modulated** infiltration
    (`fine_infil_factor` — coarse infiltrates more than fine, dual lithology).
    """
    from gospl.tools.constants import ICE_COVER_MIN

    m = _gw_model("minimal_gw.yml")
    try:
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()                       # builds rcvID / rainVal
        own = m.inIDs == 1

        m.updateGroundwater()
        R0 = m.rechargeL.getArray().copy()

        # (slope) infiltration falls on slopes — never increases, strictly less
        # where the terrain is steep and there was recharge.
        m.gwInfilSlopeRef = 0.01
        m.updateGroundwater()
        Rs = m.rechargeL.getArray().copy()
        m.gwInfilSlopeRef = 0.0
        slope = m._surfaceSlope()
        steep = own & (slope > 0.0) & (R0 > 0.0)
        assert steep.any(), "no sloped recharge cells — test not exercised"
        assert (Rs[steep] < R0[steep]).all(), "slope did not reduce infiltration"
        assert (Rs <= R0 + 1.0e-12).all()

        # (subglacial) synthesise ice + meltwater on a few subaerial land cells:
        # the rain path is gated off under ice, subglacial adds frac·imr/area.
        sub = np.zeros(m.lpoints, dtype=bool)
        sub[m.seaID] = True
        if getattr(m, "pitIDs", None) is not None:
            sub |= (m.pitIDs > -1) & (m.lFill > m.hLocal.getArray())
        land = np.where(~sub & own)[0][:5]
        assert len(land) >= 3
        m.iceOn = True
        m.iceHL = m.hLocal.duplicate()
        ice = np.zeros(m.lpoints)
        ice[land] = 10.0 * ICE_COVER_MIN
        m.iceHL.setArray(ice)
        m.iceMeltRiverL = m.hLocal.duplicate()
        imr = np.zeros(m.lpoints)
        imr[land] = 1.0e6                       # m³/yr glacial meltwater
        m.iceMeltRiverL.setArray(imr)
        m.gwSubglacial = 0.5
        try:
            m.updateGroundwater()
            Rg = m.rechargeL.getArray()
            assert (Rg[land] > 0.0).all(), "no subglacial recharge under ice"
            assert np.allclose(Rg[land], 0.5 * imr[land] / m.larea[land])
        finally:
            m.iceHL.destroy()
            m.iceMeltRiverL.destroy()
            m.iceOn = False
            m.gwSubglacial = 0.0
    finally:
        m.destroy()

    # (lithology) coarse infiltrates more than fine — needs dual lithology.
    mc = _gw_model("minimal_gw_combo.yml")
    try:
        mc.tEnd = mc.tNow + 0.5 * mc.dt
        mc.runProcesses()
        mc.updateGroundwater()
        Rc0 = mc.rechargeL.getArray().copy()
        mc.gwFineInfilFactor = 0.2
        mc.updateGroundwater()
        Rcf = mc.rechargeL.getArray().copy()
        own = mc.inIDs == 1
        fine = own & (mc._surfaceComposition() < 1.0) & (Rc0 > 0.0)
        assert fine.any(), "no fine-bearing recharge cells — test not exercised"
        assert (Rcf[fine] < Rc0[fine]).all(), "fine did not reduce infiltration"
    finally:
        mc.destroy()


def test_watertable_solve():
    """
    Protects (water-table + duricrust, **Phase 2** — DESIGN_WATERTABLE_DURICRUST.md
    §2/§3): the implicit Dupuit–Boussinesq head solve `(I + (Δt/S)·L(T))h = h_old +
    (Δt/S)·R` (Picard on `T(h)`, seepage clip). The solved head must sit between
    the two free boundaries — the impermeable base `z_bed = z − aquifer_base` and
    the seepage surface `z` — be finite, and (on a permeable aquifer) produce a
    NON-trivial water table below the surface (so the interior elliptic solve is
    actually exercised, not fully Dirichlet-pinned).

    Called directly after a step so head / `wtDepth` are consistent with the
    current surface (erosion later in the step lowers `z`). Analytic Dupuit-profile
    validation + np=1-vs-2 invariance are a follow-up increment.
    """
    import os
    from gospl.model import Model

    fx = str(FIXTURES_DIR)
    if not os.path.exists(os.path.join(fx, "minimal_gw.yml")):
        pytest.skip("minimal_gw.yml fixture not present")
    cwd = os.getcwd()
    os.chdir(fx)
    try:
        m = Model("minimal_gw.yml", verbose=False, showlog=False)
    finally:
        os.chdir(cwd)

    try:
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()
        m.updateGroundwater()                 # solve on the current surface

        z = m.hLocal.getArray()
        h = m.headL.getArray()
        wt = m.wtDepth
        zbed = z - float(m.gwAquiferBase)
        own = m.inIDs == 1

        # Bounded head between the two free boundaries.
        assert np.isfinite(h).all(), "non-finite head"
        assert (h <= z + 1.0e-6).all(), "head above the surface (seepage clip failed)"
        assert (h >= zbed - 1.0e-6).all(), "head below the aquifer base"
        # wtDepth = z − h, within [0, aquifer_base].
        assert np.allclose(wt, z - h), "wtDepth != z − h"
        assert (wt >= -1.0e-6).all() and (wt <= m.gwAquiferBase + 1.0e-6).all()
        # Non-trivial water table (interior solve exercised, not all pinned).
        assert (wt[own] > 0.01).any(), "water table fully saturated — solve not exercised"
    finally:
        m.destroy()


def test_watertable_steady():
    """
    Protects (water-table + duricrust, **Phase 2 analytic-ish validation**): with
    the surface and recharge held fixed, repeated implicit head solves converge to
    a **steady water table** — the quasi-steady Dupuit fixed point `∇·(T∇h)+R=0`
    (each `updateGroundwater` uses the previous head as `h_old`, so iterating is
    an outer relaxation onto steady state). The head change between successive
    solves must collapse toward zero.

    A full analytic Dupuit-parabola benchmark (flat hillslope fixture + drain,
    compared to `h² = h_L² + (R/K)·x·(L−x)`) is a `benchmarks/`-suite follow-up.
    """
    import os
    from gospl.model import Model

    fx = str(FIXTURES_DIR)
    if not os.path.exists(os.path.join(fx, "minimal_gw.yml")):
        pytest.skip("minimal_gw.yml fixture not present")
    cwd = os.getcwd()
    os.chdir(fx)
    try:
        m = Model("minimal_gw.yml", verbose=False, showlog=False)
    finally:
        os.chdir(cwd)

    try:
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()                      # fix a surface / recharge

        deltas = []
        prev = None
        for _ in range(55):
            m.updateGroundwater()             # outer relaxation onto steady state
            h = m.headL.getArray().copy()
            if prev is not None:
                deltas.append(float(np.max(np.abs(h - prev))))
            prev = h

        deltas = np.array(deltas)
        assert np.isfinite(deltas).all() and (deltas >= 0).all()
        # Contracting fixed-point iteration: the head-change per solve collapses
        # toward zero. Geometric decay (ratio ~0.96/iter on this coarse pitted
        # mesh — full mm-level convergence takes ~60 iters), so we assert a clear
        # multi-fold contraction over the relaxation, not full convergence.
        assert deltas[-1] < 0.3 * deltas[0], (
            f"head not converging to steady state: first Δ={deltas[0]:.3g}, "
            f"last Δ={deltas[-1]:.3g}"
        )
        assert deltas[-5:].mean() < deltas[:5].mean(), "no convergence trend"
        # The water-table depth is bounded by the aquifer thickness (z - z_bed):
        # a water table cannot sit below its aquifer base (the head can drain
        # below it via the min_sat_thickness floor, but the REPORTED wtDepth is
        # clamped). So 0 <= wtDepth <= z - z_bed everywhere.
        z = m.hLocal.getArray()
        aq = np.maximum(z - m._gwZbed(z), 0.0)
        assert (m.wtDepth >= -1.0e-9).all(), "negative water-table depth"
        assert (m.wtDepth <= aq + 1.0e-6).all(), "wtDepth exceeds the aquifer thickness"
    finally:
        m.destroy()


def test_duricrust_forms_at_fringe():
    """
    Protects (water-table + duricrust, **Phase 3** — DESIGN_WATERTABLE_DURICRUST.md
    §3 step 5): the capillary-fringe duricrust forms **selectively at the fringe**.
    With a synthetic water-table depth, `_updateDuricrust` must grow `duriH` where
    `wt ≈ fringe_depth` (Gaussian favourability `Φ→1`) and leave it ~0 where the
    table is far from the fringe (`Φ→0`); the induration `duriF = duriH/duriH_max`
    stays in [0,1] and drives `duriKarmor = 1 − armor_max·duriF`.
    """
    m = _gw_model("minimal_gw.yml")
    try:
        assert m.duriOn, "expected duricrust on"
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()                       # populate rain / drainage

        z = m.hLocal.getArray()
        m.gwZlast = z.copy()                   # zero the incision term
        n = m.lpoints
        # Crust forms only on SUBAERIAL land (no capillary fringe under the sea /
        # a lake), so restrict the fringe/far node sets to subaerial nodes.
        sa = np.where(m._subaerialMask())[0]
        assert sa.size > 20, "fixture has too few subaerial nodes to test"
        # Synthetic water table: most nodes far from the fringe (Φ≈0); a subset
        # sitting exactly at the fringe depth (Φ=1).
        m.wtDepth = np.full(n, 50.0)
        fringe = sa[::7]                       # subaerial, at the fringe depth
        m.wtDepth[fringe] = m.duriFringeDepth
        far = np.setdiff1d(sa, fringe)         # subaerial, away from the fringe
        m.duriHL.set(0.0)

        for _ in range(200):
            m._updateDuricrust()

        duriH = m.duriHL.getArray()
        assert np.isfinite(duriH).all()
        assert (duriH >= 0.0).all() and (duriH <= m.duriMaxThick + 1e-9).all()
        # Crust grows at the fringe, essentially none far from it.
        assert duriH[fringe].min() > 0.0, "no crust formed at the fringe"
        assert duriH[far].max() < 0.01 * duriH[fringe].mean(), (
            "crust formed away from the capillary fringe"
        )
        # Induration + armor multiplier are consistent and in range.
        assert np.allclose(m.duriF, duriH / m.duriMaxThick)
        assert (m.duriF >= 0.0).all() and (m.duriF <= 1.0).all()
        assert np.allclose(m.duriKarmor, 1.0 - m.duriArmorMax * m.duriF)
        assert (m.duriKarmor <= 1.0).all() and (m.duriKarmor > 0.0).all()
    finally:
        m.destroy()


def test_duricrust_discharge_gate():
    """
    Protects the absolute-accumulation gate (DESIGN_WATERTABLE_DURICRUST.md §3b):
    `_dischargeWeight` restricts crust formation to groundwater DISCHARGE zones.

    - Gate OFF (default) → weight ≡ 1 everywhere (backwards-compatible: the crust
      forms wherever the fringe favourability `Φ` is non-zero).
    - Gate ON → `G = (−div q)⁺ / ((−div q)⁺ + R)` ∈ [0,1]: 0 at a recharge
      (divergent, `div q > 0`) node, →1 at a strongly convergent discharge node.
    """
    m = _gw_model("minimal_gw.yml")
    try:
        n = m.lpoints
        R = np.full(n, 0.2)                      # uniform recharge (m/yr)
        # Synthetic lateral divergence: recharge (divergent), neutral, and strong
        # convergence (discharge) populations.
        divq = np.zeros(n)
        divq[0::3] = +0.5                         # recharge / divergent  → G = 0
        divq[1::3] = 0.0                          # neutral               → G = 0
        divq[2::3] = -100.0                       # strong convergence    → G → 1

        # Gate OFF → all ones regardless of divq.
        m.duriDischargeGate = False
        assert np.allclose(m._dischargeWeight(divq, R), 1.0)

        # Gate ON → discharge-only weight in [0,1].
        m.duriDischargeGate = True
        G = m._dischargeWeight(divq, R)
        assert np.isfinite(G).all()
        assert (G >= 0.0).all() and (G <= 1.0).all()
        assert np.allclose(G[0::3], 0.0), "recharge (divergent) node must not form crust"
        assert np.allclose(G[1::3], 0.0), "neutral node must not form crust"
        assert (G[2::3] > 0.99).all(), "strong discharge node should form crust (G→1)"
    finally:
        m.destroy()


def test_duricrust_soilfree():
    """
    Protects (Phase 3): the duricrust ships **soil-independent** — it forms on a
    groundwater run with NO `soil:` block (`cptSoil` False). The default proxy
    weathering supply `Ψ = max(0, rain−evap)^p` needs no regolith, so the fringe
    crust must still grow.
    """
    m = _gw_model("minimal_gw_nosoil.yml")
    try:
        assert m.duriOn and not getattr(m, "cptSoil", False), (
            "fixture must be groundwater+duricrust with soil OFF"
        )
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()

        z = m.hLocal.getArray()
        m.gwZlast = z.copy()
        m.wtDepth = np.full(m.lpoints, m.duriFringeDepth)   # all at the fringe
        m.duriHL.set(0.0)
        for _ in range(200):
            m._updateDuricrust()

        duriH = m.duriHL.getArray()
        assert np.isfinite(duriH).all()
        assert duriH.max() > 0.0, "no duricrust formed on a soil-free run"
    finally:
        m.destroy()


def test_duricrust_weathering_rate():
    """
    Protects (Phase 3, §3a Level A): the explicit `weathering: mode: rate`
    (Maher–Chamberlain) supply `W = R·C_eq·(1 − exp(−Dw/(R·L)))·…` **responds to
    the recharge `R`** the head solve computes — a monotone, finite increase with
    `R` (the control the climate proxy lacks), and `W=0` where `R=0`.
    """
    m = _gw_model("minimal_gw.yml")
    try:
        # Switch the supply to the Level-A explicit rate.
        m.duriWeatherMode = "rate"
        m.duriWeatherCeq = 1.0
        m.duriWeatherDw = 1.0
        m.duriWeatherL = 20.0
        m.duriWeatherEa = 0.0
        m.duriWeatherability = 1.0
        m.prodSoil = None                      # disable the regolith cap for the test
        n = m.lpoints

        def W_for(rval):
            m.rechargeL.setArray(np.full(n, rval, dtype=float))
            return m._weatheringSupply()

        W0 = W_for(0.0)
        W_lo = W_for(0.05)
        W_hi = W_for(0.5)

        assert np.isfinite(W_lo).all() and np.isfinite(W_hi).all()
        assert (W0 == 0.0).all(), "W must be 0 where recharge is 0"
        assert (W_lo > 0.0).all(), "no weathering at positive recharge"
        assert (W_hi > W_lo).all(), "W must increase with recharge R"
    finally:
        m.destroy()


def test_duricrust_armors_K():
    """
    Protects (water-table + duricrust, **Phase 4** — DESIGN_WATERTABLE_DURICRUST.md
    §5): the single erodibility hook `_surfaceArmoringK` cuts `K` at an indurated
    cell so it erodes ≪ a bare one, composing multiplicatively in `_surfaceLithoK`
    with **no branching** in the eroders. And it is a byte-identical **no-op** when
    the duricrust is off (returns the scalar 1.0).
    """
    m = _gw_model("minimal_gw.yml")
    try:
        assert m.duriOn, "expected duricrust on"
        n = m.lpoints
        armor = float(m.duriArmorMax)

        # No-op when off: the armoring factor is the scalar 1.0 and _surfaceLithoK
        # is exactly the lithology-only multiplier.
        m.duriOn = False
        assert m._surfaceArmoringK() == 1.0
        base_litho = m._surfaceLithoK().copy()

        # Turn it on with a half-indurated pattern (duriF = 1 on a subset, 0 else).
        m.duriOn = True
        m.duriF = np.zeros(n, dtype=np.float64)
        hard = np.arange(0, n, 3)
        m.duriF[hard] = 1.0
        bare = np.setdiff1d(np.arange(n), hard)

        arm = m._surfaceArmoringK()
        assert isinstance(arm, np.ndarray)
        assert np.allclose(arm[hard], 1.0 - armor)     # fully indurated
        assert np.allclose(arm[bare], 1.0)             # bare unchanged

        # Effective erodibility K = surfaceK * surfaceLithoK: hard cells reduced
        # by exactly armor_max relative to their bare-litho value.
        K_on = m._surfaceK() * m._surfaceLithoK()
        K_bare_litho = m._surfaceK() * base_litho
        assert np.allclose(K_on[hard], (1.0 - armor) * K_bare_litho[hard])
        assert np.allclose(K_on[bare], K_bare_litho[bare])
        # Indurated cells are strictly less erodible than bare (where K>0).
        pos = K_bare_litho > 0
        assert (K_on[hard[pos[hard]]] < K_bare_litho[hard[pos[hard]]]).all()
    finally:
        m.destroy()


def test_groundwater_baseflow_conserves():
    """
    Protects (water-table + duricrust, **Phase 5** — DESIGN_WATERTABLE_DURICRUST.md
    §3 step 7): the opt-in baseflow closure conserves water — the seepage-return
    discharge `Σ baseflow` accounts to `Σ(R·A) − ΔS` and, relaxed to the
    quasi-steady water table (`ΔS→0`), matches the total recharge, so river
    discharge stays `≈ rain − evap`.
    """
    from mpi4py import MPI

    m = _gw_model("minimal_gw.yml")
    try:
        assert m.gwConserveBaseflow, "fixture must have conserve_baseflow on"
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()
        for _ in range(100):                   # relax onto steady state (ΔS→0)
            m.updateGroundwater()

        own = m.inIDs == 1
        A = m.larea
        Vrech = MPI.COMM_WORLD.allreduce(
            float((m.rechargeL.getArray() * A)[own].sum()), op=MPI.SUM
        )
        Vbase = MPI.COMM_WORLD.allreduce(
            float(m.baseflowL.getArray()[own].sum()), op=MPI.SUM
        )
        assert Vrech > 0.0, "no recharge — test not exercised"
        assert abs(Vbase - Vrech) < 0.02 * Vrech, (
            f"baseflow not conserving recharge: Σbase={Vbase:.4g} vs "
            f"Σrech={Vrech:.4g}"
        )
    finally:
        m.destroy()


def test_groundwater_baseflow_reinjection():
    """
    Protects (deferred item #4 — DESIGN_WATERTABLE_DURICRUST.md §3 step 7): with
    `conserve_baseflow`, the baseflow is **re-injected into the surface-flow
    source** (`applyForces` builds `bL = rain·A − recharge·A + baseflow`). The
    infiltrated recharge leaves surface runoff and returns at the seepage nodes,
    so the source is redistributed but **globally near-neutral** (Σ recharge ≈ Σ
    baseflow). Applied once per step (the single `bL` reset), so it does not
    double-count across the two per-step `flowAccumulation` calls.
    """
    m = _gw_model("minimal_gw.yml")
    try:
        assert m.gwConserveBaseflow, "fixture must have conserve_baseflow on"
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()
        for _ in range(40):                    # populate a steady recharge/baseflow
            m.updateGroundwater()

        own = m.inIDs == 1
        A = m.larea
        # Re-injection ON vs OFF via the single per-step bL reset.
        m.gwConserveBaseflow = True
        m.applyForces()
        bl_on = m.bL.getArray().copy()
        m.gwConserveBaseflow = False
        m.applyForces()
        bl_off = m.bL.getArray().copy()
        m.gwConserveBaseflow = True

        assert np.allclose(bl_off, m.rainVal * A), "off path != raw runoff source"
        assert np.isfinite(bl_on).all() and (bl_on >= 0.0).all()
        assert not np.allclose(bl_on[own], bl_off[own]), "re-injection not applied"
        # Globally near-neutral (recharge removed ≈ baseflow added at steady state).
        tot_on = float(bl_on[own].sum())
        tot_off = float(bl_off[own].sum())
        assert abs(tot_on - tot_off) < 0.05 * tot_off, (
            f"re-injection not water-neutral: on={tot_on:.4g} off={tot_off:.4g}"
        )
    finally:
        m.destroy()


def test_groundwater_lake_exchange():
    """
    Protects (deferred item #2 — DESIGN_WATERTABLE_DURICRUST.md §15): the opt-in
    lake ↔ aquifer VOLUME coupling. When `lake_exchange` is on, the signed
    across-bed groundwater flux (`gwLakeFlux`, `∇·(T∇h)·A = −(L·h)·A`) is computed
    each head solve — positive where the aquifer discharges into the surface
    (gaining), negative where it leaks (losing) — and fed into the per-lake fill
    budget in `_distributeDownstream`. It is a strict **no-op when off** (default):
    `gwLakeFlux` stays zero and the lake cascade is byte-identical.
    """
    # OFF (default): the flux is never computed — stays zero.
    m0 = _gw_model("minimal_gw.yml")
    try:
        assert not m0.gwLakeExchange, "lake_exchange must default off"
        m0.tEnd = m0.tNow + m0.dt
        m0.runProcesses()
        assert (m0.gwLakeFlux == 0.0).all(), "flux computed while lake_exchange off"
    finally:
        m0.destroy()

    # ON: signed across-bed flux computed and the cascade coupling runs.
    m = _gw_model("minimal_gw.yml")
    try:
        m.gwLakeExchange = True                # opt-in the volume coupling
        m.tEnd = m.tNow + 2 * m.dt
        m.runProcesses()                       # step 1: head solve populates flux
        m.runProcesses()                       # step 2: cascade consumes it
        f = m.gwLakeFlux
        assert np.isfinite(f).all()
        assert (f != 0.0).any(), "no across-bed groundwater flux computed"
        # A non-trivial water table both gains (discharge>0) and loses (leak<0).
        assert (f > 0.0).any() and (f < 0.0).any()
        assert np.isfinite(m.gwLakeInflow)     # diagnostic accumulated, run stable
    finally:
        m.destroy()


def test_groundwater_from_soil():
    """
    Protects (Phase 5, §8 soil coupling): `aquifer_base: from_soil` ties the
    aquifer floor to the bedrock elevation `z_bed = lHbed − bedrock_depth`
    (permeable regolith over impermeable bedrock). In a **depositional basin**
    the porous sediment fill IS the aquifer, so the base deepens to the bottom of
    the (non-sentinel) stratigraphic pile — the base is the deeper of
    `lHbed − bedrock_depth` and `z − Σ sediment`. The head stays bounded between
    that floor and the surface.
    """
    m = _gw_model("minimal_gw.yml")
    try:
        assert getattr(m, "cptSoil", False), "fixture must track soil (lHbed)"
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()                       # sets lHbed = z − Lsoil

        z = m.hLocal.getArray()
        m.gwAquiferBase = "from_soil"
        m.gwBedrockDepth = 2.0
        lo = int(getattr(m, "bedrockLay", 0))
        top = m.stratStep + 1

        # Base = deeper of (lHbed − bedrock_depth) and the sediment-pile bottom.
        zbed = m._gwZbed(z)
        exp = m.lHbed.getArray() - 2.0
        if top > lo:
            exp = np.minimum(exp, z - m.stratH[:, lo:top].sum(axis=1))
        assert np.allclose(zbed, exp), "from_soil z_bed formula mismatch"

        # Depositional basin: a thick porous fill drives the base to the pile
        # bottom (well below the thin-regolith lHbed).
        m.stratH[:, lo] = 40.0
        zbed_b = m._gwZbed(z)
        sed = m.stratH[:, lo : m.stratStep + 1].sum(axis=1)
        deep = (z - sed) < (m.lHbed.getArray() - 2.0)
        assert deep.any(), "basin fill did not deepen the base — test not exercised"
        assert np.allclose(zbed_b[deep], (z - sed)[deep]), (
            "basin aquifer base != bottom of the sediment pile"
        )

        m.updateGroundwater()                  # solve on the bedrock floor
        h = m.headL.getArray()
        assert np.isfinite(h).all()
        assert (h >= m._gwZbed(z) - 1.0e-6).all(), "head below the aquifer floor"
        assert (h <= z + 1.0e-6).all(), "head above the surface (seepage failed)"
    finally:
        m.destroy()


def test_duricrust_regolith_limited():
    """
    Protects (Phase 5, §8 soil coupling): when soil is tracked, duricrust
    formation is **regolith-supply-limited** — capped by `_regolithSupplyRate`
    (`prodSoil·rain`), so chemical crust growth cannot outpace physical regolith
    production. Where the cap binds, the crust grows strictly less than the
    unlimited (soil-off) formation.
    """
    m = _gw_model("minimal_gw.yml")
    try:
        assert m.cptSoil and m.duriOn
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()

        z = m.hLocal.getArray()
        m.wtDepth = np.full(m.lpoints, m.duriFringeDepth)   # Φ=1 everywhere
        supply_unlim = m.duriFormRate * 1.0 * 1.0           # Φ=Ψ=1 (proxy, rain=1)
        reg = m._regolithSupplyRate()
        binds = reg < supply_unlim
        assert binds.any(), "regolith cap never binds — test not exercised"

        # Regolith-limited formation (cptSoil on).
        m.duriHL.set(0.0)
        m.gwZlast = z.copy()
        m._updateDuricrust()
        dH_lim = m.duriHL.getArray().copy()

        # Unlimited formation (soil coupling off).
        saved = m.cptSoil
        m.cptSoil = False
        m.duriHL.set(0.0)
        m.gwZlast = z.copy()
        m._updateDuricrust()
        dH_unlim = m.duriHL.getArray().copy()
        m.cptSoil = saved

        # Formation is subaerial-only (Φ gated off the sea / lakes), so compare on
        # subaerial nodes; marine/ponded nodes are 0 in both.
        sa = m._subaerialMask()
        assert np.allclose(dH_lim[sa], m.dt * np.minimum(supply_unlim, reg)[sa])
        assert np.allclose(dH_unlim[sa], m.dt * supply_unlim)
        assert (dH_lim[binds & sa] < dH_unlim[binds & sa]).all(), "regolith cap not applied"
    finally:
        m.destroy()


def test_duricrust_strata_exhumation():
    """
    Protects (water-table + duricrust, **Phase 6** — DESIGN_WATERTABLE_DURICRUST.md
    §9): the per-layer induration archive `stratDuri` records the crust and a
    **buried crust re-armors on re-exposure**. `_recordInduration` syncs the live
    `duriF` with the top non-empty layer: while buried (fresh layers on top,
    `stratDuri=0`) the surface stays weak; once erosion exhumes the indurated
    layer, `duriF` (and `duriKarmor`) re-arm to the preserved value.
    """
    m = _gw_model("minimal_gw.yml")
    try:
        assert m.gwOn and m.stratNb > 0, "need groundwater + stratigraphy"
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()
        assert m.stratDuri is not None, "induration archive not allocated"

        n = m.lpoints
        # Build a 3-layer column everywhere: deep indurated crust (layer 0),
        # then two fresh (uncemented) layers on top.
        m.stratStep = 2
        m.stratH[:] = 0.0
        m.stratH[:, 0] = 5.0            # buried, indurated
        m.stratH[:, 1] = 3.0
        m.stratH[:, 2] = 2.0            # fresh surface
        m.stratDuri[:] = 0.0
        m.stratDuri[:, 0] = 1.0         # relict crust locked in the record
        # Live crust thickness > 0 so the write-down reaches the top layer at every
        # node (the write-down spans `duriH` below the surface; the gated run may
        # leave duriH=0 on marine nodes, which is not what this archive test probes).
        m.duriHL.set(1.0)

        # Buried: the exposed top layer is uncemented → surface stays weak.
        m.duriF = np.full(n, 0.1)
        m._recordInduration()
        assert (m.duriF < 0.5).all(), "buried crust wrongly armored the surface"
        # Write-down recorded the (weak) live crust into the exposed top layer.
        assert np.allclose(m.stratDuri[:, 2], 0.1)
        # The buried relict is preserved.
        assert np.allclose(m.stratDuri[:, 0], 1.0)

        # Exhume: erode the two upper layers so layer 0 becomes the surface.
        m.stratH[:, 1] = 0.0
        m.stratH[:, 2] = 0.0
        m.duriF = np.full(n, 0.1)
        m._recordInduration()
        assert np.allclose(m.duriF, 1.0), "exhumed crust did not re-arm the surface"
        assert np.allclose(m.duriKarmor, 1.0 - m.duriArmorMax)
    finally:
        m.destroy()


def test_duricrust_multilayer_record():
    """
    Protects (§9 formation depth range): a crust of thickness `duriH` is recorded
    into **every stratigraphic layer within `duriH` below the surface**, not just
    the top layer — so a thick crust spanning several thin layers is preserved
    over its full thickness (and re-arms the surface across that whole span on
    exhumation). Layers deeper than `duriH` are untouched.
    """
    m = _gw_model("minimal_gw.yml")
    try:
        assert m.gwOn and m.stratNb > 0
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()
        assert m.stratDuri is not None

        n = m.lpoints
        # 4-layer column; depth-to-top from the surface: L3=0, L2=1, L1=2, L0=3.
        m.stratStep = 3
        m.stratH[:] = 0.0
        m.stratH[:, 0] = 2.0
        m.stratH[:, 1] = 1.0
        m.stratH[:, 2] = 1.0
        m.stratH[:, 3] = 1.0            # surface layer
        m.stratDuri[:] = 0.0
        m.duriHL.setArray(np.full(n, 2.5))     # crust 2.5 m thick from the surface
        m.duriF = np.full(n, 0.7)

        m._recordInduration()
        d = m.stratDuri
        # Within 2.5 m of the surface: L3 (0), L2 (1), L1 (2) — all recorded.
        assert np.allclose(d[:, 3], 0.7)
        assert np.allclose(d[:, 2], 0.7)
        assert np.allclose(d[:, 1], 0.7), "crust not recorded across its full depth"
        # L0 top is 3 m down (> 2.5) — below the crust, untouched.
        assert np.allclose(d[:, 0], 0.0), "induration written below the crust depth"
    finally:
        m.destroy()


def test_groundwater_dual_provenance_combo():
    """
    Protects (water-table + duricrust, **§10 compatibility**): the feature runs
    with BOTH dual lithology AND in-model provenance on. The induration archive
    `stratDuri` is one more independent, intensive per-layer field alongside the
    fine pile (`stratHf`/`phiF`) and the provenance partition (`stratP`); an
    end-to-end run completes and every conservation invariant still holds —
    Σ-over-classes == `stratH`, fine ≤ total, induration ∈ [0,1].
    """
    m = _gw_model("minimal_gw_combo.yml")
    try:
        assert m.gwOn and m.duriOn and m.stratLith and m.provOn and m.stratNb > 0, (
            "combo fixture must enable gw+duricrust+dual+provenance+strata"
        )
        assert m.stratDuri is not None and m.stratHf is not None
        assert m.stratP is not None

        while m.tNow < m.tEnd:
            m.runProcesses()

        own = m.inIDs == 1
        top = m.stratStep + 1
        # Provenance partition: Σ over classes == layer thickness.
        psum = m.stratP[:, :top, :].sum(axis=2)
        assert np.allclose(psum[own], m.stratH[own, :top], atol=1.0e-6, rtol=1.0e-4), (
            "provenance Σ-over-classes drifted from stratH"
        )
        # Dual lithology: fine bulk never exceeds the layer total.
        assert (m.stratHf[:, :top] <= m.stratH[:, :top] + 1.0e-9).all()
        # Induration archive: finite, in [0,1], and armor multiplier consistent.
        assert np.isfinite(m.stratDuri).all()
        assert (m.stratDuri >= 0.0).all() and (m.stratDuri <= 1.0 + 1.0e-9).all()
        assert np.allclose(m.duriKarmor, 1.0 - m.duriArmorMax * m.duriF)
    finally:
        m.destroy()


def test_groundwater_restart(tmp_path, monkeypatch):
    """
    Protects (water-table + duricrust, **Phase 7 restart**): the water table
    `head` and duricrust thickness `duriH` are model memory (they integrate over
    My) and MUST survive restart. A run restarted from step 1 restores `headL`/
    `duriHL` from the output HDF5 (and rebuilds `wtDepth`/`duriF`/`duriKarmor`),
    matching the values written at that step.
    """
    import os
    import shutil
    from gospl.model import Model

    fx = str(FIXTURES_DIR)
    if not os.path.exists(os.path.join(fx, "minimal_gw.yml")):
        pytest.skip("minimal_gw.yml fixture not present")
    monkeypatch.chdir(tmp_path)
    for f in ("mesh.npz", "soiltemp.npz"):
        shutil.copy(os.path.join(fx, f), tmp_path / f)
    base = open(os.path.join(fx, "minimal_gw.yml")).read()
    (tmp_path / "gw.yml").write_text(base)

    # Full run, then read back head/duriH written at step 1.
    m = Model("gw.yml", verbose=False, showlog=False)
    m.runProcesses()
    m.destroy()
    import h5py

    with h5py.File(tmp_path / "gw_out" / "h5" / "gospl.1.p0.h5", "r") as f:
        h1 = np.array(f["wtable"])[:, 0].copy()
        d1 = np.array(f["duricrust"])[:, 0].copy()

    # Restart from step 1 and check the state was restored (not re-initialised).
    (tmp_path / "gwr.yml").write_text(base.replace("start: 0.", "start: 0.\n    rstep: 1"))
    mr = Model("gwr.yml", verbose=False, showlog=False)
    try:
        assert np.allclose(mr.headL.getArray(), h1, atol=1.0e-4), "head not restored"
        assert np.allclose(mr.duriHL.getArray(), d1, atol=1.0e-6), "duriH not restored"
        assert np.allclose(
            mr.wtDepth, mr.hLocal.getArray() - mr.headL.getArray(), atol=1.0e-6
        )
        assert np.allclose(mr.duriKarmor, 1.0 - mr.duriArmorMax * mr.duriF)
    finally:
        mr.destroy()
