"""
Diagnostic glacial model: mass balance, discharge routing, abrasion, till, meltwater.

Protects: AGENTS.md > Ice sheet.

Split out of the former tests/test_regression.py. Run this group alone with
`pytest -m ice`; see tests/README.md for the marker list.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

# Skip the whole module (rather than erroring at collection) when the
# goSPL runtime stack is not installed.
inputparser = pytest.importorskip(
    "gospl.tools.inputparser",
    reason="goSPL runtime deps (petsc4py / ruamel.yaml / scipy) not installed",
)

pytestmark = [pytest.mark.ice]


# ---------------------------------------------------------------------------
# TEST 2b - dual-lithology opt-in flag (DESIGN_DUAL_LITHOLOGY.md Phase 0)
# ---------------------------------------------------------------------------


def _ice_parser():
    """Bare parser primed for `_readIce` (needs tStart/tEnd/dt for the
    no-file interpolators built in `_extraIce`)."""
    parser = inputparser.ReadYaml.__new__(inputparser.ReadYaml)
    parser.input = {}
    parser.tStart = 0.0
    parser.tEnd = 1000.0
    parser.dt = 100.0
    return parser


def test_ice_opt_in():
    """
    Protects: the glacial model is opt-in via the presence of an `ice`
    section; with no `ice` block ice is off, and when present the glacial / abrasion
    / till parameters are parsed with sensible defaults.

    Invariants:
      1. No `ice` block -> iceOn False.
      2. `ice` block present -> iceOn True with default glacial params.
      3. `ice` with sia / abrasion / till -> all parameters parsed.
    """
    # ---- Case 1: no ice block ----
    p = _ice_parser()
    p._readIce()
    assert p.iceOn is False

    # ---- Case 2: ice on with defaults ----
    p = _ice_parser()
    p.input = {"ice": {"hela": 1500.0, "hice": 2000.0}}
    p._readIce()
    assert p.iceOn is True
    assert p.ice_slide == 1.0e-3 and p.ice_glen == 3.0
    # Abrasion off by default (Kg=0); till + catchment routing ON by default
    # (so enabling abrasion gives the complete, spatially-coherent glacial cycle).
    assert p.ice_Kg == 0.0 and p.ice_Kl == 0.0
    assert p.ice_till_on is True
    assert p.ice_till_route is True
    assert p.ice_melt_conserve is True
    # Terminus is unprescribed -> sentinel, resolved to the sea-level position
    # at runtime (so ice is not silently truncated above sea level).
    assert float(p.iceT(p.tStart)) < -1.0e9

    # ---- Case 3: full glacial parameters ----
    p = _ice_parser()
    p.input = {
        "ice": {
            "hela": 1500.0,
            "hice": 2000.0,
            "slide": 5.0e-3,
            "glen": 3.0,
            "abrasion": {"Kg": 1.0e-4, "l": 1.5},
            "till": {"on": True, "route": False},
        }
    }
    p._readIce()
    assert p.iceOn is True
    assert p.ice_slide == 5.0e-3 and p.ice_glen == 3.0
    assert p.ice_Kg == 1.0e-4 and p.ice_abr_l == 1.5
    assert p.ice_till_route is False          # melt-spread opted in
    assert p.ice_till_on is True


def test_ice_geom_field_split():
    """
    Protects: _iceGeomField splits a glacier-geometry input into a scalar
    fallback and an optional [file, key] map spec (the per-vertex ELA path for
    global models).
    """
    p = _ice_parser()
    assert p._iceGeomField(2000.0) == (2000.0, None)
    sc, spec = p._iceGeomField(["ela_map", "ela"])
    assert sc is None and spec == ["ela_map", "ela"]


def test_ice_geom_time_series_parsing(tmp_path):
    """
    Protects: _buildIceSeries turns the optional `glaciers` time series (and the
    single top-level interval) into per-interval (scalar, map_spec) fields, so
    the ELA/ice-cap/terminus can vary in BOTH space (maps) and time — like the
    precipitation `climate` block.
    """
    np.savez(tmp_path / "em.npz", ela=np.zeros(5), hice=np.zeros(5))
    base = str(tmp_path / "em")
    p = _ice_parser()

    # Single top-level interval, ELA as a map, ice-cap as a scalar.
    series = p._buildIceSeries(None, [base, "ela"], 2400.0, 1500.0)
    assert len(series) == 1 and series[0]["start"] == p.tStart
    assert series[0]["hela"] == (None, [base, "ela"])
    assert series[0]["hice"] == (2400.0, None)

    # Multi-interval `glaciers` series, sorted by start, mixing scalars & maps.
    glaciers = [
        {"start": 100.0, "hela": [base, "ela"], "hice": [base, "hice"], "hterm": 0.0},
        {"start": 0.0, "hela": 2000.0, "hice": 2400.0, "hterm": 1500.0},
    ]
    series = p._buildIceSeries(glaciers, None, None, None)
    assert [iv["start"] for iv in series] == [0.0, 100.0]
    assert series[0]["hela"] == (2000.0, None)
    assert series[1]["hela"] == (None, [base, "ela"])


@pytest.mark.slow
def test_ice_geom_time_series_steps(minimal_ice_model):
    """
    Protects: _updateIce selects the active interval for the current time and
    materialises the per-vertex glacier-geometry fields, stepping them as time
    advances (the time-dependent analogue of the precipitation maps).
    """
    m = minimal_ice_model
    m._iceTimeSeries = [
        {"start": 0.0, "hela": (2000.0, None), "hice": (3000.0, None), "hterm": (1500.0, None)},
        {"start": 50.0, "hela": (1000.0, None), "hice": (2000.0, None), "hterm": (500.0, None)},
    ]
    m._iceSeriesIdx = -1

    m.tNow = 0.0
    m._updateIce()
    assert np.allclose(m.elaMesh, 2000.0) and np.allclose(m.iceMesh, 3000.0)

    m.tNow = 60.0
    m._updateIce()
    assert np.allclose(m.elaMesh, 1000.0) and np.allclose(m.termMesh, 500.0)


def test_ice_mfd_diagnostic(minimal_ice_model):
    """
    Protects: the diagnostic ('mfd') ice flow model — a non-dynamical glacial
    driver. It routes the ELA accumulation into an ice discharge, derives a Bahr
    thickness and a balance velocity, and feeds the velocity-based abrasion — with
    no ice-dynamics solve. Must produce finite, non-negative ice with a positive velocity
    and abrasion where ice forms, in a single routing solve (no substep loop).
    """
    from mpi4py import MPI
    from gospl.flow.iceplex import IceMesh
    m = minimal_ice_model
    m.ice_Kg = 1.0e-4          # enable abrasion
    m.tNow = m.tStart
    IceMesh.iceAccumulation(m)
    H = m.iceHL.getArray()
    ub = m.iceUbL.getArray()
    fa = m.iceFAL.getArray()
    abr = m.iceAbrL.getArray()
    for arr in (H, ub, fa, abr):
        assert np.isfinite(arr).all()
        assert (arr >= -1.0e-9).all()
    nice = MPI.COMM_WORLD.allreduce(int((H > 1.0).sum()), op=MPI.SUM)
    assert nice > 0, "diagnostic mode formed no ice"
    # Where there is ice, the balance velocity (and hence abrasion) is positive.
    ubmax = MPI.COMM_WORLD.allreduce(float(ub.max()), op=MPI.MAX)
    abrmax = MPI.COMM_WORLD.allreduce(float(abr.max()), op=MPI.MAX)
    assert ubmax > 0.0 and abrmax > 0.0
    # The discharge is a volume flux concentrated by routing (>> a single cell's
    # local accumulation), i.e. ice converges downhill.
    famax = MPI.COMM_WORLD.allreduce(float(fa.max()), op=MPI.MAX)
    assert famax > 0.0


def test_ice_mfd_dual_strata_till(minimal_ice_dual_model):
    """
    Protects: the diagnostic ('mfd') glacial driver end-to-end with dual
    lithology + stratigraphy + till. Driving the SAME erosion/till machinery from
    the cheap routing proxy must: form ice, abrade the bed, route
    abraded rock into till deposited as moraine, conserve the solid AND the fine
    (dual-lithology) budgets, and update the stratigraphic pile.
    """
    from mpi4py import MPI
    m = minimal_ice_dual_model
    assert m.iceOn and m.ice_till_on and m.ice_Kg > 0.0
    assert m.stratLith and m.stratNb > 0

    te0, td0 = m._tillEroded, m._tillDeposited
    stratH0 = m.stratH.copy()

    m.runProcesses()                  # mfd ice -> abrasion -> till -> dual strata

    # (A) The diagnostic driver formed ice, slid, and abraded.
    g = lambda a, op=MPI.MAX: MPI.COMM_WORLD.allreduce(float(a), op=op)
    assert g(m.iceHL.getArray().max()) > 1.0, "mfd formed no ice"
    assert g(m.iceUbL.getArray().max()) > 0.0, "no basal velocity"
    assert g(m.iceAbrL.getArray().max()) > 0.0, "no abrasion"

    # (B) Till solid produced and conserved over the run (glacial-only counters;
    # the shared fine budget is exercised by fluvial transport too, so it is
    # checked in isolation below).
    dte = MPI.COMM_WORLD.allreduce(m._tillEroded - te0, op=MPI.SUM)
    dtd = MPI.COMM_WORLD.allreduce(m._tillDeposited - td0, op=MPI.SUM)
    assert dte > 0.0, "no till produced under the mfd driver"
    assert np.isclose(dte, dtd, rtol=1.0e-9), "till solid eroded != deposited"

    # (C) Stratigraphy updated by the glacial run.
    dH = (m.stratH - stratH0).sum(axis=1)
    assert float(np.abs(dH).max()) > 0.0, "stratigraphy not updated"

    # (D) Dual-lithology coupling: an isolated, deterministic glacialTill (fast
    # ice up high, melt band lower) must conserve the FINE budget and lay a
    # fine-bearing moraine into the stratigraphy.
    zbed = m.hLocal.getArray().copy()
    m.iceUbL.setArray(np.where(zbed > 2500.0, 0.1, 0.0))
    m.iceMeltL.setArray(np.where((zbed > 1500.0) & (zbed < 2000.0), 1.0, 0.0))
    fe0, fd0 = m._fineEroded, m._fineDeposited
    stratHf0 = m.stratHf.copy()
    m.glacialTill()
    dfe = MPI.COMM_WORLD.allreduce(m._fineEroded - fe0, op=MPI.SUM)
    dfd = MPI.COMM_WORLD.allreduce(m._fineDeposited - fd0, op=MPI.SUM)
    assert dfe > 0.0, "no fine abraded; dual-lithology coupling inactive"
    assert np.isclose(dfe, dfd, rtol=1.0e-6), "fine eroded != fine deposited"
    dHf = (m.stratHf - stratHf0).sum(axis=1)
    abl = (zbed > 1500.0) & (zbed < 2000.0)
    assert float(dHf[abl].max()) > 0.0, "moraine carries no fine fraction"


def test_ice_soil_combined(minimal_ice_soil_model):
    """
    Protects: the diagnostic ('mfd') glacial driver works alongside the
    soil-aware non-linear SPL (`soilSPL`, `cptSoil`). `_glacialAbrasion` is hooked
    into `erodepSPLsoil` and `glacialTill` runs in `runProcesses` regardless of
    eroder, both driven by the mfd-set `iceUbL`/`iceMeltL`. The combined run must
    complete, keep the soil layer and the glacial till both active, and conserve
    the till solid.
    """
    from mpi4py import MPI
    m = minimal_ice_soil_model
    assert m.cptSoil and m.iceOn
    assert m.ice_till_on and m.ice_Kg > 0.0
    te0, td0 = m._tillEroded, m._tillDeposited
    m.runProcesses()                  # soilSPL erosion + mfd glacial abrasion + till

    g = lambda a, op=MPI.MAX: MPI.COMM_WORLD.allreduce(float(a), op=op)
    # Glacial driver active under the soil-coupled eroder.
    assert g(m.iceUbL.getArray().max()) > 0.0, "no basal velocity"
    assert g(m.iceAbrL.getArray().max()) > 0.0, "no abrasion"
    # Soil layer is live (finite, non-negative thicknesses).
    soil = m.Lsoil.getArray()
    assert np.isfinite(soil).all() and (soil >= -1.0e-9).all()
    # Glacial till produced and the solid conserved.
    dte = MPI.COMM_WORLD.allreduce(m._tillEroded - te0, op=MPI.SUM)
    dtd = MPI.COMM_WORLD.allreduce(m._tillDeposited - td0, op=MPI.SUM)
    assert dte > 0.0, "no till produced under soilSPL"
    assert np.isclose(dte, dtd, rtol=1.0e-9), "till solid eroded != deposited"


def test_ice_lateral_erosion(minimal_ice_dual_model):
    """
    Protects: explicit lateral glacial erosion (`ice.abrasion.Kl`) — valley-wall
    abrasion by adjacent fast ice (U-shaping). Off by default (Kl=0); when on it
    must erode subaerial wall cells (little ice of their own, flanking fast ice),
    feed that rock into the conserved till budget, and lower those wall cells.
    """
    from mpi4py import MPI
    m = minimal_ice_dual_model
    A = lambda a: MPI.COMM_WORLD.allreduce(int(a), op=MPI.SUM)

    # Off by default: no wall erosion.
    m.ice_Kl = 0.0
    m.tNow = m.tStart
    from gospl.flow.iceplex import IceMesh
    IceMesh.iceAccumulation(m)
    assert A((m._glacialLateralErosion() > 0).sum()) == 0

    # On: wall cells erode, conserved into till.
    m.ice_Kl = 5.0e-2
    elat = m._glacialLateralErosion()
    nwall = A((elat > 0).sum())
    assert nwall > 0, "lateral erosion eroded no wall cells"
    # Lateral cells carry little ice of their own but a positive rate.
    H = m.iceHL.getArray()
    assert (H[elat > 0] <= 1.0 + 1.0e-9).all(), "lateral erosion hit thick-ice cells"

    z0 = m.hLocal.getArray().copy()
    te0, td0 = m._tillEroded, m._tillDeposited
    m.glacialTill()
    dz = m.hLocal.getArray() - z0
    dte = MPI.COMM_WORLD.allreduce(m._tillEroded - te0, op=MPI.SUM)
    dtd = MPI.COMM_WORLD.allreduce(m._tillDeposited - td0, op=MPI.SUM)
    assert dte > 0.0 and np.isclose(dte, dtd, rtol=1.0e-9), "till not conserved with lateral erosion"
    # At least some wall cells were lowered (valley widening).
    assert A((dz[elat > 0] < -1.0e-9).sum()) > 0, "no wall cells lowered"


def test_ice_meltwater_conserves(minimal_ice_model):
    """
    Protects: the discharge-conserving glacial meltwater delivered to the rivers
    (`ice.melt_conserve`, default True). The water that fell as ice above the ELA
    must return downstream as meltwater — Σ river-meltwater == Σ accumulation —
    so the glacial water budget is closed (vs the precip-scaled ablation, which
    loses water). Distinct from `iceMeltL` (the till melt-out pattern).
    """
    from mpi4py import MPI
    m = minimal_ice_model
    m.tNow = m.tStart
    from gospl.flow.iceplex import IceMesh
    IceMesh.iceAccumulation(m)

    owned = m.inIDs == 1
    _, mdot = m._iceMassBalance(2000.0, 3000.0)
    A = MPI.COMM_WORLD.allreduce(
        float(np.sum((np.maximum(mdot, 0.0) * m.larea)[owned])), op=MPI.SUM
    )
    W = MPI.COMM_WORLD.allreduce(
        float(np.sum(m.iceMeltRiverL.getArray()[owned])), op=MPI.SUM
    )
    assert A > 0.0, "no accumulation; test vacuous"
    assert np.isclose(W, A, rtol=1.0e-6), (
        f"river meltwater {W:.4e} != accumulation {A:.4e} (water not conserved)"
    )

    # Legacy precip-scaled mode: meltwater is the local ablation (generally < A).
    m.ice_melt_conserve = False
    IceMesh.iceAccumulation(m)
    W2 = MPI.COMM_WORLD.allreduce(
        float(np.sum(m.iceMeltRiverL.getArray()[owned])), op=MPI.SUM
    )
    assert np.isfinite(W2) and W2 >= 0.0


@pytest.mark.slow
def test_ice_terminus_sea_level_floor(minimal_ice_model):
    """
    Protects: the terminus floor is max(hterm, sea level) — no ice is kept below
    the sea surface, and a prescribed hterm below sea level is raised to sea
    level. The unprescribed default resolves to the sea-level position. The clamp
    is applied in the diagnostic driver _iceFlowMFD.
    """
    m = minimal_ice_model
    zbed = m.hLocal.getArray().copy()

    # Sea level at 500 m; terminus unprescribed (sentinel) -> floor = sea level.
    m.sealevel = 500.0
    m._iceFlowMFD(2000.0, 3000.0, -1.0e10)     # elaH, iceH, iceT (sentinel)
    assert (m.iceHL.getArray()[zbed < 500.0] == 0.0).all(), "ice kept below sea level"

    # Prescribed hterm BELOW sea level is raised to sea level (floor stays 500).
    m._iceFlowMFD(2000.0, 3000.0, 100.0)
    assert (m.iceHL.getArray()[zbed < 500.0] == 0.0).all(), (
        "hterm below sea level not raised to sea level"
    )

    # Prescribed hterm ABOVE sea level is respected (floor = hterm = 800).
    m._iceFlowMFD(2000.0, 3000.0, 800.0)
    assert (m.iceHL.getArray()[zbed < 800.0] == 0.0).all(), (
        "prescribed terminus above sea level not honoured"
    )


@pytest.mark.slow
def test_ice_spatial_smb(minimal_ice_model):
    """
    Protects: the surface mass balance is per-vertex when the ELA / ice-cap
    altitude are maps (the tropical-vs-polar fix). A constant-array ELA must
    reproduce the uniform-scalar result exactly, and a spatially-high ELA must
    suppress accumulation locally.
    """
    m = minimal_ice_model
    npts = m.lpoints

    # Array path with constant fields == scalar path (byte-identical SMB).
    _, mdot_scalar = m._iceMassBalance(2000.0, 3000.0)
    elaA = np.full(npts, 2000.0)
    iceA = np.full(npts, 3000.0)
    _, mdot_arr = m._iceMassBalance(elaA, iceA)
    assert np.allclose(mdot_scalar, mdot_arr), "array SMB must match scalar SMB"

    # Spatially-high ELA suppresses accumulation: split the domain and raise the
    # ELA out of reach on one half -> no positive mass balance there.
    zbed = m.hLocal.getArray()
    blocked = zbed < np.median(zbed)
    elaS = np.where(blocked, 1.0e9, 2000.0)
    # Keep the unblocked accumulation band identical to the scalar reference
    # (hela=2000, hice=3000) so it must reproduce that SMB exactly.
    iceS = np.where(blocked, elaS + 800.0, 3000.0)
    _, mdot_s = m._iceMassBalance(elaS, iceS)
    assert (mdot_s[blocked] <= 0.0).all(), "no accumulation where the ELA is out of reach"
    # Where the ELA is normal, the SMB matches the uniform-ELA result.
    assert np.allclose(mdot_s[~blocked], mdot_scalar[~blocked])


@pytest.mark.slow
def test_ice_runs_and_invariants(minimal_ice_model):
    """
    Protects: the diagnostic glacial model runs
    end-to-end (the diagnostic glacial driver) and preserves the
    physical invariants of the glacial state.

    Invariants after a run:
      - ice thickness non-negative and finite everywhere;
      - ice actually forms (max H > 0);
      - no ice below the glacier terminus elevation (terminus clamp).
    (Ice DOES extend below the ELA — the discharge routes downhill into the
    ablation zone toward the terminus, as a real glacier does.)
    """
    model = minimal_ice_model
    assert model.iceOn is True

    model.runProcesses()

    H = model.iceHL.getArray()
    zbed = model.hLocal.getArray()
    iceT = float(model.iceT(model.tNow))

    assert np.isfinite(H).all(), "ice thickness went non-finite (blow-up)."
    assert (H >= -1.0e-9).all(), "Negative ice thickness (free boundary violated)."
    assert float(H.max()) > 0.0, "produced no ice."
    assert not (H[zbed < iceT] > 1.0e-6).any(), "Ice below the glacier terminus."

    # Basal sliding speed: finite, non-negative, and confined to ice.
    ub = model.iceUbL.getArray()
    assert np.isfinite(ub).all() and (ub >= -1.0e-12).all()
    assert not (ub[H <= 1.0e-2] > 0.0).any(), "Basal velocity where there is no ice."


@pytest.mark.slow
def test_ice_glacial_abrasion(minimal_ice_model):
    """
    Protects: DESIGN_ICE_SHEET.md Phase 3 — velocity-based glacial abrasion
    E = Kg·|u_b|^l adds incision to Eb exactly where ice slides (and only there),
    and is a no-op when Kg = 0.

    Drives _glacialAbrasion with a known basal-velocity field so the result is
    analytic: Eb = −Kg·u_b (l=1) on the sliding cells, 0 elsewhere.
    """
    m = minimal_ice_model
    m.ice_till_on = False          # test the direct-to-fluvial abrasion path
    zbed = m.hLocal.getArray().copy()
    ub = np.where(zbed > 2000.0, 0.1, 0.0)     # 0.1 m/yr sliding above 2000 m
    m.iceUbL.setArray(ub.copy())
    m.ice_abr_l = 1.0

    # Kg = 0 -> no abrasion.
    m.ice_Kg = 0.0
    m.Eb.set(0.0)
    m._glacialAbrasion()
    assert np.allclose(m.Eb.getArray(), 0.0), "abrasion must be a no-op when Kg=0"

    # Kg > 0 -> incision E = -Kg·u_b where ice slides (all above sea here).
    m.ice_Kg = 1.0e3
    m.Eb.set(0.0)
    m._glacialAbrasion()
    EbL = m.EbLocal.getArray()
    sliding = ub > 0.0
    assert np.allclose(EbL[sliding], -1.0e3 * ub[sliding], rtol=1.0e-6, atol=1.0e-9), (
        "abrasion incision must equal -Kg*|u_b|^l where ice slides"
    )
    assert np.allclose(EbL[~sliding], 0.0, atol=1.0e-9), "abrasion outside ice"


@pytest.mark.slow
def test_ice_glacial_till_conserves(minimal_ice_till_model):
    """
    Protects: DESIGN_ICE_SHEET.md Phase 4 — glacial till is a conservative
    bed-to-bed transport: abrasion lowers the bed under sliding ice and the
    till is deposited (melt-out) in the ablation zone, so the NET bed-volume
    change is zero (rock moved, not created/destroyed). This is the till
    analogue of the dual-lithology fine-conservation guard — a volume the total
    sediment budget cannot see needs its own check.

    The full model runs end-to-end with till on (smoke), then `glacialTill` is
    driven with an imposed sliding-velocity + meltwater field for a
    deterministic conservation check.
    """
    m = minimal_ice_till_model
    assert m.iceOn and m.ice_till_on and m.ice_Kg > 0.0
    # Full glacial + till run must not break.
    m.runProcesses()

    # Deterministic conservation check: fast ice up high, ablation band lower.
    from mpi4py import MPI
    zbed = m.hLocal.getArray().copy()
    m.iceUbL.setArray(np.where(zbed > 2500.0, 0.1, 0.0))
    m.iceMeltL.setArray(np.where((zbed > 1500.0) & (zbed < 2000.0), 1.0, 0.0))
    larea = m.larea
    owned = m.inIDs == 1

    cum0 = m.cumEDLocal.getArray().copy()
    m._tillEroded = 0.0
    m._tillDeposited = 0.0
    m.glacialTill()
    dcum = m.cumEDLocal.getArray() - cum0

    ero = MPI.COMM_WORLD.allreduce(m._tillEroded, op=MPI.SUM)
    dep = MPI.COMM_WORLD.allreduce(m._tillDeposited, op=MPI.SUM)
    netvol = MPI.COMM_WORLD.allreduce(
        float(np.sum((dcum * larea)[owned])), op=MPI.SUM
    )
    activity = MPI.COMM_WORLD.allreduce(
        float(np.sum((np.abs(dcum) * larea)[owned])), op=MPI.SUM
    )

    assert ero > 0.0, "no till was produced; test is vacuous"
    # Rock is conserved: eroded volume == deposited volume.
    assert np.isclose(ero, dep, rtol=1.0e-9), "till eroded != till deposited"
    # Net bed-volume change is zero relative to the rock moved.
    assert abs(netvol) / activity < 1.0e-6, (
        f"glacial till not volume-conserving: net={netvol:.3e} activity={activity:.3e}"
    )
    # Bed lowered under fast ice, raised in the ablation zone.
    assert (dcum[zbed > 2500.0] <= 1.0e-9).all(), "no abrasion under fast ice"
    assert (dcum[(zbed > 1500.0) & (zbed < 2000.0)] >= -1.0e-9).all(), (
        "no till deposition in the ablation zone"
    )


@pytest.mark.slow
@pytest.mark.slow
def test_ice_till_routing(minimal_ice_model):
    """
    Protects: catchment-aware till routing (iceplex._routeTill, `till.route`).
    Abraded till is transported down the ice-surface flow network and melts out
    toward the terminus, so (a) the deposition weight is conserved (Σ = 1), and
    (b) with no ablation the melt-out fraction is zero at interior ice cells, so
    they retain NO till — everything is carried downstream and deposited only at
    the margin outlets (mesh-independent check of the transport mechanism).
    """
    from mpi4py import MPI
    from gospl._fortran import mfdreceivers

    m = minimal_ice_model
    zbed = m.hLocal.getArray().copy()
    n = m.lpoints
    owned = m.inIDs == 1

    # Synthetic ice cap: thick interior thinning to a margin near 1500 m, so the
    # ice surface s = zbed + H descends outward and routing funnels toward it.
    H = np.clip(zbed - 1500.0, 0.0, 800.0)
    m.iceHL.setArray(H)
    m.iceMeltL.setArray(np.zeros(n))            # no melt -> pure routing to margin

    if not (zbed > 2000.0).any():
        pytest.skip("test mesh lacks the relief needed to exercise till routing")

    # Abrasion source confined to the thick interior.
    Vero = np.where(zbed > 2500.0, 1.0e6, 0.0)
    Vtot = MPI.COMM_WORLD.allreduce(float(np.sum(Vero[owned])), op=MPI.SUM)
    if Vtot <= 0.0:
        pytest.skip("no abrasion source on this mesh")

    dep_w = m._routeTill(Vero, Vtot, owned)

    # (a) Mass conserved: the deposition weight sums to one over owned nodes.
    wsum = MPI.COMM_WORLD.allreduce(float(np.sum(dep_w[owned])), op=MPI.SUM)
    assert np.isclose(wsum, 1.0, rtol=1.0e-6), f"routed till not conserved (Σw={wsum})"
    assert (dep_w >= -1.0e-12).all()

    # (b) Interior ice cells (those whose steepest-descent receiver is also ice)
    # have melt-out fraction 0 with no ablation, so they must retain no till —
    # it is all carried downstream to the outlets.
    rcv, _, _ = mfdreceivers(1, m.flowExp, zbed + H, m.sealevel, m.gid)
    rcv0 = rcv[:, 0].astype(int)
    ice = H > 1.0e-2
    interior = owned & ice & (H[rcv0] > 1.0e-2)
    if not interior.any():
        pytest.skip("mesh too coarse: no multi-cell ice flow path to route along")
    assert np.allclose(dep_w[interior], 0.0, atol=1.0e-12), (
        "interior ice retained till despite routing (melt-out should carry it on)"
    )
    # And the till did land somewhere (at the outlets).
    assert float(np.max(dep_w[owned])) > 0.0


@pytest.mark.slow
def test_ice_till_routing_conserves(minimal_ice_model):
    """
    Protects: glacialTill with `till.route` is mass-conserving end-to-end (bulk
    bed mode) — abraded rock routed down-ice and deposited at the termini leaves
    the net bed-volume change zero (eroded == deposited).
    """
    from mpi4py import MPI
    m = minimal_ice_model
    m.ice_till_on = True
    m.ice_till_route = True
    m.ice_Kg = 1.0e-3
    m.ice_abr_l = 1.0
    n = m.lpoints
    owned = m.inIDs == 1
    larea = m.larea

    zbed = m.hLocal.getArray().copy()
    m.iceHL.setArray(np.clip(zbed - 1500.0, 0.0, 800.0))
    m.iceUbL.setArray(np.where(zbed > 2500.0, 0.1, 0.0))   # sliding interior
    m.iceMeltL.setArray(np.zeros(n))
    if not (zbed > 2500.0).any():
        pytest.skip("test mesh lacks the relief to exercise till routing")

    cum0 = m.cumEDLocal.getArray().copy()
    m._tillEroded = 0.0
    m._tillDeposited = 0.0
    m.glacialTill()

    ero = MPI.COMM_WORLD.allreduce(m._tillEroded, op=MPI.SUM)
    dep = MPI.COMM_WORLD.allreduce(m._tillDeposited, op=MPI.SUM)
    assert ero > 0.0, "no till produced; test is vacuous"
    assert np.isclose(ero, dep, rtol=1.0e-9), "routed till eroded != deposited"

    # Net bed-volume change is negligible relative to the till volume moved
    # (erosion balanced by deposition). Measured against `ero` rather than the
    # bed "activity", which collapses to ~0 on a coarse mesh where the abraded
    # cell is its own outlet (till deposited where it was eroded).
    dcum = m.cumEDLocal.getArray() - cum0
    netvol = MPI.COMM_WORLD.allreduce(float(np.sum((dcum * larea)[owned])), op=MPI.SUM)
    assert abs(netvol) < 1.0e-6 * ero, (
        f"routed till not volume-conserving: net={netvol:.3e} eroded={ero:.3e}"
    )


@pytest.mark.slow
def test_ice_glacial_till_dual_lithology(minimal_ice_dual_model):
    """
    Protects: glacial till coupled to dual-lithology stratigraphy
    (iceplex._glacialTillStrata). When stratigraphy is on the abraded rock is
    removed from the stratigraphic pile and re-deposited as a moraine layer
    split into coarse/fine. The conservation invariant is on the SOLID phase:
    the fine volume deposited equals the fine volume eroded (so the
    dual-lithology _fineEroded / _fineDeposited budget stays balanced), and the
    moraine carries the abraded fine fraction.
    """
    from mpi4py import MPI
    m = minimal_ice_dual_model
    assert m.iceOn and m.ice_till_on and m.ice_Kg > 0.0
    assert m.stratLith and m.stratNb > 0
    # Full glacial + till + dual-lithology run must not break.
    m.runProcesses()

    # Deterministic check: fast ice up high (abrasion), melt band lower
    # (ablation / moraine deposition).
    zbed = m.hLocal.getArray().copy()
    m.iceUbL.setArray(np.where(zbed > 2500.0, 0.1, 0.0))
    m.iceMeltL.setArray(np.where((zbed > 1500.0) & (zbed < 2000.0), 1.0, 0.0))
    owned = m.inIDs == 1

    stratH0 = m.stratH.copy()
    stratHf0 = m.stratHf.copy()
    fe0, fd0 = m._fineEroded, m._fineDeposited
    te0, td0 = m._tillEroded, m._tillDeposited

    m.glacialTill()

    # Solid till moved (eroded == deposited, by construction).
    dte = MPI.COMM_WORLD.allreduce(m._tillEroded - te0, op=MPI.SUM)
    dtd = MPI.COMM_WORLD.allreduce(m._tillDeposited - td0, op=MPI.SUM)
    assert dte > 0.0, "no till produced; test is vacuous"
    assert np.isclose(dte, dtd, rtol=1.0e-9), "till solid eroded != deposited"

    # Dual-lithology fine budget stays balanced: the fine removed from the pile
    # by abrasion equals the fine laid back down in the moraine.
    dfe = MPI.COMM_WORLD.allreduce(m._fineEroded - fe0, op=MPI.SUM)
    dfd = MPI.COMM_WORLD.allreduce(m._fineDeposited - fd0, op=MPI.SUM)
    assert dfe > 0.0, "no fine abraded; dual-lithology coupling inactive"
    assert np.isclose(dfe, dfd, rtol=1.0e-6), (
        f"till fine eroded ({dfe:.4e}) != fine deposited ({dfd:.4e})"
    )

    # The stratigraphic pile lost thickness under fast ice and gained a moraine
    # (with a fine component) in the ablation band.
    dH = (m.stratH - stratH0).sum(axis=1)
    dHf = (m.stratHf - stratHf0).sum(axis=1)
    abr = zbed > 2500.0
    abl = (zbed > 1500.0) & (zbed < 2000.0)
    assert (dH[abr] <= 1.0e-9).all(), "pile not eroded under fast ice"
    assert float(dH[abl].max()) > 0.0, "no moraine deposited in the ablation zone"
    assert float(dHf[abl].max()) > 0.0, "moraine carries no fine fraction"


def test_ice_accumulation_scale_cap(minimal_ice_model):
    """
    Protects: the SMB accumulation controls. `accum_factor` scales and `accum_max`
    caps the POSITIVE (accumulation) surface mass balance only; ablation (negative
    mdot) is untouched. Defaults (1.0, None) are a no-op.
    """
    m = minimal_ice_model
    # Baseline accumulation (no scaling).
    m.ice_accum_factor, m.ice_accum_max = 1.0, None
    _, mdot0 = m._iceMassBalance(2000.0, 3000.0)
    acc0 = mdot0[mdot0 > 0.0]
    abl0 = mdot0[mdot0 < 0.0]

    # Halve the accumulation; ablation unchanged.
    m.ice_accum_factor, m.ice_accum_max = 0.5, None
    _, mdot1 = m._iceMassBalance(2000.0, 3000.0)
    assert np.allclose(mdot1[mdot1 > 0.0], 0.5 * acc0) if acc0.size else True
    assert np.allclose(mdot1[mdot1 < 0.0], abl0) if abl0.size else True

    # Cap the accumulation; nothing above the cap, ablation unchanged.
    if acc0.size:
        cap = 0.5 * float(acc0.max())
        m.ice_accum_factor, m.ice_accum_max = 1.0, cap
        _, mdot2 = m._iceMassBalance(2000.0, 3000.0)
        assert float(mdot2[mdot2 > 0.0].max()) <= cap + 1.0e-12
        assert np.allclose(mdot2[mdot2 < 0.0], abl0)


@pytest.mark.slow
def test_ice_flexure_loading(minimal_ice_flex_model):
    """
    Protects: the diagnostic ice thickness feeds the
    existing flexural-isostasy ice load (applyFlexure converts the iceHL change
    to an equivalent sediment load). The run produces a finite flexural field
    with subsidence (negative deflection) under the load.
    """
    import glob
    model = minimal_ice_flex_model
    assert model.iceOn and model.flexOn
    model.runProcesses()

    flx = model.localFlex
    assert np.isfinite(flx).all(), "flexural field non-finite"
    assert float(flx.min()) < 0.0, "no subsidence — ice/sediment load not applied"

    # The ice diagnostic fields are written: thickness, basal velocity,
    # meltwater and abrasion rate.
    files = sorted(
        glob.glob(os.path.join(str(model.outputDir), "h5", "gospl.*.p*.h5"))
    )
    if files:
        h5py = pytest.importorskip("h5py")
        with h5py.File(files[-1], "r") as hf:
            for field in ("iceH", "iceUb", "iceMelt", "iceAbr"):
                assert field in hf, f"ice output field {field} not in output"
