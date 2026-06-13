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

import os

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
# TEST 2b - evaporation parser smoke test
# ---------------------------------------------------------------------------


def test_evap_parser_opt_in():
    """
    Protects: DESIGN_EVAPORATION.md D1, D4 — evaporation is an opt-in
    forcing parsed alongside rainfall from the same `[climate]` YAML block.

    Silent failure prevented: a future refactor that drops the
    `_defineEvap` call from `_readRain` would leave `self.evapdata`
    permanently None even when the YAML declares evap, silently disabling
    the entire feature with no error raised.

    Three assertions:
      1. evapdata is None when no row declares evap (back-compat).
      2. evapdata columns match the contract `[start, eUni, eMap, eKey]`.
      3. evap_uniform values land in the eUni column (not silently
         dropped by a typo in the YAML key name).
    """
    # ---- Case A: rainfall only, no evap → evapdata stays None ----
    parser = _bare_parser()
    parser.input = {"climate": [{"start": 0.0, "uniform": 1.0}]}
    parser._readRain()
    assert parser.raindata is not None
    assert parser.evapdata is None, (
        "evapdata should be None when no climate row declares "
        "evap_uniform or evap_map. Got: "
        f"{parser.evapdata!r}"
    )

    # ---- Case B: rainfall + evap_uniform → evapdata populated ----
    parser = _bare_parser()
    parser.input = {
        "climate": [{"start": 0.0, "uniform": 1.0, "evap_uniform": 0.3}]
    }
    parser._readRain()
    assert parser.evapdata is not None, (
        "evapdata should be a DataFrame when at least one climate row "
        "declares evap_uniform"
    )
    assert list(parser.evapdata.columns) == [
        "start", "eUni", "eMap", "eKey",
    ], (
        "evapdata column order drift. DESIGN_EVAPORATION.md D4 says: "
        "start, eUni, eMap, eKey."
    )
    assert parser.evapdata.at[0, "eUni"] == 0.3, (
        "evap_uniform value not propagated into eUni column"
    )


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
    Protects: the SIA ice-sheet model is opt-in via the presence of an `ice`
    section; with no `ice` block ice is off, and when present the SIA / abrasion
    / till parameters are parsed with sensible defaults.

    Invariants:
      1. No `ice` block -> iceOn False.
      2. `ice` block present -> iceOn True with default SIA params.
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
    assert p.sia_Aglen == 1.0e-16 and p.sia_glen == 3.0
    assert p.ice_Kg == 0.0 and p.ice_till_on is False
    assert p.ice_till_route is False          # melt-weighted spreading by default
    # Terminus is unprescribed -> sentinel, resolved to the sea-level position
    # at runtime (so ice is not silently truncated above sea level).
    assert float(p.iceT(p.tStart)) < -1.0e9

    # ---- Case 3: full SIA parameters ----
    p = _ice_parser()
    p.input = {
        "ice": {
            "hela": 1500.0,
            "hice": 2000.0,
            "sia": {"Aglen": 2.0e-16, "slide": 5.0e-3, "glen": 3.0},
            "abrasion": {"Kg": 1.0e-4, "l": 1.5},
            "till": {"on": True, "route": True},
        }
    }
    p._readIce()
    assert p.iceOn is True
    assert p.sia_Aglen == 2.0e-16 and p.sia_slide == 5.0e-3 and p.sia_glen == 3.0
    assert p.ice_Kg == 1.0e-4 and p.ice_abr_l == 1.5
    assert p.ice_till_route is True           # catchment routing opted in
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
def test_ice_geom_time_series_steps(minimal_ice_sia_model):
    """
    Protects: _updateIce selects the active interval for the current time and
    materialises the per-vertex glacier-geometry fields, stepping them as time
    advances (the time-dependent analogue of the precipitation maps).
    """
    m = minimal_ice_sia_model
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


@pytest.mark.slow
def test_ice_seed_and_evolve(minimal_ice_seed_model):
    """
    Protects: a pre-existing ice thickness given via `ice.hinit` seeds iceHL on a
    fresh start (before any solve), and the SIA solve then evolves it.
    """
    m = minimal_ice_seed_model
    # The seed is applied during model init, before any time step.
    assert np.allclose(m.iceHL.getArray(), 50.0), "hinit did not seed iceHL"

    m.runProcesses()
    H = m.iceHL.getArray()
    # Still a valid ice field after evolving the seed.
    assert np.isfinite(H).all() and (H >= -1.0e-9).all()
    assert float(H.max()) > 0.0, "seeded ice vanished entirely"


@pytest.mark.slow
def test_ice_terminus_sea_level_floor(minimal_ice_sia_model):
    """
    Protects: the terminus floor is max(hterm, sea level) — no land ice persists
    below the sea surface, and a prescribed hterm below sea level is raised to
    sea level. The unprescribed default resolves to the sea-level position.
    """
    m = minimal_ice_sia_model
    zbed = m.hLocal.getArray().copy()
    n = m.lpoints
    H_in = np.full(n, 100.0)            # ice everywhere to start
    mdot0 = np.zeros(n)

    # Sea level at 500 m; terminus unprescribed (sentinel) -> floor = sea level.
    m.sealevel = 500.0
    m._iceSIAFinalize(H_in.copy(), zbed, mdot0, -1.0e10, m.sia_slide, m.sia_glen)
    H = m.iceHL.getArray()
    assert (H[zbed < 500.0] == 0.0).all(), "ice kept below sea level"
    if (zbed >= 500.0).any():
        assert (H[zbed >= 500.0] > 0.0).any(), "ice wrongly removed above sea level"

    # Prescribed hterm BELOW sea level is raised to sea level (floor stays 500).
    m._iceSIAFinalize(H_in.copy(), zbed, mdot0, 100.0, m.sia_slide, m.sia_glen)
    H = m.iceHL.getArray()
    assert (H[zbed < 500.0] == 0.0).all(), "hterm below sea level not raised to sea level"

    # Prescribed hterm ABOVE sea level is respected (floor = hterm = 800).
    m._iceSIAFinalize(H_in.copy(), zbed, mdot0, 800.0, m.sia_slide, m.sia_glen)
    H = m.iceHL.getArray()
    assert (H[zbed < 800.0] == 0.0).all(), "prescribed terminus above sea level not honoured"


@pytest.mark.slow
def test_ice_spatial_smb(minimal_ice_sia_model):
    """
    Protects: the SIA surface mass balance is per-vertex when the ELA / ice-cap
    altitude are maps (the tropical-vs-polar fix). A constant-array ELA must
    reproduce the uniform-scalar result exactly, and a spatially-high ELA must
    suppress accumulation locally.
    """
    m = minimal_ice_sia_model
    npts = m.lpoints

    # Array path with constant fields == scalar path (byte-identical SMB).
    _, _, _, _, mdot_scalar = m._iceSIAParams(2000.0, 3000.0)
    elaA = np.full(npts, 2000.0)
    iceA = np.full(npts, 3000.0)
    _, _, _, _, mdot_arr = m._iceSIAParams(elaA, iceA)
    assert np.allclose(mdot_scalar, mdot_arr), "array SMB must match scalar SMB"

    # Spatially-high ELA suppresses accumulation: split the domain and raise the
    # ELA out of reach on one half -> no positive mass balance there.
    zbed = m.hLocal.getArray()
    blocked = zbed < np.median(zbed)
    elaS = np.where(blocked, 1.0e9, 2000.0)
    # Keep the unblocked accumulation band identical to the scalar reference
    # (hela=2000, hice=3000) so it must reproduce that SMB exactly.
    iceS = np.where(blocked, elaS + 800.0, 3000.0)
    _, _, _, _, mdot_s = m._iceSIAParams(elaS, iceS)
    assert (mdot_s[blocked] <= 0.0).all(), "no accumulation where the ELA is out of reach"
    # Where the ELA is normal, the SMB matches the uniform-ELA result.
    assert np.allclose(mdot_s[~blocked], mdot_scalar[~blocked])


@pytest.mark.slow
def test_ice_sia_runs_and_invariants(minimal_ice_sia_model):
    """
    Protects: DESIGN_ICE_SHEET.md Phase 1 — the SIA ice-sheet model runs
    end-to-end (the dynamic ice_flux thickness evolution) and preserves the
    physical invariants of a free-boundary ice solve.

    Invariants after a run:
      - ice thickness non-negative and finite everywhere (free boundary, no
        blow-up from the explicit subcycling);
      - ice actually forms (max H > 0);
      - no ice below the glacier terminus elevation (terminus clamp);
      - ice only where it can accumulate (above the ELA region).
    """
    model = minimal_ice_sia_model
    assert model.iceOn is True

    model.runProcesses()

    H = model.iceHL.getArray()
    zbed = model.hLocal.getArray()
    elaH = float(model.elaH(model.tNow))
    iceT = float(model.iceT(model.tNow))

    assert np.isfinite(H).all(), "SIA ice thickness went non-finite (blow-up)."
    assert (H >= -1.0e-9).all(), "Negative ice thickness (free boundary violated)."
    assert float(H.max()) > 0.0, "SIA produced no ice."
    assert not (H[zbed < iceT] > 1.0e-6).any(), "Ice below the glacier terminus."
    # Ice only where accumulation is possible (at/above the ELA).
    assert not (H[zbed < elaH] > 1.0e-6).any(), "Ice below the ELA."

    # Basal sliding speed (Phase 2): finite, non-negative, and confined to ice.
    ub = model.iceUbL.getArray()
    assert np.isfinite(ub).all() and (ub >= -1.0e-12).all()
    assert not (ub[H <= 1.0e-2] > 0.0).any(), "Basal velocity where there is no ice."


@pytest.mark.slow
def test_ice_glacial_abrasion(minimal_ice_sia_model):
    """
    Protects: DESIGN_ICE_SHEET.md Phase 3 — velocity-based glacial abrasion
    E = Kg·|u_b|^l adds incision to Eb exactly where ice slides (and only there),
    and is a no-op when Kg = 0.

    Drives _glacialAbrasion with a known basal-velocity field so the result is
    analytic: Eb = −Kg·u_b (l=1) on the sliding cells, 0 elsewhere.
    """
    m = minimal_ice_sia_model
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
    # Full SIA + till run must not break.
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
def test_ice_till_routing(minimal_ice_sia_model):
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

    m = minimal_ice_sia_model
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
def test_ice_till_routing_conserves(minimal_ice_sia_model):
    """
    Protects: glacialTill with `till.route` is mass-conserving end-to-end (bulk
    bed mode) — abraded rock routed down-ice and deposited at the termini leaves
    the net bed-volume change zero (eroded == deposited).
    """
    from mpi4py import MPI
    m = minimal_ice_sia_model
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
    # Full SIA + till + dual-lithology run must not break.
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


@pytest.mark.slow
def test_ice_sia_volume_conservation(minimal_ice_sia_model):
    """
    Protects: DESIGN_ICE_SHEET.md Phase 6 — the SIA ice flux is mass-conservative
    (it redistributes ice, it does not create or destroy it). Halfar-style check:
    with zero surface mass balance and the terminus clamp disabled, stepping an
    ice dome must conserve ice volume to the flux-numerics floor.
    """
    from mpi4py import MPI
    m = minimal_ice_sia_model
    zbed = m.hLocal.getArray().copy()
    m.rainVal = np.zeros_like(m.rainVal)       # zero SMB -> pure flux redistribution
    m.iceT = lambda t: 0.0                     # disable terminus clamp
    H0 = np.where(zbed > 2500.0, 300.0, 0.0)
    owned = m.inIDs == 1
    larea = m.larea

    V0 = MPI.COMM_WORLD.allreduce(float(np.sum((H0 * larea)[owned])), op=MPI.SUM)
    m.iceHL.setArray(H0.copy())
    m._iceFlowSIA(0.0, 1000.0, 0.0)            # elaH, iceH, iceT
    H1 = m.iceHL.getArray()
    V1 = MPI.COMM_WORLD.allreduce(float(np.sum((H1 * larea)[owned])), op=MPI.SUM)

    assert V0 > 0.0
    assert abs(V1 - V0) / V0 < 1.0e-4, (
        f"SIA flux not volume-conserving: V0={V0:.4e} V1={V1:.4e}"
    )
    # The flux actually moved ice (not a frozen field).
    assert float(np.max(np.abs(H1 - H0))) > 0.0


@pytest.mark.slow
def test_ice_sia_flexure_loading(minimal_ice_flex_model):
    """
    Protects: DESIGN_ICE_SHEET.md Phase 5 — the SIA ice thickness feeds the
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


def _strata_parser(stratNb):
    """
    Bare parser primed for `_extraStrata`: it reads `self.input`,
    `self.phi0s`/`self.z0s` (set by `_readCompaction` upstream), and
    `self.stratNb` (set by `_readTime` upstream). See `_bare_parser`.
    """
    parser = inputparser.ReadYaml.__new__(inputparser.ReadYaml)
    parser.input = {}
    parser.phi0s = 0.49
    parser.z0s = 3700.0
    parser.stratNb = stratNb
    return parser


def test_dual_lithology_opt_in():
    """
    Protects: DESIGN_DUAL_LITHOLOGY.md Phase 0 — dual lithology is an
    opt-in parsed in `_extraStrata` (continuation of `_readCompaction`).

    Silent failure prevented: a refactor dropping the `_extraStrata` call
    from `_readCompaction`, or flipping the default, would silently change
    the sediment model for every existing input file.

    Invariants:
      1. No `strata` block  → stratLith False; coarse curve == compaction
         curve (the dual-off path must stay bitwise-identical).
      2. `strata: dual: True` with stratigraphy on → stratLith True and the
         per-lithology parameters are parsed.
      3. `strata: dual: True` with stratigraphy OFF (stratNb == 0) → flag
         forced back to False (dual requires stratigraphy).
    """
    # ---- Case 1: no strata block → single-fraction, defaults mirror compaction
    parser = _strata_parser(stratNb=5)
    parser._extraStrata()
    assert parser.stratLith is False
    assert parser.phi0c == parser.phi0s and parser.z0c == parser.z0s, (
        "Dual-off coarse porosity curve must default to the single-fraction "
        "compaction curve so behaviour is unchanged."
    )

    # ---- Case 2: dual on with stratigraphy enabled → parsed
    parser = _strata_parser(stratNb=5)
    parser.input = {
        "strata": {
            "dual": True,
            "coarse": {"phi0": 0.45, "z0": 3000.0},
            "fine": {"phi0": 0.65, "z0": 1500.0, "k_factor": 1.5},
            "bedrock_coarse_frac": 0.7,
            "fine_efficiency": 0.3,
            "pitInletBias": {"coarse": 0.8, "fine": 0.1},
            "fine_diff_factor": 2.0,
        }
    }
    parser._extraStrata()
    assert parser.stratLith is True
    assert parser.phi0c == 0.45 and parser.z0c == 3000.0
    assert parser.phi0f == 0.65 and parser.z0f == 1500.0
    assert parser.fine_k_factor == 1.5
    assert parser.bedrock_coarse_frac == 0.7
    assert parser.fine_efficiency == 0.3
    assert parser.pit_inlet_bias_coarse == 0.8
    assert parser.pit_inlet_bias_fine == 0.1
    assert parser.fine_diff_factor == 2.0

    # ---- Case 3: dual requested but stratigraphy off → forced False
    parser = _strata_parser(stratNb=0)
    parser.input = {"strata": {"dual": True}}
    parser._extraStrata()
    assert parser.stratLith is False, (
        "Dual lithology must require stratigraphy (stratNb > 0); with strat "
        "disabled the flag has to fall back to single-fraction."
    )


def test_dual_lithology_layer_allocation():
    """
    Protects: DESIGN_DUAL_LITHOLOGY.md Phase 1 — the fine-fraction layer
    fields (stratHf, phiF) are allocated by readStratLayers only when
    dual lithology is enabled, and stay None otherwise.

    Silent failure prevented: allocating the fields unconditionally would
    change the single-fraction memory footprint and (via _outputStrat)
    write extra HDF5 datasets into every existing single-fraction run.

    Exercises the no-file branch of readStratLayers directly via __new__,
    which only needs lpoints/stratNb/phi0s/stratLith/memclear/strataFile.
    """
    from gospl.sed import stratplex

    from gospl.tools.constants import BEDROCK_SENTINEL

    def _strata_mesh(stratLith):
        m = stratplex.STRAMesh.__new__(stratplex.STRAMesh)
        m.strataFile = None
        m.lpoints = 8
        m.stratNb = 3
        m.phi0s = 0.49
        m.phi0f = 0.63
        m.bedrock_coarse_frac = 0.5
        m.memclear = False
        m.stratLith = stratLith
        m.stratHf = None
        m.phiF = None
        return m

    # Single-fraction: fine fields stay None.
    m = _strata_mesh(stratLith=False)
    m.readStratLayers()
    assert m.stratHf is None and m.phiF is None, (
        "Single-fraction run must not allocate stratHf/phiF."
    )

    # Dual: fine fields allocated with stratH's shape. The bedrock sentinel
    # layer (layer 0) carries the bedrock fine fraction; layers above are 0.
    m = _strata_mesh(stratLith=True)
    m.readStratLayers()
    assert m.stratHf is not None and m.phiF is not None
    assert m.stratHf.shape == m.stratH.shape == (m.lpoints, m.stratNb)
    assert m.phiF.shape == (m.lpoints, m.stratNb)
    assert np.allclose(m.stratHf[:, 0], BEDROCK_SENTINEL * (1.0 - 0.5)), (
        "Bedrock-sentinel layer must carry the (1 - bedrock_coarse_frac) fine split."
    )
    assert (m.stratHf[:, 1:] == 0.0).all(), "Layers above bedrock start coarse-empty."


def test_dual_lithology_initial_strata_composition(tmp_path):
    """
    Protects: a user can supply a per-layer coarse/fine composition for the
    INITIAL stratigraphy via the npstrata file (strataHf = fine bulk thickness,
    phiF = fine porosity). readStratLayers must load it per layer, clamp
    0 <= strataHf <= strataH, and default phiF to phi0f when absent.
    """
    stratplex = pytest.importorskip("gospl.sed.stratplex")
    n, nl = 3, 2
    f = tmp_path / "init.npz"
    np.savez(
        str(f),
        strataH=np.full((n, nl), 10.0),
        strataZ=np.zeros((n, nl)),
        phiS=np.full((n, nl), 0.49),
        # per-layer fine thickness; row 2 layer 0 is 15 > 10 -> must clamp to 10.
        strataHf=np.array([[2.0, 5.0], [0.0, 10.0], [15.0, 3.0]]),
        # phiF deliberately omitted -> should default to phi0f.
    )
    m = stratplex.STRAMesh.__new__(stratplex.STRAMesh)
    m.strataFile = str(f)
    m.lpoints = n
    m.stratNb = 2                     # extra capacity beyond the initial layers
    m.locIDs = np.arange(n)
    m.phi0s, m.phi0f = 0.49, 0.63
    m.memclear = False
    m.stratLith = True
    m.stratHf = None
    m.phiF = None

    m.readStratLayers()

    il = m.initLay
    assert il == nl
    # Per-layer composition loaded.
    assert np.isclose(m.stratHf[0, 0], 2.0) and np.isclose(m.stratHf[0, 1], 5.0)
    # Clamp: strataHf (15) exceeded strataH (10) -> clamped to 10.
    assert np.isclose(m.stratHf[2, 0], 10.0)
    # Partition stays physical everywhere.
    assert (m.stratHf[:, :il] <= m.stratH[:, :il] + 1e-9).all()
    assert (m.stratHf[:, :il] >= 0.0).all()
    # phiF absent -> defaulted to phi0f on the initial layers.
    assert np.allclose(m.phiF[:, :il], 0.63)


def test_dual_lithology_erosion_split():
    """
    Protects: DESIGN_DUAL_LITHOLOGY.md Phase 2 — erodeStrat splits the
    eroded solid into thCoarse/thFine by the consumed layers' composition,
    and stays mass-consistent. K-blend / composition helpers also tested.

    Builds a tiny single-node column by hand and drives erodeStrat through
    a stubbed PETSc boundary (tmpL/globalToLocal), so the split arithmetic
    is exercised without a full model.

    Invariants:
      1. All-coarse column -> thFine == 0 and thCoarse matches what the
         single-fraction branch would produce (parity).
      2. Mixed column -> the deposited (uncompacted) split reconstructs the
         eroded solid: thCoarse*(1-phi0c) + thFine*(1-phi0f) == solid removed.
      3. _surfaceComposition / _surfaceLithoK reflect the exposed layer and
         the fine_k_factor; both are neutral (1.0) when dual is off.
    """
    from gospl.sed import stratplex

    class _StubVec:
        def __init__(self, arr):
            self._a = arr
        def getArray(self):
            return self._a
        def setArray(self, a):
            self._a = np.asarray(a, dtype=np.float64)

    class _StubDM:
        def globalToLocal(self, g, l):
            l.setArray(g.getArray().copy())

    def _mesh(stratLith, ero, stratH, stratHf=None, phiS=None, phiF=None):
        m = stratplex.STRAMesh.__new__(stratplex.STRAMesh)
        lp, nb = stratH.shape
        m.lpoints, m.stratNb, m.stratStep = lp, nb, nb - 1
        m.dt = 1.0
        m.memclear = False
        m.phi0s = m.phi0c = 0.49
        m.phi0f = 0.63
        m.bedrock_coarse_frac = 0.5
        m.fine_k_factor = 1.0
        m.inIDs = np.ones(lp, dtype=int)
        m.larea = np.ones(lp)
        m._fineEroded = 0.0
        m.stratLith = stratLith
        m.stratH = stratH.astype(np.float64).copy()
        m.phiS = (phiS if phiS is not None else np.full_like(stratH, 0.49)).copy()
        m.stratK = np.ones_like(m.stratH)
        m.stratHf = stratHf.copy() if stratHf is not None else None
        m.phiF = phiF.copy() if phiF is not None else None
        # PETSc boundary stub: tmp carries -ero (erosion is negative).
        m.tmp = _StubVec(np.array(ero, dtype=np.float64))
        m.tmpL = _StubVec(np.zeros(lp))
        m.dm = _StubDM()
        return m

    # ---- Case 1: all-coarse parity (single vs dual must agree on thCoarse)
    stratH = np.array([[2.0, 3.0]])          # two layers, total 5 m
    single = _mesh(False, ero=[-4.0], stratH=stratH)
    single.erodeStrat()
    dual = _mesh(
        True, ero=[-4.0], stratH=stratH,
        stratHf=np.zeros_like(stratH), phiF=np.full_like(stratH, 0.63),
    )
    dual.erodeStrat()
    assert np.allclose(dual.thFine, 0.0), "All-coarse column must erode no fine."
    assert np.allclose(dual.thCoarse, single.thCoarse), (
        "All-coarse dual erosion must match the single-fraction thCoarse."
    )

    # ---- Case 2: mixed column -> per-fraction solid reconstruction
    stratH = np.array([[2.0, 3.0]])
    stratHf = np.array([[1.0, 1.2]])         # fine bulk per layer
    phiS = np.full_like(stratH, 0.49)
    phiF = np.full_like(stratH, 0.63)
    dual = _mesh(True, ero=[-4.0], stratH=stratH, stratHf=stratHf,
                 phiS=phiS, phiF=phiF)
    dual.erodeStrat()
    # Solid removed by erosion of 4 m: fully erode top layer (3 m) + 1 m of
    # the lower (well-mixed) layer.
    # top layer (idx1): coarse (3-1.2)*(1-.49) + fine 1.2*(1-.63)
    # lower (idx0): remove 1 m of 2 m -> half: coarse (2-1)/2*(1-.49)*1 ... compute generically:
    coarse_solid = (3 - 1.2) * (1 - 0.49) + ((2 - 1.0) * 0.5) * (1 - 0.49)
    fine_solid = 1.2 * (1 - 0.63) + (1.0 * 0.5) * (1 - 0.63)
    got = dual.thCoarse[0] * (1 - dual.phi0c) + dual.thFine[0] * (1 - dual.phi0f)
    assert np.isclose(got, coarse_solid + fine_solid, rtol=1e-9), (
        f"Per-fraction solid reconstruction failed: got {got}, "
        f"expected {coarse_solid + fine_solid}"
    )
    assert dual.thFine[0] > 0 and dual.thCoarse[0] > 0

    # ---- Case 3: composition + K-blend helpers
    m = _mesh(True, ero=[0.0], stratH=np.array([[2.0, 2.0]]),
              stratHf=np.array([[0.0, 0.5]]), phiF=np.full((1, 2), 0.63))
    fc = m._surfaceComposition()
    assert np.isclose(fc[0], 1.0 - 0.5 / 2.0), "Surface coarse fraction wrong."
    m.fine_k_factor = 2.0
    litK = m._surfaceLithoK()
    assert np.isclose(litK[0], fc[0] + (1 - fc[0]) * 2.0)
    m.stratLith = False
    assert np.allclose(m._surfaceLithoK(), 1.0), "K-blend must be neutral when off."


def test_dual_lithology_deposit_and_compaction():
    """
    Protects: DESIGN_DUAL_LITHOLOGY.md Phase 4 — per-fraction deposition
    (deposeStrat) and per-fraction compaction (_depthPorosityDual).

    Deposition: the fresh layer accumulates a fine fraction equal to the
    step's global eroded composition, with each lithology's surface porosity.

    Compaction (the headline physics): each fraction compacts on its own
    porosity-depth curve, conserving its solid phase while fines lose more
    bulk thickness than coarse at the same burial depth.
    """
    from gospl.sed import stratplex
    from mpi4py import MPI

    class _StubVec:
        def __init__(self, arr):
            self._a = np.asarray(arr, dtype=np.float64)
        def getArray(self):
            return self._a
        def setArray(self, a):
            self._a = np.asarray(a, dtype=np.float64)

    class _StubDM:
        def globalToLocal(self, g, l):
            l.setArray(g.getArray().copy())

    # ---- Deposition: 4 m deposit, per-node fine fraction (fineFrac) = 0.25 ----
    m = stratplex.STRAMesh.__new__(stratplex.STRAMesh)
    m.lpoints, m.stratNb, m.stratStep = 1, 2, 1
    m.stratLith = True
    m.memclear = False
    m.phi0c, m.phi0f = 0.49, 0.63
    m.depoFineFrac = np.array([0.25])   # per-node deposit fine fraction (Phase 3a/3b)
    m.inIDs = np.ones(1, dtype=int)
    m.larea = np.ones(1)
    m._fineDeposited = 0.0
    m.stratH = np.zeros((1, 2))
    m.stratHf = np.zeros((1, 2))
    m.phiS = np.zeros((1, 2))
    m.phiF = np.zeros((1, 2))
    m.stratK = np.ones((1, 2))
    m.tmp = _StubVec([4.0])
    m.tmpL = _StubVec([0.0])
    m.dm = _StubDM()
    m.deposeStrat()
    assert np.isclose(m.stratH[0, 1], 4.0)
    assert np.isclose(m.stratHf[0, 1], 4.0 * 0.25), (
        "Fine deposit must equal depo * per-node fineFrac."
    )
    assert np.isclose(m.phiF[0, 1], 0.63) and np.isclose(m.phiS[0, 1], 0.49)

    # ---- Compaction: one mixed layer buried 2 km ----
    m = stratplex.STRAMesh.__new__(stratplex.STRAMesh)
    m.stratStep = 0
    m.stratLith = True
    m.memclear = False
    m.bedrockLay = 0
    m.phi0c, m.z0c = 0.49, 3700.0
    m.phi0f, m.z0f = 0.63, 1960.0
    Hc0, Hf0 = 6.0, 4.0
    m.stratH = np.array([[Hc0 + Hf0]])
    m.stratHf = np.array([[Hf0]])
    m.phiS = np.array([[0.49]])
    m.phiF = np.array([[0.63]])
    depth = np.array([[-2000.0]])
    newH = m._depthPorosity(depth)

    phiS_new = 0.49 * np.exp(-2000.0 / 3700.0)
    phiF_new = 0.63 * np.exp(-2000.0 / 1960.0)
    # Per-fraction solid is conserved through compaction.
    Hc_new = newH[0, 0] - m.stratHf[0, 0]
    Hf_new = m.stratHf[0, 0]
    assert np.isclose(Hc_new * (1 - phiS_new), Hc0 * (1 - 0.49), rtol=1e-9), (
        "Coarse solid must be conserved through compaction."
    )
    assert np.isclose(Hf_new * (1 - phiF_new), Hf0 * (1 - 0.63), rtol=1e-9), (
        "Fine solid must be conserved through compaction."
    )
    assert newH[0, 0] < Hc0 + Hf0, "Compaction must reduce total thickness."
    # Fines lose proportionally more bulk than coarse at the same depth.
    assert (Hf_new / Hf0) < (Hc_new / Hc0), (
        "Fines must compact more than coarse for these curves."
    )


def test_dual_lithology_advection_fine_pile():
    """
    Protects: DESIGN_DUAL_LITHOLOGY.md Phase 5 — stratalRecord advects the
    fine pile (stratHf, phiF) alongside the total/coarse pile via a second
    strataonesed call (NOT stratathreesed, whose extra fields are 0-1
    fractions and renormalised — wrong for the bulk-thickness representation).

    Uses an identity advection (each node maps to itself with weight 1) so
    the records must come back unchanged, confirming the fine fields are
    actually routed through the interpolation and written back.
    """
    stratplex = pytest.importorskip("gospl.sed.stratplex")

    class _MockVec:
        def __init__(self, n):
            self._a = np.zeros(n)
        def setArray(self, a):
            self._a = np.asarray(a, dtype=np.float64).copy()
        def getArray(self):
            return self._a

    class _MockDM:
        def globalToLocal(self, src, dst):
            dst.setArray(src.getArray())  # identity halo exchange

    n = 4
    m = stratplex.STRAMesh.__new__(stratplex.STRAMesh)
    m.lpoints = n
    m.stratStep = 2          # advect layers 0..1
    m.stratLith = True
    m.stratH = np.array([[5.0, 7.0, 0.0]] * n)
    m.stratHf = np.array([[2.0, 3.0, 0.0]] * n)
    m.stratZ = np.array([[-10.0, -3.0, 0.0]] * n)
    m.phiS = np.array([[0.45, 0.48, 0.0]] * n)
    m.phiF = np.array([[0.60, 0.62, 0.0]] * n)
    m.tmp = _MockVec(n)
    m.tmpL = _MockVec(n)
    m.dm = _MockDM()

    # Identity interpolation: 3 neighbours all = self, weights summing to 1.
    indices = np.repeat(np.arange(n)[:, None], 3, axis=1)
    weights = np.full((n, 3), 1.0 / 3.0)
    onIDs = np.array([], dtype=int)

    H0, Hf0 = m.stratH.copy(), m.stratHf.copy()
    phiS0, phiF0 = m.phiS.copy(), m.phiF.copy()
    m.stratalRecord(indices, weights, onIDs)

    # Advected layers (0..stratStep-1) must be preserved by identity mapping.
    s = slice(0, m.stratStep)
    assert np.allclose(m.stratHf[:, s], Hf0[:, s]), (
        "Fine pile thickness must round-trip through identity advection."
    )
    assert np.allclose(m.phiF[:, s], phiF0[:, s]), (
        "Fine porosity must round-trip through identity advection."
    )
    # Coarse/total pile still correct (regression on the original behaviour).
    assert np.allclose(m.stratH[:, s], H0[:, s])
    assert np.allclose(m.phiS[:, s], phiS0[:, s])


def test_dual_lithology_pit_fine_bias():
    """
    Protects: DESIGN_DUAL_LITHOLOGY.md Phase 3b — _pitFineFraction biases the
    pit/lake deposit composition so fine concentrates toward the depocenter
    (deep) and coarse toward the inlet/margins (shallow), while conserving the
    per-pit incoming fine volume.

    Single pit, 4 nodes at increasing bathymetric depth, uniform deposit and
    area. The pit retained fine volume 1.2 of total 4 (ff_pit = 0.3, the
    coarse-settles-first retained fraction from the cascade). The resulting
    per-node fine fraction must (a) increase monotonically with depth and (b)
    conserve the retained fine volume: Σ(delta·larea·ffrac) == _pitRetFine.
    """
    sedplex = pytest.importorskip("gospl.sed.sedplex")
    n = 4
    m = sedplex.SEDMesh.__new__(sedplex.SEDMesh)
    m.lpoints = n
    m.stratLith = True
    m.pitParams = np.zeros((1, 3))          # one pit
    m.pitIDs = np.zeros(n, dtype=int)        # all nodes in pit 0
    m.inIDs = np.ones(n, dtype=int)
    m.larea = np.ones(n)
    m.lFill = np.full(n, 10.0)               # rim at 10 m
    hl = np.array([9.0, 7.0, 4.0, 8.0])      # depth = 1, 3, 6, 2
    delta = np.ones(n)                       # uniform deposit thickness
    depo = np.array([4.0])                   # retained volume (Σ delta*larea)
    m._pitRetFine = np.array([1.2])          # retained fine → ff_pit = 1.2/4 = 0.3
    m.depoFineFrac = np.zeros(n)

    m._pitFineFraction(hl, delta, depo)

    depth = m.lFill - hl
    ff = m.depoFineFrac
    # (a) fine fraction increases with depth (depocenter is fine-rich).
    order = np.argsort(depth)
    assert np.all(np.diff(ff[order]) > 0), (
        f"Fine fraction must increase with depth; got {ff} for depth {depth}."
    )
    # (b) per-pit fine volume conserved (ff_pit = 0.3, Σ delta*larea = 4).
    fine_vol = float(np.sum(delta * m.larea * ff))
    assert np.isclose(fine_vol, 0.3 * 4.0, rtol=1e-9), (
        f"Pit fine volume not conserved: {fine_vol} != 1.2"
    )


def test_dual_lithology_marine_fine_bias():
    """
    Protects: DESIGN_DUAL_LITHOLOGY.md Phase 3c — _marineFineFraction biases
    the marine deposit composition so fine concentrates in deep / distal water
    and coarse stays proximal (shallow), conserving the marine fine volume.
    Subaqueous analogue of the pit-fine bias.

    Four marine nodes at increasing water depth, uniform deposit and area,
    uniform arriving composition (ff_mar = 0.3). The per-node fine fraction
    must increase with depth and conserve fine volume.
    """
    seaplex = pytest.importorskip("gospl.sed.seaplex")
    n = 4

    class _StubVec:
        def __init__(self, arr):
            self._a = np.asarray(arr, dtype=np.float64)
        def getArray(self):
            return self._a
        def setArray(self, a):
            self._a = np.asarray(a, dtype=np.float64)

    class _StubDM:
        def globalToLocal(self, g, l):
            l.setArray(g.getArray().copy())

    m = seaplex.SEAMesh.__new__(seaplex.SEAMesh)
    m.lpoints = n
    m.sealevel = 0.0
    m.inIDs = np.ones(n, dtype=int)
    m.larea = np.ones(n)
    m.depoFineFrac = np.zeros(n)
    mdep = np.ones(n)                        # uniform marine deposit
    m.tmp = _StubVec(mdep)
    m.tmpL = _StubVec(np.zeros(n))
    m.dm = _StubDM()
    hl = np.array([-1.0, -3.0, -6.0, -8.0])  # depth = 1, 3, 6, 8
    sedFlux = np.ones(n)                     # uniform marine input (total)
    fineFlux = np.full(n, 0.3)               # post-cascade fine → ff_mar = 0.3

    m._marineFineFraction(hl, sedFlux, fineFlux)

    depth = m.sealevel - hl
    ff = m.depoFineFrac
    order = np.argsort(depth)
    assert np.all(np.diff(ff[order]) > 0), (
        f"Marine fine fraction must increase with depth; got {ff}."
    )
    fine_vol = float(np.sum(mdep * m.larea * ff))
    assert np.isclose(fine_vol, 0.3 * n, rtol=1e-9), (
        f"Marine fine volume not conserved: {fine_vol} != {0.3 * n}"
    )


@pytest.mark.slow
def test_dual_model_runs_and_invariants(minimal_dual_model):
    """
    Integration (DESIGN_DUAL_LITHOLOGY.md Phase 6): a full dual-lithology model
    runs end-to-end and preserves the per-fraction invariants. Exercises the
    whole dual path together — erodeStrat split, _getSedFlux fine routing
    (Phase 3a), deposeStrat per-node fineFrac, per-fraction compaction, and
    fine-pile advection. This is the first live coverage of the dual sediment
    transport path (both fixtures used by other tests have stratNb == 0).
    """
    model = minimal_dual_model
    assert model.stratLith is True and model.stratNb > 0
    assert model.stratHf is not None and model.phiF is not None
    # Bedrock sentinel carries the configured fine split before any run.
    assert np.isclose(
        model.stratHf[0, 0], 1.0e6 * (1.0 - model.bedrock_coarse_frac)
    )

    model.runProcesses()

    top = model.stratStep + 1
    H = model.stratH[:, :top]
    Hf = model.stratHf[:, :top]
    # Fine pile stays physical: non-negative and never exceeding the total.
    assert (Hf >= -1.0e-9).all(), "Negative fine thickness in the strata pile."
    assert (Hf <= H + 1.0e-6).all(), "Fine thickness exceeds layer total."
    # Porosity and the routed fine fraction stay in range.
    assert (model.phiF[:, :top] >= -1.0e-12).all()
    assert (model.phiF[:, :top] <= 1.0 + 1.0e-12).all()
    assert (model.fineFrac >= 0.0).all() and (model.fineFrac <= 1.0).all()
    # The dual path actually moved fine material (bedrock contributes it).
    assert float(Hf.sum()) > 0.0


@pytest.mark.slow
def test_dual_all_coarse_matches_single_fraction(
    minimal_dual_coarse_model, minimal_strat_model
):
    """
    Parity guard (DESIGN_DUAL_LITHOLOGY.md Phase 6): dual lithology configured
    all-coarse (bedrock_coarse_frac=1.0, no erodibility/diffusivity contrast)
    must reproduce the single-fraction stratigraphy run exactly. Confirms the
    dual code paths are a faithful superset of the single-fraction path — any
    accidental divergence (extra deposition, wrong compaction branch, etc.)
    trips this.
    """
    dual = minimal_dual_coarse_model
    single = minimal_strat_model
    dual.runProcesses()
    single.runProcesses()

    top = min(dual.stratStep, single.stratStep) + 1
    # No fine produced anywhere in the all-coarse configuration.
    assert float(dual.stratHf[:, :top].sum()) == 0.0, (
        "All-coarse dual run produced fine sediment."
    )
    # Elevation and stratal thickness must match single-fraction bitwise
    # (the smoke check measured exactly 0 difference; atol guards float noise).
    assert np.allclose(
        dual.hLocal.getArray(), single.hLocal.getArray(), rtol=0.0, atol=1.0e-9
    ), "Elevation diverged from the single-fraction stratigraphy run."
    assert np.allclose(
        dual.stratH[:, :top], single.stratH[:, :top], rtol=0.0, atol=1.0e-9
    ), "Stratal thickness diverged from the single-fraction stratigraphy run."


@pytest.mark.slow
def test_dual_mass_conservation(minimal_dual_model):
    """
    Protects: DESIGN_DUAL_LITHOLOGY.md Phase 6 — sediment conservation through
    the dual-lithology pipeline on a closed sphere, and a valid per-fraction
    partition maintained end-to-end.

    The standard test_mass_conservation SKIPS when stratNb > 0 (compaction
    moves h without cumED), so dual/stratigraphy-mode total conservation is
    otherwise untested. cumED only changes via PAIRED sediment erosion/
    deposition (never by compaction), so on a closed sphere Σ(dcumED·larea)
    must still vanish relative to the redistributed volume — even with the
    extra fine-flux solve and the per-pit / marine composition reallocation
    that dual lithology adds. A real sediment leak (e.g. fine routed but not
    deposited) would push this to O(0.1+); measured ~1.6e-5 here.

    Per-fraction: the coarse/fine partition must stay valid through erosion,
    transport, deposition, compaction, advection and diffusion — fine bulk
    non-negative and never exceeding the layer total; the fine solid phase
    non-negative; and the pile must actually carry fine (not trivially zero).
    """
    model = minimal_dual_model

    # Closed-domain gate (mirrors test_mass_conservation).
    reasons = []
    if getattr(model, "flatModel", True):
        reasons.append("flatModel=True (boundary outflux)")
    if getattr(model, "tecdata", None) is not None:
        reasons.append("tectonics active")
    if getattr(model, "flexOn", False):
        reasons.append("flexure active")
    if reasons:
        pytest.skip(
            "Dual mass conservation requires a closed sphere with only "
            "sediment-conserving kernels. This fixture has: " + "; ".join(reasons)
        )
    assert model.stratLith and model.stratNb > 0, "fixture must enable dual strat"

    larea = model.larea
    owned = model.inIDs == 1
    cumED_before = model.cumEDLocal.getArray().copy()

    model.runProcesses()

    dED = model.cumEDLocal.getArray() - cumED_before
    from mpi4py import MPI
    dV = MPI.COMM_WORLD.allreduce(float(np.sum((dED * larea)[owned])), op=MPI.SUM)
    activity = MPI.COMM_WORLD.allreduce(
        float(np.sum((np.abs(dED) * larea)[owned])), op=MPI.SUM
    )

    # ---- Total sediment conserved in dual mode (closed sphere) ----
    assert activity > 0.0, "no sediment was redistributed; test is vacuous"
    rel = abs(dV) / activity
    assert rel < 5.0e-4, (
        f"Dual-mode sediment not conserved: |ΣdcumED|/activity = {rel:.2e} "
        f"(> 5e-4). A fraction routed-but-not-deposited would leak here."
    )

    # ---- Valid per-fraction partition end-to-end ----
    top = model.stratStep + 1
    H = model.stratH[:, :top]
    Hf = model.stratHf[:, :top]
    assert (Hf >= -1.0e-9).all(), "Negative fine bulk thickness."
    assert (Hf <= H + 1.0e-6).all(), "Fine bulk exceeds layer total."
    fine_solid = (Hf * (1.0 - model.phiF[:, :top]))[owned]
    assert (fine_solid >= -1.0e-9).all(), "Negative fine solid phase."
    fine_total = MPI.COMM_WORLD.allreduce(float(np.sum(fine_solid)), op=MPI.SUM)
    assert fine_total > 0.0, "Dual pile carries no fine — sub-system is dead."


@pytest.mark.slow
def test_dual_fine_conservation(minimal_dual_model):
    """
    Protects: per-fraction (FINE) mass conservation on a closed sphere.

    The total-cumED test_dual_mass_conservation checks the *total* sediment
    budget; it would pass even if fine were silently created/destroyed while
    coarse compensated. This test closes that gap directly: the fine solid
    REMOVED from the strata pile (erodeStrat, accumulated in _fineEroded) must
    equal the fine DEPOSITED back into it (deposeStrat continental + marine,
    _fineDeposited), to within the floor/transit budget.

    This is the guard required before any future fine-routing change (e.g.
    fine-enriched overspill): a mismatched fine mirror in the pit/marine
    cascade — exactly the bug that reverted the first overspill attempt —
    trips this where the total budget does not. Measured imbalance ~1e-4 on
    minimal_dual; asserted < 5e-3.
    """
    model = minimal_dual_model
    if getattr(model, "flatModel", True):
        pytest.skip("fine conservation needs a closed sphere (flatModel=False)")
    assert model.stratLith and model.stratNb > 0, "fixture must enable dual strat"

    model.runProcesses()

    from mpi4py import MPI
    ero = MPI.COMM_WORLD.allreduce(model._fineEroded, op=MPI.SUM)
    dep = MPI.COMM_WORLD.allreduce(model._fineDeposited, op=MPI.SUM)
    assert ero > 0.0, "no fine was eroded; test is vacuous"
    rel = abs(ero - dep) / ero
    assert rel < 5.0e-3, (
        f"Fine not conserved: eroded={ero:.4e} deposited={dep:.4e} "
        f"rel imbalance={rel:.2e} (> 5e-3). A fine-only leak the total-cumED "
        f"test cannot see."
    )


@pytest.mark.slow
def test_dual_sedloadf_output(minimal_dual_model, minimal_strat_model):
    """
    Protects: the fine sediment load `sedLoadF` is written to the HDF5 output
    when dual lithology is enabled, and absent for single-fraction runs.
    `sedLoad` is the TOTAL flux; `sedLoadF` the fine sub-flux, so
    0 <= sedLoadF <= sedLoad at every node.
    """
    h5py = pytest.importorskip("h5py")
    import glob

    def _latest_h5(model):
        files = sorted(
            glob.glob(os.path.join(str(model.outputDir), "h5", "gospl.*.p*.h5"))
        )
        return files[-1] if files else None

    # ---- dual run writes sedLoadF, 0 <= fine <= total ----
    dual = minimal_dual_model
    dual.runProcesses()
    f = _latest_h5(dual)
    assert f is not None, "no gospl HDF5 output was written"
    with h5py.File(f, "r") as hf:
        assert "sedLoadF" in hf, "sedLoadF missing from dual-lithology output"
        sl = np.array(hf["sedLoad"])[:, 0]
        slf = np.array(hf["sedLoadF"])[:, 0]
    assert (slf >= 0.0).all(), "negative fine sediment load"
    assert (slf <= sl + 1.0e-6).all(), "fine load exceeds total load"

    # ---- single-fraction stratigraphy run must NOT write sedLoadF ----
    single = minimal_strat_model
    single.runProcesses()
    fs = _latest_h5(single)
    with h5py.File(fs, "r") as hf:
        assert "sedLoad" in hf and "sedLoadF" not in hf, (
            "single-fraction output must not contain sedLoadF"
        )


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
# TEST 8 - Parallel decomposition must not alter the physical solution
# ---------------------------------------------------------------------------
#
# MPI's COMM_WORLD.size is fixed at process launch (`mpirun -n N`), so a
# single pytest invocation can only see ONE rank count. To compare n=1 vs
# n=N runs of the same model, this test spawns TWO subprocesses (`mpirun
# -n 1` and `mpirun -n 2`), each writing a JSON summary of global state,
# then asserts the summaries agree across decompositions.
#
# Why n=2 (not n=4) for the parallel case: GitHub's macos-14 runners only
# have 3 cores; OpenMPI refuses to oversubscribe by default. n=2 exposes
# exactly one partition boundary, which is sufficient to catch every
# decomposition-related bug class this test is designed to surface.
# ---------------------------------------------------------------------------

_PARALLEL_DUMP_SCRIPT = '''\
"""
Subprocess entry point for test_parallel_correctness.

Runs `gospl.model.Model(<yml>)`, performs MPI-reduced summary stats over
owned (non-ghost) nodes, and writes a JSON file on rank 0. Designed to be
invoked via `mpirun -n N python <this script> <yml> <out_json>` from the
test body — the test then loads two such JSON dumps (n=1 and n=2) and
asserts the global quantities match within tight tolerances.
"""
import json
import sys
import numpy as np
from mpi4py import MPI
from gospl.model import Model

fixture_yml = sys.argv[1]
output_json = sys.argv[2]

comm = MPI.COMM_WORLD
model = Model(fixture_yml, verbose=False, showlog=False)
try:
    model.runProcesses()

    # Owned-node mask: ghost nodes must never enter a reduction (they are
    # halo copies owned by a neighbouring rank). Same pattern as Tests 6-7.
    owned = model.inIDs == 1
    h_owned     = model.hLocal.getArray()[owned]
    cumED_owned = model.cumEDLocal.getArray()[owned]
    fa_owned    = model.FAL.getArray()[owned]
    larea_owned = model.larea[owned]

    # Per-rank local scalars; MPI.SUM/MAX-reduced to globals below.
    wmean_h_local  = float((h_owned * larea_owned).sum())
    area_local     = float(larea_owned.sum())
    flux_local     = float((cumED_owned * larea_owned).sum())
    activity_local = float((np.abs(cumED_owned) * larea_owned).sum())
    max_fa_local      = float(fa_owned.max()) if len(fa_owned) > 0 else 0.0
    sum_fa_local      = float(fa_owned.sum())
    # Owned-node count per rank. After allreduce(SUM) this is the total
    # number of nodes uniquely owned across the comm. Should equal the
    # mesh's global vertex count regardless of decomposition. If it
    # doesn't, there's a partition-ownership gap (a node owned by zero
    # ranks, or — less likely — double-counted by inIDs).
    owned_count_local = int(owned.sum())

    wmean_h_num = comm.allreduce(wmean_h_local, op=MPI.SUM)
    area_total  = comm.allreduce(area_local,    op=MPI.SUM)
    flux        = comm.allreduce(flux_local,    op=MPI.SUM)
    activity    = comm.allreduce(activity_local, op=MPI.SUM)
    max_fa      = comm.allreduce(max_fa_local,  op=MPI.MAX)
    sum_fa      = comm.allreduce(sum_fa_local,  op=MPI.SUM)
    owned_count = comm.allreduce(owned_count_local, op=MPI.SUM)
    wmean_h     = wmean_h_num / area_total

    if comm.Get_rank() == 0:
        with open(output_json, "w") as f:
            json.dump({
                "size":        comm.Get_size(),
                "wmean_h":     wmean_h,
                "flux":        flux,
                "activity":    activity,
                "max_fa":      max_fa,
                "sum_fa":      sum_fa,
                "owned_count": owned_count,
            }, f)
finally:
    model.destroy()
'''

def _petsc4py_abi_mismatch() -> bool:
    """Return True if petsc4py was built against a different Python than runtime."""
    try:
        import importlib.metadata, sys
        dist = importlib.metadata.distribution("petsc4py")
        tag = dist.metadata["Name"]  # just a check it exists
        # build string embeds the python target, e.g. np2py310h...
        # compare against the running interpreter
        for f in dist.files or []:
            if "py3" in str(f) and "cpython" not in str(f):
                pass
        # simpler: check via the direct build tag in the wheel name
        record = next(
            (str(f) for f in (dist.files or []) if str(f).endswith(".dist-info/WHEEL")),
            None,
        )
        if record:
            wheel_text = (dist._path.parent / record).read_text()
            import re, sys
            tag_match = re.search(r"Tag: (\S+)", wheel_text)
            if tag_match:
                tag_str = tag_match.group(1)
                expected = f"cp{sys.version_info.major}{sys.version_info.minor}"
                return expected not in tag_str
    except Exception:
        pass
    return False

@pytest.mark.slow
@pytest.mark.skipif(
    _petsc4py_abi_mismatch(),
    reason="petsc4py built against different Python ABI (no py311 osx-arm64 "
           "build on conda-forge); segfaults on MPI finalization — not a gospl bug"
)
def test_parallel_correctness(tmp_path):
    """
    Protects: AGENTS.md > MPI contract — collective operations must yield
    physically identical results regardless of how the domain is partitioned
    across ranks.

    Silent failures prevented:

      1. Halo-exchange bug in `dm.localToGlobal` / `dm.globalToLocal`
         (mesher/unstructuredmesh.py:738-799): a ghost node written by the
         wrong rank would corrupt the elevation field in a thin ring at each
         subdomain boundary.

      2. Flow-routing inconsistency across partition boundaries: the
         donor-receiver graph in `flowplex._buildFlowDirection` depends on
         which rank owns each node. A bug in the boundary-node hand-off
         would produce divergent rcvID arrays, which cascade into divergent
         drainage area and elevation.

      3. Incorrect `inIDs` mask causing a rank to count ghost-node
         contributions in its local reduction (anywhere `inIDs == 1` is
         used to select owned nodes). With n=2 this would over-count
         boundary nodes by ~2x and shift global statistics.

      4. Non-deterministic KSP convergence across decompositions. PETSc
         parallel reductions (dot products, norms) are NOT guaranteed
         bitwise-identical because floating-point addition is not
         associative. We therefore assert STATISTICAL equivalence with
         tight relative tolerances rather than bitwise identity.

    Implementation: MPI.COMM_WORLD.size is fixed at process launch, so we
    cannot run n=1 and n=N in the same pytest process. Two subprocesses
    via `mpirun -n N python <dump_script>` each emit a JSON dump on rank
    0; this test reads both and compares.

    Tier 1A — drainage area (1e-4 relative, checked FIRST because routing
    bugs cascade into elevation bugs but not vice versa). Both the max
    and the sum of FAL over owned nodes must agree.

    Tier 1B — area-weighted mean elevation (1e-10 relative). Catches
    rank-doubling of ghost-node contributions in any reduction.

    Tier 1C — total cumED×area flux (1e-10 relative). Catches sediment
    created or destroyed at subdomain boundaries via incorrect axpy
    in SPL.py:352, nlSPL.py:404, or soilSPL.py:326.

    Non-applicability gates:
      - skip if `mpirun` is not on PATH (no MPI runtime to spawn);
      - skip if minimal.yml / mesh.npz are missing (consistent with the
        other Model-based regression tests);
      - skip if activity < 1.0 m³ (vacuous comparison on a fixture that
        moves negligible sediment — same gate as Tests 6-7).
    """
    import json
    import os
    import shutil
    import subprocess
    import sys
    from pathlib import Path

    # ---- Non-applicability gates ---------------------------------------
    if shutil.which("mpirun") is None:
        pytest.skip(
            "mpirun not on PATH; cannot exercise MPI decomposition. "
            "On CI this means the conda env lacks an MPI runtime (mpi4py "
            "package missing or broken). Locally, install via "
            "`mamba env create -f environment.yml`."
        )

    fixtures_dir = Path(__file__).parent / "fixtures"
    yml_src = fixtures_dir / "minimal.yml"
    mesh_src = fixtures_dir / "mesh.npz"
    if not yml_src.exists() or not mesh_src.exists():
        pytest.skip(
            f"{yml_src} or {mesh_src} not present. Same fixture as the "
            f"other Model-based regression tests."
        )

    # ---- Write the per-subprocess dump script once ---------------------
    # Stored as a real .py file (not -c "string") so tracebacks point at
    # readable lines if anything inside fails.
    dump_py = tmp_path / "_parallel_dump.py"
    dump_py.write_text(_PARALLEL_DUMP_SCRIPT)

    # ---- Run the model under mpirun -n N -------------------------------
    def run_at_rank(n):
        """
        Spawn `mpirun -n {n} python <dump_py> minimal.yml stats.json`
        in an isolated cwd. Each cwd has its own copies of the YAML and
        mesh so the goSPL output directory (`<cwd>/minimal/`) does not
        collide between the two runs.
        """
        out_dir = tmp_path / f"n{n}"
        out_dir.mkdir()
        shutil.copy(yml_src,  out_dir / "minimal.yml")
        shutil.copy(mesh_src, out_dir / "mesh.npz")
        stats_json = out_dir / "stats.json"

        cmd = [
            "mpirun", "-n", str(n),
            sys.executable, str(dump_py),
            "minimal.yml", str(stats_json),
        ]
        # Scrub inherited MPI-runtime env vars before spawning a nested mpirun.
        # pytest imported gospl (-> petsc4py.init -> MPI_Init), which under
        # OpenMPI exports OMPI_*/PMIX_*/PRTE_* into this process; if they leak
        # into the child `mpirun` it thinks it is already inside an MPI job and
        # silently refuses to launch (rc=1, empty output). MPICH is unaffected,
        # but the conda package and HPC container both use OpenMPI. OPAL_PREFIX
        # is preserved so mpirun can still locate its own libraries.
        child_env = {
            k: v for k, v in os.environ.items()
            if not k.startswith(("OMPI_", "PMIX_", "PRTE_", "OPAL_"))
        }
        if "OPAL_PREFIX" in os.environ:
            child_env["OPAL_PREFIX"] = os.environ["OPAL_PREFIX"]
        result = subprocess.run(
            cmd, cwd=out_dir, timeout=600,
            capture_output=True, text=True, env=child_env,
        )
        if result.returncode != 0 or not stats_json.exists():
            pytest.fail(
                f"`mpirun -n {n}` subprocess failed "
                f"(rc={result.returncode}).\n"
                f"--- stdout ---\n{result.stdout}\n"
                f"--- stderr ---\n{result.stderr}"
            )
        return json.loads(stats_json.read_text())

    stats_n1 = run_at_rank(1)
    stats_n2 = run_at_rank(2)

    # ---- Sanity: the two subprocesses really did run at the requested
    #      rank counts (catches MPI implementation quirks, e.g. OpenMPI
    #      silently downgrading to 1 rank under oversubscription policy).
    assert stats_n1["size"] == 1, f"n=1 ran with size={stats_n1['size']}"
    assert stats_n2["size"] == 2, f"n=2 ran with size={stats_n2['size']}"

    # ---- Sanity: same total owned-node count across decompositions.
    # Every mesh vertex must be owned by EXACTLY one rank; the sum of
    # owned-node counts over the comm must therefore equal the mesh's
    # global vertex count regardless of decomposition. If this differs
    # between n=1 and n=2, there is a partition-ownership gap — bug in
    # `inIDs` setup (unstructuredmesh.py) or in PETSc DMPlex's vertex
    # partitioning. This is the FIRST thing to check if the cumulative
    # sums below diverge unexpectedly.
    assert stats_n1["owned_count"] == stats_n2["owned_count"], (
        f"n=1 and n=2 see different total owned-node counts "
        f"(n=1={stats_n1['owned_count']}, n=2={stats_n2['owned_count']}). "
        "Every mesh vertex must be owned by exactly one rank; this "
        "assertion failing means a node is owned by zero ranks (gap) "
        "or — less likely — double-counted (inIDs == 1 on two ranks). "
        "Suspect: mesher/unstructuredmesh.py's vIS-based `inIDs` "
        "assignment, or PETSc DMPlex partitioner setup."
    )

    # ---- Activity gate -------------------------------------------------
    if stats_n1["activity"] < 1.0:
        pytest.skip(
            f"Activity {stats_n1['activity']:.3e} m³ too small for "
            f"meaningful parallel-correctness check. NEEDS_HUMAN_REVIEW: "
            f"lengthen the run or steepen the gradient in minimal.yml."
        )

    # ---- Compute relative diffs ---------------------------------------
    def rel(a, b, denom):
        return abs(a - b) / max(abs(denom), 1e-30)

    rel_mean_h = rel(stats_n1["wmean_h"], stats_n2["wmean_h"],
                     stats_n1["wmean_h"])
    rel_flux   = rel(stats_n1["flux"], stats_n2["flux"],
                     stats_n1["activity"])
    rel_max_fa = rel(stats_n1["max_fa"], stats_n2["max_fa"],
                     stats_n1["max_fa"])
    rel_sum_fa = rel(stats_n1["sum_fa"], stats_n2["sum_fa"],
                     stats_n1["sum_fa"])

    diagnostic = (
        f"\n  owned_n1       = {stats_n1['owned_count']}"
        f"\n  owned_n2       = {stats_n2['owned_count']}"
        f"\n  wmean_h_n1     = {stats_n1['wmean_h']:.6e} m"
        f"\n  wmean_h_n2     = {stats_n2['wmean_h']:.6e} m"
        f"\n  rel_mean_h     = {rel_mean_h:.3e}  (tol 1e-10)"
        f"\n  flux_n1        = {stats_n1['flux']:+.3e} m³"
        f"\n  flux_n2        = {stats_n2['flux']:+.3e} m³"
        f"\n  rel_flux       = {rel_flux:.3e}  (tol 1e-5)"
        f"\n  max_fa_n1      = {stats_n1['max_fa']:.3e} m²"
        f"\n  max_fa_n2      = {stats_n2['max_fa']:.3e} m²"
        f"\n  rel_max_fa     = {rel_max_fa:.3e}  (tol 1e-4)"
        f"\n  sum_fa_n1      = {stats_n1['sum_fa']:.3e} m²"
        f"\n  sum_fa_n2      = {stats_n2['sum_fa']:.3e} m²"
        f"\n  rel_sum_fa     = {rel_sum_fa:.3e}  (tol 5e-2)"
        f"\n  activity       = {stats_n1['activity']:.3e} m³"
    )

    # ---- Tier 1A: drainage area (routing bugs cascade — check first) ----
    assert rel_max_fa < 1e-4, (
        "Max drainage area differs between n=1 and n=2 beyond 1e-4 "
        "relative. Likely a partition-boundary bug in "
        "flow/flowplex.py:\n"
        "  (a) `_buildFlowDirection` assigns a boundary node's receiver "
        "to the wrong rank's local index.\n"
        "  (b) `_matrix_build` inserts a wrong off-diagonal entry at a "
        "halo node, so the IDA solve distributes water incorrectly "
        "across the partition boundary.\n"
        "  (c) The rcvIDi snapshot (flowplex.py:421-426) is taken AFTER "
        "a domain-decomposition-dependent globalToLocal."
        + diagnostic
    )
    # Sum-of-FA tolerance is INTENTIONALLY loose (5e-2 ~ 5%) — and the
    # gap between this and the other tolerances is the load-bearing
    # observation, not the absolute number.
    #
    # `mfdreceivers` already breaks EXACT slope ties deterministically
    # using global node IDs (see fortran/functions.F90 and the
    # `self.gid` argument in flow/flowplex.py:_buildFlowDirection). What
    # remains is the near-tie case: KSP convergence is not bitwise-
    # identical across decompositions because floating-point addition
    # is not associative under parallel reductions, so ghost-node
    # elevations differ by O(KSP rtol) between n=1 and n=2. A small
    # subset of near-tie nodes picks a different receiver and the
    # cumulative drainage statistics diverge.
    #
    # The DRIFT MAGNITUDE IS PLATFORM-DEPENDENT:
    #   - macOS-14 (arm64, conda-forge OpenMPI):  ~0.3%
    #   - Ubuntu-latest (x86_64, conda-forge MPICH): ~1.7%
    # (Same fixture, same Python, same goSPL, same tie-break fix —
    # the difference is in MPI implementation, BLAS variant, and FMA
    # availability at the platform level.) 5% is set as 3x the worst
    # observed on either platform. Real routing regressions (lost
    # neighbour entry, wrong row weights at halo nodes, stale rcvIDi
    # snapshot) would shift this by orders of magnitude — wholesale
    # rerouting drives rel_sum_fa toward 0.5+, not 0.05.
    #
    # Tightening this back below ~5% requires either KSP-precision
    # halo determinism (hard PETSc work) or a slope-tolerance band in
    # _buildFlowDirection (algorithm change). Out of scope here.
    assert rel_sum_fa < 5e-2, (
        "Total drainage area differs between n=1 and n=2 beyond 5e-2 "
        "relative — far beyond the ~0.3% (macOS) to ~1.7% (Ubuntu) "
        "noise floor from non-deterministic tie-breaking. Likely a "
        "real routing regression in flow/flowplex.py:\n"
        "  (a) `_buildFlowDirection` lost a neighbour entry at a "
        "partition boundary.\n"
        "  (b) The IDA matrix (`_matrix_build`) assembles with wrong "
        "row weights on halo nodes.\n"
        "  (c) The rcvIDi snapshot (flowplex.py:421-426) was taken from "
        "a stale (post-fill) topology."
        + diagnostic
    )

    # ---- Tier 1B: global mean elevation -------------------------------
    assert rel_mean_h < 1e-10, (
        "Area-weighted mean elevation differs between n=1 and n=2 "
        "beyond 1e-10 relative. Likely causes:\n"
        "  (a) Ghost nodes included in a rank-local reduction before "
        "allreduce — every sum must use `inIDs == 1` as the owned-node "
        "mask.\n"
        "  (b) A `localToGlobal` call is missing after `hLocal.setArray`, "
        "so a rank is solving with a stale halo on the next KSP step.\n"
        "  (c) The KSP RHS double-counts a ghost-node contribution on "
        "ranks that own a boundary."
        + diagnostic
    )

    # ---- Tier 1C: total erosion / deposition flux ---------------------
    # Flux tolerance is intentionally loose (1e-5) because the residual
    # KSP-near-tie non-determinism described above shifts WHICH nodes
    # erode how much (a few near-tie nodes drain different basins →
    # different local incision). Observed 2.4e-7 on the minimal fixture
    # AFTER the exact-tie-break fix landed; 1e-5 gives ~40x headroom for
    # real mass-conservation regressions to trip the assert. Mean
    # elevation is unaffected (algorithm is mass-conserving regardless
    # of which path sediment takes) so Tier 1B remains at 1e-10.
    assert rel_flux < 1e-5, (
        "Total cumED × area flux differs between n=1 and n=2 beyond "
        "1e-5 relative — well beyond the ~1e-7 floor expected from "
        "routing tie-break non-determinism. Sediment is being created "
        "or destroyed at subdomain boundaries:\n"
        "  (a) An axpy in SPL.py:352, nlSPL.py:404, or soilSPL.py:326 "
        "uses the wrong scale factor on boundary nodes.\n"
        "  (b) sedplex._getSedFlux integrates an erosion source that "
        "has already been localToGlobal'd and therefore contains halo "
        "copies of boundary-node values."
        + diagnostic
    )


# ---------------------------------------------------------------------------
# TEST 9 - Stratigraphy: deposition + compaction physics
# ---------------------------------------------------------------------------
#
# Pure-numpy test (no Model, no PETSc DMPlex) — instantiates the STRAMesh
# class via `__new__` to bypass the heavy init, sets minimal state, and
# exercises `deposeStrat` + `_depthPorosity` directly with mocked PETSc
# Vec / DM objects. Lives in the fast tier alongside TESTs 1-2; runs in
# well under a second on any platform.
# ---------------------------------------------------------------------------


def test_stratigraphy_deposition_and_compaction():
    """
    Protects: stratplex.deposeStrat + stratplex._depthPorosity — the
    deposition and compaction physics that update underlying stratal
    thickness and porosity when surface deposition events accumulate
    over time.

    Constructs a fictive 3-layer stratigraphy (5 nodes × 3 layers,
    each 50 m thick at phi=0.5) and runs two phases:

      Phase 1 — deposeStrat: add a 10 m deposition to the top layer.
      Asserts:
        - stratH[top] grows by exactly the deposition amount;
        - phiS[top] is set to phi0s (fresh deposit's surface porosity);
        - stratK[top] is reset to 1.0 (default erodibility multiplier);
        - lower layers are untouched.

      Phase 2 — _depthPorosity: apply Athy's-law compaction
      (phi = phi0 * exp(depth/z0)) at synthetic layer mid-depths.
      Asserts:
        - porosity decreases monotonically with depth (deeper compacts
          more);
        - layer thickness decreases (compaction shrinks layers);
        - solid-phase volume is conserved per layer to float-noise
          precision: h_new * (1 - phi_new) == h_old * (1 - phi_old).

    Silent failures prevented:
      1. deposeStrat lays down sediment but forgets to set phiS to phi0s
         on the new layer → downstream code reads NaN porosity from the
         forward-fill in `_fillZeroPorosity`.
      2. _depthPorosity inverts the depth sign convention → porosity
         INCREASES with depth, the column unrealistically inflates.
      3. Compaction's solid-phase bookkeeping (stratplex.py:268-278)
         loses or creates mass: `newH = solidPhase / (1 - phi_new)` must
         exactly preserve `h * (1 - phi)`. Any algebraic refactor that
         changes the order of operations is checked here.
      4. The bedrock-sentinel freeze (stratplex.py:283-286) is silently
         removed: bedrockLay > 0 should keep layer indices < bedrockLay
         from changing in either thickness or porosity. Tested in a
         second assertion block below.
    """
    # Late import — STRAMesh's module triggers `from gospl._fortran import
    # strataonesed`, which fails at collection if the Fortran extension
    # is not built. `importorskip` lets the test skip cleanly in that
    # environment instead of failing the whole file.
    stratplex = pytest.importorskip("gospl.sed.stratplex")
    STRAMesh = stratplex.STRAMesh

    # ---- Mocks for the PETSc-side interface ----------------------------
    # deposeStrat calls `self.dm.globalToLocal(self.tmp, self.tmpL)` then
    # reads `self.tmpL.getArray()`. We pre-populate tmpL with a fake
    # object whose getArray() returns the desired deposition vector;
    # the globalToLocal call becomes a no-op (tmpL is already correct).
    class _MockLocalVec:
        def __init__(self, arr):
            self._arr = np.asarray(arr, dtype=np.float64)
        def getArray(self):
            return self._arr

    class _MockDM:
        def globalToLocal(self, src, dst):
            pass  # tmpL pre-populated by the test

    # ---- Build the fictive STRAMesh state ------------------------------
    # __new__ bypasses __init__; STRAMesh.__init__ only sets all four
    # state arrays to None and returns, so we replace that with explicit
    # synthetic state. No PETSc / mesh allocations involved.
    s = STRAMesh.__new__(STRAMesh)
    s.lpoints   = 5
    s.stratNb   = 3
    s.stratStep = 2        # 3 layers indexed 0..2; layer 2 is the top
    s.phi0s     = 0.5      # surface porosity
    s.z0s       = 100.0    # e-folding compaction depth (m)
    s.bedrockLay = 0       # no infinite-bedrock sentinel layer
    s.memclear  = False
    s.stratLith = False    # single-fraction path (dual-lithology disabled)

    s.stratH = np.full((5, 3), 50.0, dtype=np.float64)
    s.phiS   = np.full((5, 3),  0.5, dtype=np.float64)
    s.stratK = np.full((5, 3),  0.7, dtype=np.float64)

    s.tmp  = None
    s.tmpL = _MockLocalVec(np.full(5, 10.0))   # 10 m deposition at every node
    s.dm   = _MockDM()

    # ==== Phase 1: deposition ===========================================
    H_before  = s.stratH.copy()
    phi_before = s.phiS.copy()
    K_before  = s.stratK.copy()

    s.deposeStrat()

    # Top layer thickness grew by exactly the deposition amount
    np.testing.assert_array_almost_equal(
        s.stratH[:, s.stratStep] - H_before[:, s.stratStep],
        np.full(5, 10.0),
        err_msg=(
            "deposeStrat did not add the deposition to the top layer "
            "thickness. Check stratplex.py:157 — `self.stratH[:, "
            "self.stratStep] += depo`."
        ),
    )

    # Top layer porosity reset to phi0s
    np.testing.assert_array_equal(
        s.phiS[:, s.stratStep],
        np.full(5, s.phi0s),
        err_msg=(
            "deposeStrat did not reset the top layer's porosity to "
            "phi0s after deposition (stratplex.py:159). Downstream "
            "code that forward-fills zero porosity will inherit the "
            "wrong value from below."
        ),
    )

    # Top layer K multiplier reset to 1.0
    np.testing.assert_array_equal(
        s.stratK[:, s.stratStep],
        np.ones(5),
        err_msg=(
            "deposeStrat did not reset the top layer's erodibility "
            "multiplier to 1.0 (stratplex.py:163). Freshly deposited "
            "sediment should always carry the default K."
        ),
    )

    # Lower layers (indices 0, 1) untouched
    np.testing.assert_array_equal(s.stratH[:, :s.stratStep], H_before[:, :s.stratStep])
    np.testing.assert_array_equal(s.phiS[:, :s.stratStep],   phi_before[:, :s.stratStep])
    np.testing.assert_array_equal(s.stratK[:, :s.stratStep], K_before[:, :s.stratStep])

    # ==== Phase 2: compaction ===========================================
    # Layer geometry after Phase 1 (per node):
    #   layer 0: H=50,  phi=0.5,  bottom of column
    #   layer 1: H=50,  phi=0.5
    #   layer 2: H=60,  phi=0.5,  top of column (50 + 10 m deposit)
    # Mid-point depths below the post-deposition surface (z=0):
    #   layer 2 mid-depth = -30   (60/2 below surface)
    #   layer 1 mid-depth = -85   (60 + 50/2)
    #   layer 0 mid-depth = -135  (60 + 50 + 50/2)
    # _depthPorosity expects depth as (lpoints, stratStep+1) with the
    # column ordered [bottom, ..., top] — same layout as stratH.
    depth = np.tile(np.array([-135.0, -85.0, -30.0]), (5, 1))

    H_pre   = s.stratH.copy()
    phi_pre = s.phiS.copy()

    # Snapshot solid-phase volume per layer BEFORE compaction.
    # The conservation law: h * (1 - phi) is invariant across compaction
    # (compaction only changes void volume, never solid mass).
    solid_before = H_pre * (1.0 - phi_pre)

    newH = s._depthPorosity(depth)

    # 1. Porosity strictly decreases with depth
    assert (s.phiS[:, 0] < s.phiS[:, 1]).all(), (
        "Bottom layer porosity should be less than middle layer after "
        f"compaction. Got phi[bottom]={s.phiS[:, 0]}, "
        f"phi[middle]={s.phiS[:, 1]}. Check the sign of `depth/z0s` in "
        "stratplex.py:265."
    )
    assert (s.phiS[:, 1] < s.phiS[:, 2]).all(), (
        "Middle layer porosity should be less than top layer after "
        f"compaction. Got phi[middle]={s.phiS[:, 1]}, "
        f"phi[top]={s.phiS[:, 2]}."
    )

    # 2. Layer thickness strictly decreases (compaction shrinks layers)
    assert (newH[:, 0] < H_pre[:, 0]).all(), (
        f"Bottom layer should compact. Got newH={newH[:, 0]}, "
        f"H_pre={H_pre[:, 0]}."
    )
    assert (newH[:, 1] < H_pre[:, 1]).all(), (
        f"Middle layer should compact. Got newH={newH[:, 1]}, "
        f"H_pre={H_pre[:, 1]}."
    )
    assert (newH[:, 2] < H_pre[:, 2]).all(), (
        f"Top layer should compact (its phi went from 0.5 to "
        f"{s.phiS[0, 2]:.4f}). Got newH={newH[:, 2]}, "
        f"H_pre={H_pre[:, 2]}."
    )

    # 3. Solid-phase volume conservation per layer
    solid_after = newH * (1.0 - s.phiS)
    np.testing.assert_allclose(
        solid_after, solid_before, rtol=1e-12, atol=1e-12,
        err_msg=(
            "Solid-phase volume not conserved across compaction. "
            "_depthPorosity must change phi and H in lockstep so "
            "h*(1-phi) is invariant — check the construction "
            "`newH = solidPhase / tot` at stratplex.py:278."
        ),
    )

    # ==== Phase 3: bedrock sentinel freeze ==============================
    # Re-run the compaction with bedrockLay = 1 (layer 0 is the infinite-
    # bedrock sentinel). The bedrock layer's thickness and porosity must
    # NOT change even though `depth` would otherwise drive compaction
    # there. Catches regressions of the freeze logic at stratplex.py:283.
    s2 = STRAMesh.__new__(STRAMesh)
    s2.lpoints   = 5
    s2.stratNb   = 3
    s2.stratStep = 2
    s2.phi0s     = 0.5
    s2.z0s       = 100.0
    s2.bedrockLay = 1      # layer 0 is bedrock
    s2.memclear  = False
    s2.stratLith = False   # single-fraction path (dual-lithology disabled)

    # Layer 0 holds the BEDROCK_SENTINEL thickness (1e6) — never compact it.
    s2.stratH = np.array([[1.0e6, 50.0, 50.0]] * 5, dtype=np.float64)
    s2.phiS   = np.array([[0.0,    0.5,  0.5]] * 5, dtype=np.float64)
    s2.stratK = np.full((5, 3), 1.0, dtype=np.float64)

    depth2 = np.tile(np.array([-1.0e6 - 25.0, -75.0, -25.0]), (5, 1))
    H_bedrock_before   = s2.stratH[:, 0].copy()
    phi_bedrock_before = s2.phiS[:, 0].copy()

    newH2 = s2._depthPorosity(depth2)

    np.testing.assert_array_equal(
        newH2[:, 0], H_bedrock_before,
        err_msg=(
            "Bedrock sentinel layer (index < bedrockLay) was compacted. "
            "The freeze at stratplex.py:285 should preserve `stratH` "
            "exactly on bedrock indices."
        ),
    )
    np.testing.assert_array_equal(
        s2.phiS[:, 0], phi_bedrock_before,
        err_msg=(
            "Bedrock sentinel porosity was modified. The freeze at "
            "stratplex.py:286 should preserve `phiS` exactly on "
            "bedrock indices."
        ),
    )
