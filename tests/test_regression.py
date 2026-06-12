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
            "fine": {"phi0": 0.65, "z0": 1500.0},
            "bedrock_coarse_frac": 0.7,
            "fine_efficiency": 0.3,
            "pitInletBias": {"coarse": 0.8, "fine": 0.1},
            "Dc": 0.01,
            "Df": 0.05,
        }
    }
    parser._extraStrata()
    assert parser.stratLith is True
    assert parser.phi0c == 0.45 and parser.z0c == 3000.0
    assert parser.phi0f == 0.65 and parser.z0f == 1500.0
    assert parser.bedrock_coarse_frac == 0.7
    assert parser.fine_efficiency == 0.3
    assert parser.pit_inlet_bias_coarse == 0.8
    assert parser.pit_inlet_bias_fine == 0.1
    assert parser.Dc == 0.01 and parser.Df == 0.05

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

    def _strata_mesh(stratLith):
        m = stratplex.STRAMesh.__new__(stratplex.STRAMesh)
        m.strataFile = None
        m.lpoints = 8
        m.stratNb = 3
        m.phi0s = 0.49
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

    # Dual: fine fields allocated with the same shape as stratH, all-coarse.
    m = _strata_mesh(stratLith=True)
    m.readStratLayers()
    assert m.stratHf is not None and m.phiF is not None
    assert m.stratHf.shape == m.stratH.shape == (m.lpoints, m.stratNb)
    assert m.phiF.shape == (m.lpoints, m.stratNb)
    assert (m.stratHf == 0.0).all(), "Initial dual column should be all-coarse."


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
