"""
Input-parser regression tests (no mesh, no PETSc DMPlex).

Protects: Forcing DataFrame layout contract and YAML parsing helpers (AGENTS.md).

Split out of the former tests/test_regression.py. Run this group alone with
`pytest -m parser`; see tests/README.md for the marker list.
"""

from __future__ import annotations

import pandas as pd
import pytest

# Skip the whole module (rather than erroring at collection) when the
# goSPL runtime stack is not installed.
inputparser = pytest.importorskip(
    "gospl.tools.inputparser",
    reason="goSPL runtime deps (petsc4py / ruamel.yaml / scipy) not installed",
)

pytestmark = [pytest.mark.parser]


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
