"""
The per-step JSON-Lines run summary (gospl/tools/runsummary.py,
`Model(..., summary=path)` / `gospl --summary path`).

Agents and scripts read this file instead of parsing verbose logs, so its
schema is an interface: the keys asserted here are the ones documented in the
module docstring and in AGENTS.md > Working loop.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys

import pytest

from _helpers import FIXTURES_DIR, MPI_TIMEOUT, mpi_child_env

pytest.importorskip(
    "gospl.tools.inputparser",
    reason="goSPL runtime deps (petsc4py / ruamel.yaml / scipy) not installed",
)

pytestmark = pytest.mark.tools

STEP_KEYS = {"step", "t", "dt", "wall_s", "phases", "elev", "FA_max", "sealevel",
             "nonfinite", "volume", "ksp", "events"}


def _read(path):
    return [json.loads(line) for line in open(path)]


def _run_serial(tmp_path, summary):
    from gospl.model import Model

    cwd = os.getcwd()
    os.chdir(FIXTURES_DIR)
    try:
        m = Model("minimal.yml", verbose=False, showlog=False, summary=summary)
        try:
            m.runProcesses()
        finally:
            m.destroy()
    finally:
        os.chdir(cwd)
    return m


def test_summary_records_every_step(tmp_path):
    out = tmp_path / "run.jsonl"
    _run_serial(tmp_path, str(out))
    recs = _read(out)
    start, steps, end = recs[0], recs[1:-1], recs[-1]

    assert start["start"] is True and start["nranks"] == 1 and start["mpoints"] > 0
    assert end["end"] is True and end["steps"] == len(steps)
    n_expected = round((start["tEnd"] - start["tStart"]) / start["dt"])
    assert len(steps) == n_expected

    for i, r in enumerate(steps):
        assert STEP_KEYS <= set(r), STEP_KEYS - set(r)
        assert r["step"] == i
        assert r["t"] == pytest.approx(start["tStart"] + i * start["dt"])
        assert r["elev"]["min"] <= r["elev"]["mean"] <= r["elev"]["max"]
        assert r["nonfinite"] == {"elev": 0, "cumED": 0, "FA": 0}
        assert r["volume"]["eroded"] <= 0.0 <= r["volume"]["deposited"]
        assert r["volume"]["net"] == pytest.approx(
            r["volume"]["eroded"] + r["volume"]["deposited"])
        # the main discharge solve runs every step and is counted
        assert r["ksp"]["fatal"]["solves"] >= 1
        assert {e["kind"] for e in r["events"]} <= {
            "flow_cascade", "flow_ksp_fallback_failed", "flow_ksp_exact_rescue",
            "ksp_stale_guess_reset",
            "soil_snes_fallback",
            "soil_substep", "nlspl_snes_fallback", "nlspl_snes_failed"}

    # minimal.yml is a closed sphere: the per-step budget closes to the
    # documented floor (DEPOSIT_FLOOR + pit residue, ~1e-4 of activity).
    for r in steps:
        activity = -r["volume"]["eroded"] + r["volume"]["deposited"]
        if activity > 0:
            assert abs(r["volume"]["net"]) / activity < 1e-3


def test_summary_off_adds_nothing(tmp_path):
    m = _run_serial(tmp_path, None)
    assert m._summary is None
    assert getattr(m, "_summaryEvents", None) is None


@pytest.mark.mpi
def test_summary_parallel_matches_serial(tmp_path):
    """The step_end reductions are collective: np=2 must finish (no hang) and
    agree with np=1 to the platform KSP noise floor."""
    if shutil.which("mpirun") is None:
        pytest.skip("mpirun not on PATH")
    for f in ("minimal.yml", "mesh.npz"):       # minimal.yml reads mesh.npz
        shutil.copy(FIXTURES_DIR / f, tmp_path / f)

    outs = {}
    for n in (1, 2):
        out = tmp_path / f"run{n}.jsonl"
        res = subprocess.run(
            ["mpirun", "-n", str(n), sys.executable, "-m", "gospl", "-i",
             "minimal.yml", "--summary", str(out)],
            cwd=tmp_path, timeout=MPI_TIMEOUT, capture_output=True, text=True,
            env=mpi_child_env())
        assert res.returncode == 0, res.stdout[-2000:] + res.stderr[-2000:]
        outs[n] = _read(out)[1:-1]

    assert len(outs[1]) == len(outs[2])
    for a, b in zip(outs[1], outs[2]):
        assert b["nonfinite"] == {"elev": 0, "cumED": 0, "FA": 0}
        assert b["elev"]["max"] == pytest.approx(a["elev"]["max"], rel=1e-4)
        assert b["elev"]["mean"] == pytest.approx(a["elev"]["mean"], rel=1e-6)
