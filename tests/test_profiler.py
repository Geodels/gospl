"""
Unit tests for the wall-clock phase profiler (gospl/tools/profiler.py).

These guard the no-op-when-disabled contract and the accumulation / cross-rank
reduction math. They run single-rank (the cross-rank Allreduce in a size-1
communicator is the identity, so the reported min == mean == max), which is
enough to lock down the reporting arithmetic; true multi-rank imbalance is
exercised by the scaling harness on HPC, not in the fast tier.
"""

import json

import pytest

from gospl.tools.profiler import Profiler


def test_disabled_is_noop(tmp_path):
    """A disabled profiler records nothing and never communicates."""
    p = Profiler(enabled=False)
    with p.phase("flow"):
        pass
    p.start("erosion")
    p.stop("erosion")
    assert p.times == {}
    assert p.counts == {}
    assert p.reduce(total_wall=1.0) is None
    assert p.report(outputDir=str(tmp_path), total_wall=1.0) is None
    # No file written when disabled.
    assert not (tmp_path / "profile.json").exists()


def test_accumulation_counts_and_order():
    """Repeated phases accumulate time and call counts; order is first-seen."""
    p = Profiler(enabled=True)
    for _ in range(3):
        with p.phase("flow"):
            pass
    with p.phase("erosion"):
        pass
    assert p.counts["flow"] == 3
    assert p.counts["erosion"] == 1
    assert p._order == ["flow", "erosion"]
    # Manual accumulation path (start/stop) charges the same bucket.
    p.start("flow")
    p.stop("flow")
    assert p.counts["flow"] == 4


def test_reduce_math_single_rank():
    """On one rank, min == mean == max == sum and imbalance == 1.0."""
    p = Profiler(enabled=True)
    # Inject deterministic timings instead of sleeping (fast + stable).
    p._accumulate("flow", 0.4)
    p._accumulate("flow", 0.2)
    p._accumulate("erosion", 0.1)
    stats = p.reduce(total_wall=1.0)
    assert stats["__meta__"]["nranks"] == 1
    assert stats["__meta__"]["total_wall"] == pytest.approx(1.0)
    f = stats["flow"]
    assert f["calls"] == 2
    assert f["sum"] == pytest.approx(0.6)
    assert f["mean"] == pytest.approx(0.6)
    assert f["min"] == pytest.approx(0.6)
    assert f["max"] == pytest.approx(0.6)
    assert f["imbalance"] == pytest.approx(1.0)
    # A never-recorded phase has zero mean and zero (not NaN) imbalance.
    p2 = Profiler(enabled=True)
    p2._accumulate("x", 0.0)
    assert p2.reduce()["x"]["imbalance"] == 0.0


def test_report_writes_json(tmp_path):
    """report() writes a machine-readable profile.json consumed by the harness."""
    p = Profiler(enabled=True)
    p._accumulate("flow", 0.5)
    p._accumulate("sea", 0.3)
    p.report(outputDir=str(tmp_path), total_wall=1.0, to_stdout=False)
    path = tmp_path / "profile.json"
    assert path.exists()
    data = json.loads(path.read_text())
    assert set(data.keys()) == {"__meta__", "flow", "sea"}
    assert data["__meta__"]["nranks"] == 1
    assert data["flow"]["mean"] == pytest.approx(0.5)
