"""
gospl-inspect (gospl/analyse/runinspect.py): output-directory and run-summary
inspection used by agents and by `scripts/budget.py`.
"""

from __future__ import annotations

import json
import os
import shutil

import numpy as np
import pytest

from _helpers import FIXTURES_DIR

pytest.importorskip("h5py")
pytest.importorskip(
    "gospl.tools.inputparser",
    reason="goSPL runtime deps (petsc4py / ruamel.yaml / scipy) not installed",
)

pytestmark = pytest.mark.analyse


@pytest.fixture(scope="module")
def run_dir(tmp_path_factory):
    """minimal.yml run in a private directory: (output dir, summary file)."""
    from gospl.model import Model

    d = tmp_path_factory.mktemp("inspect")
    for f in ("minimal.yml", "mesh.npz"):
        shutil.copy(FIXTURES_DIR / f, d / f)
    cwd = os.getcwd()
    os.chdir(d)
    try:
        m = Model("minimal.yml", verbose=False, showlog=False,
                  summary=str(d / "run.jsonl"))
        try:
            m.runProcesses()
        finally:
            m.destroy()
    finally:
        os.chdir(cwd)
    return d / "minimal", d / "run.jsonl"


def test_inspect_output_dir(run_dir):
    from gospl.analyse import runinspect as ri

    out = ri.Output(str(run_dir[0]))
    assert out.input_ids                       # mesh found through the copied YAML
    assert out.npoints == len(np.load(FIXTURES_DIR / "mesh.npz")["v"])
    # barycentric areas tile the sphere (R = 6378137 m in mesh.npz)
    r = np.linalg.norm(out.coords, axis=1).mean()
    assert out.area.sum() == pytest.approx(4 * np.pi * r**2, rel=2e-2)

    stats = ri.field_stats(out, out.steps[-1])
    assert {"elev", "erodep", "FA"} <= set(stats)
    assert all(v["nonfinite"] == 0 for v in stats.values())

    rows = ri.budget(out)
    assert len(rows) == len(out.steps) - 1
    for r_ in rows:                            # closed sphere: net ~ 0
        assert abs(r_["net_rel"]) < 5e-3

    sp = ri.spikes(out)
    assert len(sp) == len(out.steps) - 1 and all(s["max_abs_dz"] >= 0 for s in sp)

    same = ri.compare(out, ri.Output(str(run_dir[0])), out.steps[-1])
    assert all(v["max_abs_diff"] == 0.0 for v in same.values())

    assert ri.main([str(run_dir[0]), "--json"]) == 0


def test_inspect_without_mesh_reassembles(run_dir, tmp_path):
    """No findable mesh: de-duplicated reassembly, same node count."""
    from gospl.analyse import runinspect as ri

    copy = tmp_path / "out"
    shutil.copytree(run_dir[0], copy)
    for y in copy.glob("*.yml"):
        y.unlink()                              # nothing to locate the mesh with
    out = ri.Output(str(copy))
    assert not out.input_ids
    assert out.npoints == ri.Output(str(run_dir[0])).npoints


def test_inspect_summary_flags_issues(run_dir, tmp_path):
    from gospl.analyse import runinspect as ri

    clean = ri.inspect_summary(str(run_dir[1]))
    assert clean["finished"] and clean["issues"] == []

    recs = [json.loads(l) for l in open(run_dir[1])]
    bad = recs[:3]
    bad[2] = dict(bad[2], nonfinite={"elev": 3, "cumED": 0, "FA": 0},
                  events=[{"kind": "flow_ksp_fallback_failed", "outcome": "ponded",
                           "nbad": 210}])
    p = tmp_path / "bad.jsonl"
    p.write_text("".join(json.dumps(r) + "\n" for r in bad))   # no end record
    rep = ri.inspect_summary(str(p))
    kinds = {i["issue"] for i in rep["issues"]}
    assert not rep["finished"]
    assert {"nonfinite", "flow_ksp_fallback_failed"} <= kinds
    assert ri.main(["--summary", str(p)]) == 1                 # non-finite -> rc 1
