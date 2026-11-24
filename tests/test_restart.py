"""
A run restarted from an output step reproduces the continuous run.

Restart had three faults (2026-11-24): the stratigraphic layer counter and the
next-layer time were one interval ahead (a restarted run skipped one layer and
its compaction), the first sea level was taken one step ahead, and the
surface fields read back were stored in single precision (the 1-ulp tilts of
the depression fill were lost and flow paths switched on flat ground).
"""
from __future__ import annotations

import os
import re
import shutil

import numpy as np
import pytest

from _helpers import FIXTURES_DIR

pytest.importorskip(
    "gospl._fortran", reason="goSPL Fortran extension not built"
)

pytestmark = [pytest.mark.sediment]


def _case(tmp_path, name, extra="", restart=False):
    src = FIXTURES_DIR / "minimal_dual.yml"
    if not src.exists() or not (FIXTURES_DIR / "mesh.npz").exists():
        pytest.skip("minimal_dual.yml / mesh.npz not present")
    if not (tmp_path / "mesh.npz").exists():
        shutil.copy(FIXTURES_DIR / "mesh.npz", tmp_path / "mesh.npz")
        # a sea level rising 0.3 m/yr: an off-by-one-step level matters
        np.savetxt(tmp_path / "sea.csv", np.array([[0.0, -100.0], [400.0, 20.0]]),
                   delimiter=",", fmt="%.1f")
    text = re.sub(r"dir: '[^']*'", "dir: 'out'", src.read_text())
    text = text.replace("    end: 100\n", "    end: 200\n" + ("    rstep: 2\n" if restart else ""))
    text = text.replace("    strat: 50\n", "    strat: 25\n")         # 2 layers per output
    text = text.replace("    position: -100.\n", "    position: -100.\n    curve: 'sea.csv'\n")
    text = text.replace("output:", extra + "output:", 1)
    (tmp_path / f"{name}.yml").write_text(text)
    return f"{name}.yml"


def _run(tmp_path, yml):
    from gospl.model import Model

    m = Model(yml, verbose=False, showlog=False)
    try:
        m.runProcesses()
        top = m.stratStep + 1
        out = {"h": m.hLocal.getArray().copy(), "H": m.stratH[:, :top].copy(),
               "ed": m.cumEDLocal.getArray().copy(), "top": top}
    finally:
        m.destroy()
    return out


@pytest.mark.slow
def test_restart_reproduces_continuous_run(tmp_path, monkeypatch):
    """Restarted from the output at 100 yr, the run ends (200 yr) with the
    elevation, cumulative erosion-deposition and stratigraphy of the run that
    did not stop (to the solver tolerances), with dual lithology, two layers
    per output and a sea-level curve."""
    extra = ""
    monkeypatch.chdir(tmp_path)
    a = _run(tmp_path, _case(tmp_path, "full", extra))
    b = _run(tmp_path, _case(tmp_path, "rst", extra, restart=True))
    assert a["top"] == b["top"]
    scale = max(np.abs(a["h"]).max(), 1.0)
    np.testing.assert_allclose(b["h"], a["h"], rtol=0.0, atol=1.0e-6 * scale)
    np.testing.assert_allclose(b["ed"], a["ed"], rtol=0.0, atol=1.0e-6 * scale)
    np.testing.assert_allclose(b["H"], a["H"], rtol=0.0, atol=1.0e-6 * scale)
