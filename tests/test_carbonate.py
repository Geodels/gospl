"""
Benthic carbonate production (opt-in ``carbonate:`` block,
docs/DESIGN_CARBONATES.md, phase C1). The production law itself is checked
against its exact solutions in ``benchmarks/test_carbonate.py``.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys

import numpy as np
import pytest

from _helpers import FIXTURES_DIR, MPI_TIMEOUT, mpi_child_env

pytest.importorskip(
    "gospl._fortran", reason="goSPL Fortran extension not built"
)

pytestmark = [pytest.mark.sediment]

CARB = "carbonate:\n    insolation: 400.\n\n"


def _case(tmp_path, name, block=CARB, dual=True):
    src = FIXTURES_DIR / "minimal_dual.yml"
    if not src.exists() or not (FIXTURES_DIR / "mesh.npz").exists():
        pytest.skip("minimal_dual.yml / mesh.npz not present")
    if not (tmp_path / "mesh.npz").exists():
        shutil.copy(FIXTURES_DIR / "mesh.npz", tmp_path / "mesh.npz")
    text = re.sub(r"dir: '[^']*'", f"dir: '{name}_out'", src.read_text())
    if not dual:
        text = text.replace("    dual: True\n", "    dual: False\n")
    text = text.replace("output:", block + "output:", 1)
    (tmp_path / f"{name}.yml").write_text(text)
    return f"{name}.yml"


def test_carbonate_parse(tmp_path, monkeypatch):
    """Default factories and the carbonate lithologies are added to the
    lithology table; inconsistent inputs are refused."""
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    m = Model(_case(tmp_path, "p"), verbose=False, showlog=False)
    try:
        assert m.carbOn
        assert [f["name"] for f in m.carbFact] == ["euphotic", "oligophotic", "aphotic"]
        assert m.litName == ["fine", "carb_sand", "carb_mud", "reef"]
        assert m.litFine == [True, False, True, False]
        assert abs(m.carbFact[0]["g"] - 500.0e-6) < 1e-18      # m/Myr -> m/yr
        assert m.carbProd is not None and m._litMulti
    finally:
        m.destroy()
    bad = [
        (dict(dual=False), "dual lithology"),
        (dict(block="carbonate:\n    factories:\n      - {name: a, growth: 1., "
                    "response: [{control: wind, x: [0, 1], f: [1, 0]}]}\n\n"), "control"),
        (dict(block="carbonate:\n    factories:\n      - {name: a, framework: 2.}\n\n"),
         "between 0 and 1"),
        (dict(block="carbonate:\n    lithologies: {limestone: {phi0: 0.3}}\n\n"),
         "no 'limestone'"),
    ]
    for kw, msg in bad:
        yml = _case(tmp_path, "bad", **kw)
        with pytest.raises(ValueError, match=msg):
            Model(yml, verbose=False, showlog=False)


@pytest.mark.slow
def test_carbonate_budget(tmp_path, monkeypatch):
    """
    Protects: what is produced is accounted for. Production = framework +
    loose; the framework is recorded as `reef` exactly; the loose sediment is
    deposited as carb_sand + carb_mud, minus what leaves through the open
    edges; the clastic fractions stay conserved; the output field exists.
    """
    from gospl.model import Model
    import h5py

    monkeypatch.chdir(tmp_path)
    m = Model(_case(tmp_path, "b"), verbose=False, showlog=False)
    try:
        m.runProcesses()
        P, F, L = m._carbProduced, m._carbFramework, m._carbLoose
        dep = dict(zip(m.litName, m._litDeposited))
        ero = dict(zip(m.litName, m._litEroded))
        assert P > 0.0 and F > 0.0 and L > 0.0
        assert abs(P - F - L) < 1e-9 * P
        assert abs(dep["reef"] - F) < 1e-9 * F
        loose = dep["carb_sand"] + dep["carb_mud"]
        assert 0.98 * L < loose <= L * (1 + 1e-9), (loose, L)
        assert abs(ero["fine"] - dep["fine"]) < 2e-3 * ero["fine"]
        assert np.all(m.carbProd >= 0.0)
    finally:
        m.destroy()
    with h5py.File(tmp_path / "b_out" / "h5" / "gospl.1.p0.h5", "r") as f:
        assert "carbProd" in f and float(np.max(f["carbProd"][:])) > 0.0


def test_carbonate_responses(tmp_path, monkeypatch):
    """The response tables multiply the production: a sedrate table cuts it
    where the clastic deposition rate is high, a depth table inside the
    integration, a habitat map of zeros switches a factory off."""
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    block = ("carbonate:\n    factories:\n"
             "      - {name: a, growth: 500., extinction: 0.1, saturation: 60., "
             "response: [{control: sedrate, x: [0., 1.0e-3], f: [1., 0.]}, "
             "{control: depth, x: [0., 10., 20.], f: [1., 1., 0.]}]}\n"
             "      - {name: b, growth: 500., habitat: ['hab', 'h']}\n\n")
    yml = _case(tmp_path, "r", block=block)
    np.savez(tmp_path / "hab.npz", h=np.zeros(np.load(tmp_path / "mesh.npz")["v"].shape[0]))
    m = Model(yml, verbose=False, showlog=False)
    try:
        m._carbSedRate[:] = 0.0
        m._carbSedRate[::2] = 2.0e-3
        st = m._carbStatic(m.hLocal.getArray())
        assert np.all(st[0][::2] == 0.0) and np.all(st[0][1::2] == 1.0)
        assert np.all(st[1] == 0.0)
        w = np.array([5.0, 15.0, 25.0])
        r = m._carbRates(w, [np.ones(3), np.ones(3)])
        assert r[0][2] == 0.0 and 0.0 < r[0][1] < r[0][0]
        assert np.all(r[1] > 0.0)
    finally:
        m.destroy()


_CHILD = '''\
import sys
import numpy as np
from mpi4py import MPI
from gospl.model import Model

m = Model(sys.argv[1], verbose=False, showlog=False)
try:
    m.runProcesses()
    own = m.inIDs == 1
    top = m.stratStep + 1
    rows = np.column_stack([m.locIDs[own].astype(float), m.hLocal.getArray()[own]]
                           + [h[own, 1:top] for h in m.stratHL])
    allr = MPI.COMM_WORLD.gather(rows, root=0)
    if MPI.COMM_WORLD.rank == 0:
        a = np.vstack(allr)
        np.save(sys.argv[2], a[np.argsort(a[:, 0])])
finally:
    m.destroy()
'''


@pytest.mark.mpi
@pytest.mark.slow
def test_carbonate_partition_invariant(tmp_path):
    """Carbonate production, framework and the routed loose sediment give the
    same elevation and lithology piles at np=1 and np=3."""
    if shutil.which("mpirun") is None:
        pytest.skip("mpirun not on PATH")
    yml = _case(tmp_path, "part")
    child = tmp_path / "_c.py"
    child.write_text(_CHILD)
    out = {}
    for n in (1, 3):
        f = tmp_path / f"c{n}.npy"
        r = subprocess.run(["mpirun", "-n", str(n), sys.executable, str(child), yml, str(f)],
                           cwd=tmp_path, capture_output=True, text=True,
                           env=mpi_child_env(), timeout=MPI_TIMEOUT)
        assert r.returncode == 0, r.stderr[-2000:]
        out[n] = np.load(f)
    a, b = out[1], out[3]
    assert np.array_equal(a[:, 0], b[:, 0])
    d = np.max(np.abs(a[:, 1:] - b[:, 1:]))
    assert d < 1.0e-5, f"np=1 vs np=3 differ by {d:.3e} m"
