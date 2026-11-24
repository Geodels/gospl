"""
Carbonates and the groundwater geochemistry (phase C6): supply-limited
burial from the marine reservoir of a species (``carbonate_supply``), and the
coastal plume of a species as a response control (``solute``).
"""
from __future__ import annotations

import re
import shutil

import numpy as np
import pytest

from _helpers import FIXTURES_DIR

pytest.importorskip(
    "gospl._fortran", reason="goSPL Fortran extension not built"
)

pytestmark = [pytest.mark.sediment]


def _case(tmp_path, name, block):
    src = FIXTURES_DIR / "minimal_gw_river2.yml"
    if not src.exists():
        pytest.skip("minimal_gw_river2.yml not present")
    for f in ("mesh.npz", "soiltemp.npz"):
        shutil.copy(FIXTURES_DIR / f, tmp_path / f)
    text = re.sub(r"dir: '[^']*'", f"dir: '{name}_out'", src.read_text())
    text = text.replace("sea:\n", "strata:\n    dual: True\n\n" + block + "sea:\n", 1)
    (tmp_path / f"{name}.yml").write_text(text)
    return f"{name}.yml"


CARB = ("carbonate:\n    factories:\n      - {name: a, growth: 500., extinction: 0.1, "
        "saturation: 60., framework: 0.5, sand: 0.5}\n\n")


@pytest.mark.slow
@pytest.mark.parametrize("yld", [1.0e6, 1.0e-6])
def test_carbonate_supply_limited_burial(tmp_path, monkeypatch, yld):
    """
    Protects: the global supply limit. Each step the reservoir of the species
    drops by (carbonate buried) / yield and never below zero; with a large
    yield the supply never limits (scale 1), with a tiny one the burial is
    what the reservoir allows (scale < 1, burial = reservoir * yield).
    """
    from gospl.model import Model
    from gospl.sed.carbplex import CarbMesh

    monkeypatch.chdir(tmp_path)
    block = CARB + "carbonate_supply:\n    species: carbonate\n    yield: %r\n\n" % yld
    m = Model(_case(tmp_path, "s", block), verbose=False, showlog=False)
    try:
        k = m._carbSupplyK
        steps = []
        orig = CarbMesh.carbSupplyLimit

        def spy(self):
            orig(self)
            steps.append(self.carbSupplyBudget.copy())
        monkeypatch.setattr(CarbMesh, "carbSupplyLimit", spy)
        try:
            res0 = []
            origP = CarbMesh.carbProduction

            def spyP(self):
                before = float(self.marineSolute[k])
                origP(self)
                res0.append((before, float(self.marineSolute[k]), float(self.carbBudget[0])))
            monkeypatch.setattr(CarbMesh, "carbProduction", spyP)
            m.runProcesses()
        finally:
            monkeypatch.undo()
        assert steps
        for before, after, buried in res0:
            assert after >= 0.0
            assert abs((before - after) - buried / yld) <= 1e-9 * max(before, 1.0) + 1e-12
        scales = [b[2] for b in steps]
        if yld > 1.0:
            assert all(s == 1.0 for s in scales)
        else:
            assert any(s < 1.0 for s in scales)
            for b, (before, after, buried) in zip(steps, res0):
                if b[2] < 1.0:
                    assert abs(buried - b[0]) <= 1e-6 * max(b[0], 1.0)
    finally:
        m.destroy()


def test_carbonate_solute_plume(tmp_path, monkeypatch):
    """The coastal plume of a species is largest at the river exits, zero
    on land, and drives the `solute` response; inconsistent couplings are
    refused."""
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    block = ("carbonate:\n    factories:\n      - {name: a, growth: 500., response: "
             "[{control: solute, species: carbonate, length: 2.e5, x: [0., 1.e-12], "
             "f: [1., 0.]}]}\n\n")
    m = Model(_case(tmp_path, "p", block), verbose=False, showlog=False)
    try:
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()                     # one step: river solute routed
        c = m._carbPlume("carbonate", 2.0e5)
        h = m.hLocal.getArray()
        assert np.all(c >= 0.0)
        exits = m.marineSoluteInput > 0.0
        if exits.any():
            assert c[exits].max() >= c[~exits & (h < m.sealevel)].max() - 1e-30
        st = m._carbStatic(h)
        assert np.allclose(st[0], np.interp(c, [0.0, 1.0e-12], [1.0, 0.0]))
    finally:
        m.destroy()
    with pytest.raises(ValueError, match="species"):
        Model(_case(tmp_path, "b1", CARB + "carbonate_supply:\n    yield: 1.\n\n"),
              verbose=False, showlog=False)
    with pytest.raises(ValueError, match="not found"):
        Model(_case(tmp_path, "b2", CARB + "carbonate_supply:\n    species: gypsum\n"
                    "    yield: 1.\n\n"), verbose=False, showlog=False)


@pytest.mark.slow
def test_marine_reservoir_restart(tmp_path, monkeypatch):
    """The marine reservoir (drawn down by the carbonate burial) is written
    in the outputs and restored on restart."""
    from gospl.model import Model
    import h5py

    monkeypatch.chdir(tmp_path)
    block = CARB + "carbonate_supply:\n    species: carbonate\n    yield: 1.0e-6\n\n"
    yml = _case(tmp_path, "r", block)
    m = Model(yml, verbose=False, showlog=False)
    try:
        m.runProcesses()
    finally:
        m.destroy()
    with h5py.File(tmp_path / "r_out" / "h5" / "gospl.1.p0.h5", "r") as f:
        saved = np.array(f["marineSolute"])
    text = (tmp_path / yml).read_text().replace("    start: 0.\n", "    start: 0.\n    rstep: 1\n", 1)
    (tmp_path / "rr.yml").write_text(text)
    mr = Model("rr.yml", verbose=False, showlog=False)
    try:
        assert np.allclose(mr.marineSolute, saved)
    finally:
        mr.destroy()


_CHILD = """\
import sys
import numpy as np
from mpi4py import MPI
from gospl.model import Model

m = Model(sys.argv[1], verbose=False, showlog=False)
try:
    m.runProcesses()
    own = m.inIDs == 1
    rows = np.column_stack([m.locIDs[own].astype(float), m.hLocal.getArray()[own]])
    allr = MPI.COMM_WORLD.gather(rows, root=0)
    if MPI.COMM_WORLD.rank == 0:
        a = np.vstack(allr)
        np.save(sys.argv[2], a[np.argsort(a[:, 0])])
        np.save(sys.argv[2] + ".res.npy", np.asarray(m.marineSolute))
finally:
    m.destroy()
"""


@pytest.mark.mpi
@pytest.mark.slow
def test_carbonate_supply_partition_invariant(tmp_path):
    """Supply-limited burial (a decision on global sums) gives the same
    elevation and reservoir at np=1 and np=3."""
    import subprocess
    import sys

    from _helpers import MPI_TIMEOUT, mpi_child_env

    if shutil.which("mpirun") is None:
        pytest.skip("mpirun not on PATH")
    block = CARB + "carbonate_supply:\n    species: carbonate\n    yield: 1.0e-6\n\n"
    yml = _case(tmp_path, "pp", block)
    # without the soil model: with it, carbonate runs on this fixture differ
    # by ~7e-4 m between np=1 and np=3 (open item, see gospl/sed/AGENTS.md),
    # which is not what this test is about
    text = re.sub(r"\nsoil:\n(    .*\n)*", "\n", (tmp_path / yml).read_text())
    (tmp_path / yml).write_text(text)
    child = tmp_path / "_c.py"
    child.write_text(_CHILD)
    out, res = {}, {}
    for n in (1, 3):
        f = tmp_path / f"c{n}.npy"
        r = subprocess.run(["mpirun", "-n", str(n), sys.executable, str(child), yml, str(f)],
                           cwd=tmp_path, capture_output=True, text=True,
                           env=mpi_child_env(), timeout=MPI_TIMEOUT)
        assert r.returncode == 0, r.stderr[-2000:]
        out[n] = np.load(f)
        res[n] = np.load(str(f) + ".res.npy")
    assert np.array_equal(out[1][:, 0], out[3][:, 0])
    assert np.max(np.abs(out[1][:, 1] - out[3][:, 1])) < 1.0e-5
    # the reservoir of the supply species (index 0, `carbonate`); the other
    # species (silica, with an in-transit decay) is not touched by C6
    assert np.allclose(res[1][0], res[3][0], rtol=1e-6, atol=1e-6 * max(abs(res[1][0]), 1.0))
