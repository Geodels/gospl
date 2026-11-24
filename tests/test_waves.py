"""
Wave energy and shadows (``waves:``, docs/DESIGN_CARBONATES.md 2.10, phase
C7). The exact solutions are checked in ``benchmarks/test_waves.py``.
"""
from __future__ import annotations

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


def _flat_case(path, waves, carb=""):
    dx, lx, ly = 500.0, 30.0e3, 15.0e3
    nx, ny = int(round(lx / dx)) + 1, int(round(ly / dx)) + 1
    X, Y = np.meshgrid(np.arange(nx) * dx, np.arange(ny) * dx)
    x, y = X.ravel(), Y.ravel()
    cells = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            a = j * nx + i
            cells += [[a, a + 1, a + nx + 1], [a, a + nx + 1, a + nx]]
    z = -30.0 + 2.0e-3 * x                       # shoaling towards the east
    z[(np.abs(x - 8.0e3) <= 1.5e3) & (np.abs(y - 7.5e3) <= 2.0e3)] = 20.0   # island
    np.savez(path / "flatw.npz", v=np.column_stack([x, y, np.zeros_like(x)]),
             c=np.array(cells, dtype=np.int64), z=z)
    (path / "w.yml").write_text(
        "name: waves\n"
        "domain:\n    npdata: ['flatw','v','c','z']\n    flowdir: 1\n    bc: 'wwww'\n"
        "    seadepo: True\n"
        "time:\n    start: 0.\n    end: 20.\n    tout: 10.\n    dt: 10.\n    strat: 10.\n"
        "spl:\n    K: 1.0e-20\n    d: 0.\n    m: 0.5\n"
        "diffusion:\n    hillslopeKa: 0.\n    hillslopeKm: 0.\n"
        "sea:\n    position: 0.\n"
        "climate:\n  - start: 0.\n    uniform: 1.\n"
        "strata:\n    dual: True\n" + waves + carb +
        "output:\n    dir: 'w_out'\n")
    return "w.yml"


WAVES = ("waves:\n    decay: 5.e4\n    shoal_depth: 20.\n    shoal_length: 2.e3\n"
         "    climate:\n      - {from: 250., weight: 3.}\n      - {from: 300., weight: 1.}\n")


def test_waves_parse_and_response(tmp_path, monkeypatch):
    """The climate weights are normalised; land gets no energy, the lee of
    the island less than the open water; the `wave` response scales a
    factory; inconsistent inputs are refused."""
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    carb = ("carbonate:\n    factories:\n      - {name: reef, growth: 500., "
            "response: [{control: wave, x: [0., 0.5], f: [0., 1.]}]}\n")
    m = Model(_flat_case(tmp_path, WAVES, carb), verbose=False, showlog=False)
    try:
        assert m.waveOn and np.allclose([w for _, w in m.waveClimate], [0.75, 0.25])
        m.sealevel = 0.0
        E = m.waveEnergy()
        h = m.hLocal.getArray()
        x, y = m.lcoords[:, 0], m.lcoords[:, 1]
        own = m.inIDs == 1
        assert np.all(E[own & (h >= 0.0)] == 0.0)
        assert np.all((E >= 0.0) & (E <= 1.0))
        lee = own & (np.abs(y - 7.5e3) < 1.0e3) & (x > 10.0e3) & (x < 12.0e3) & (h < 0.0)
        opn = own & (np.abs(y - 7.5e3) > 5.0e3) & (x > 10.0e3) & (x < 12.0e3) & (h < 0.0)
        if lee.any() and opn.any():
            assert E[lee].mean() < 0.5 * E[opn].mean()
        st = m._carbStatic(h)
        assert np.allclose(st[0], np.interp(m.waveE, [0.0, 0.5], [0.0, 1.0]))
    finally:
        m.destroy()
    for waves, msg in (("waves:\n    climate: []\n", "non-empty"),
                       ("waves:\n    climate:\n      - {weight: 1.}\n", "from"),
                       (WAVES.replace("decay: 5.e4", "decay: -1."), "positive")):
        with pytest.raises(ValueError, match=msg):
            Model(_flat_case(tmp_path, waves), verbose=False, showlog=False)
    with pytest.raises(ValueError, match="waves"):
        Model(_flat_case(tmp_path, "", carb), verbose=False, showlog=False)


_CHILD = '''\
import sys
import numpy as np
from mpi4py import MPI
from gospl.model import Model

m = Model(sys.argv[1], verbose=False, showlog=False)
try:
    m.sealevel = 0.0
    E = m.waveEnergy()
    own = m.inIDs == 1
    rows = np.column_stack([m.locIDs[own].astype(float), E[own]])
    allr = MPI.COMM_WORLD.gather(rows, root=0)
    if MPI.COMM_WORLD.rank == 0:
        a = np.vstack(allr)
        np.save(sys.argv[2], a[np.argsort(a[:, 0])])
finally:
    m.destroy()
'''


@pytest.mark.mpi
@pytest.mark.slow
def test_waves_partition_invariant(tmp_path):
    """The wave energy (two directions, an island, shoaling) is the same at
    np=1 and np=3."""
    if shutil.which("mpirun") is None:
        pytest.skip("mpirun not on PATH")
    yml = _flat_case(tmp_path, WAVES)
    child = tmp_path / "_w.py"
    child.write_text(_CHILD)
    out = {}
    for n in (1, 3):
        f = tmp_path / f"w{n}.npy"
        r = subprocess.run(["mpirun", "-n", str(n), sys.executable, str(child), yml, str(f)],
                           cwd=tmp_path, capture_output=True, text=True,
                           env=mpi_child_env(), timeout=MPI_TIMEOUT)
        assert r.returncode == 0, r.stderr[-2000:]
        out[n] = np.load(f)
    a, b = out[1], out[3]
    assert np.array_equal(a[:, 0], b[:, 0])
    d = np.max(np.abs(a[:, 1] - b[:, 1]))
    assert d < 1.0e-8, f"np=1 vs np=3 wave energy differs by {d:.3e}"


def _global_case(tmp_path, waves):
    src = FIXTURES_DIR / "minimal_dual.yml"
    shutil.copy(FIXTURES_DIR / "mesh.npz", tmp_path / "mesh.npz")
    import re
    text = re.sub(r"dir: '[^']*'", "dir: 'g_out'", src.read_text())
    text = text.replace("output:", waves + "output:", 1)
    (tmp_path / "g.yml").write_text(text)
    return "g.yml"


GWAVES = ("waves:\n    decay: 1.e7\n    fetch: 5.e5\n    shoal_depth: 20.\n"
          "    shoal_length: 2.e3\n    climate:\n      - {from: zonal, weight: 1.}\n")


def test_waves_global(tmp_path, monkeypatch):
    """On a global mesh (no edges) the waves grow over the fetch: the energy
    is in [0, 1], zero on land, close to g / (g + 1/decay) in the open ocean
    far from land, and the zonal belts differ from a single bearing; a global
    mesh without fetch, and zonal belts on a flat mesh, are refused."""
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    m = Model(_global_case(tmp_path, GWAVES), verbose=False, showlog=False)
    try:
        E = m.waveEnergy()
        h = m.hLocal.getArray()
        own = m.inIDs == 1
        assert np.all((E >= 0.0) & (E <= 1.0)) and np.all(E[h >= m.sealevel] == 0.0)
        einf = (1.0 / 5.0e5) / (1.0 / 5.0e5 + 1.0e-7)
        deep = own & (h < m.sealevel - 4000.0)
        assert deep.any() and np.percentile(E[deep], 90) > 0.8 * einf
        Ez = E.copy()
        m.waveClimate = [(90.0, 1.0)]
        m._waveCoef = None
        m._waveBear = None          # cached bearings (refraction)
        m._waveT = None             # per-step cache of the field (coastal erosion)
        Eb = m.waveEnergy()
        assert np.max(np.abs(Ez - Eb)) > 1.0e-3
    finally:
        m.destroy()
    with pytest.raises(ValueError, match="fetch"):
        Model(_global_case(tmp_path, GWAVES.replace("    fetch: 5.e5\n", "")),
              verbose=False, showlog=False)
    with pytest.raises(ValueError, match="zonal"):
        Model(_flat_case(tmp_path, WAVES.replace("{from: 250., weight: 3.}",
                                                 "{from: zonal, weight: 3.}")),
              verbose=False, showlog=False)


_GCHILD = _CHILD


@pytest.mark.mpi
@pytest.mark.slow
def test_waves_global_partition_invariant(tmp_path):
    """Zonal waves with a fetch on the global mesh: the same at np=1 and np=3."""
    if shutil.which("mpirun") is None:
        pytest.skip("mpirun not on PATH")
    yml = _global_case(tmp_path, GWAVES)
    child = tmp_path / "_g.py"
    child.write_text(_GCHILD)
    out = {}
    for n in (1, 3):
        f = tmp_path / f"g{n}.npy"
        r = subprocess.run(["mpirun", "-n", str(n), sys.executable, str(child), yml, str(f)],
                           cwd=tmp_path, capture_output=True, text=True,
                           env=mpi_child_env(), timeout=MPI_TIMEOUT)
        assert r.returncode == 0, r.stderr[-2000:]
        out[n] = np.load(f)
    a, b = out[1], out[3]
    assert np.array_equal(a[:, 0], b[:, 0])
    assert np.max(np.abs(a[:, 1] - b[:, 1])) < 1.0e-8
