"""
Karst dissolution of exposed carbonate (opt-in ``karst:`` block,
docs/DESIGN_KARST.md): runoff control, focusing by the discharge,
lithology-selective removal from the stratigraphy with a residue cover, and
the dissolved carbonate in the river load of the geochemistry.
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

# bedrock: `carb` (a coarse-like carbonate lithology) for `share` of it,
# the rest coarse clastic; no fine
LIMESTONE = ("    bedrock_coarse_frac: 1.0\n    lithologies:\n"
             "      - {{name: carb, size: coarse, phi0: 0.3, bedrock_frac: {share}}}\n")
KARST = "karst:\n    rate: 1.e-4\n    lithologies: [carb]\n{extra}\n"


def _case(tmp_path, name, karst=KARST, extra="", share=1.0, quiet=True):
    src = FIXTURES_DIR / "minimal_dual.yml"
    if not src.exists() or not (FIXTURES_DIR / "mesh.npz").exists():
        pytest.skip("minimal_dual.yml / mesh.npz not present")
    if not (tmp_path / "mesh.npz").exists():
        shutil.copy(FIXTURES_DIR / "mesh.npz", tmp_path / "mesh.npz")
    text = re.sub(r"dir: '[^']*'", f"dir: '{name}_out'", src.read_text())
    text = text.replace("    bedrock_coarse_frac: 0.6\n", LIMESTONE.format(share=share))
    if quiet:
        # no river or hillslope erosion: the karst alone moves the surface
        text = text.replace("K: 4.0e-6", "K: 1.0e-30")
        text = text.replace("hillslopeKa: 0.01", "hillslopeKa: 0.0")
        text = text.replace("hillslopeKm: 0.2", "hillslopeKm: 0.0")
    if karst is not None:
        text = text.replace("output:", karst.replace("{extra}", extra) + "output:", 1)
    (tmp_path / f"{name}.yml").write_text(text)
    return f"{name}.yml"


def test_karst_parse(tmp_path, monkeypatch):
    """Defaults, lithology lookup and the refused inputs."""
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    m = Model(_case(tmp_path, "p"), verbose=False, showlog=False)
    try:
        assert m.karstOn and m.karst["lit"] == [m.litName.index("carb")]
        assert m.karst["focus"] == 0.0 and m.karst["cover"] == 0.5
        assert m.karst["area"] > 0.0 and m.karstDiss is not None
    finally:
        m.destroy()
    bad = [
        (dict(karst="karst:\n    lithologies: [chalk]\n\n"), "unknown lithologies"),
        (dict(karst="karst:\n    rate: -1.\n    lithologies: [carb]\n\n"), "positive"),
        (dict(karst="karst:\n    lithologies: [carb]\n    temperature: {x: [0.], "
                    "f: [1.]}\n\n"), "sst"),
        (dict(karst="karst:\n    lithologies: [carb]\n    species: carbonate\n\n"),
         "geochemistry"),
    ]
    for kw, msg in bad:
        with pytest.raises(ValueError, match=msg):
            Model(_case(tmp_path, "bad", **kw), verbose=False, showlog=False)


def test_karst_off_is_identical(tmp_path, monkeypatch):
    """Protects: without a `karst:` block nothing changes (no state, no
    output field)."""
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    m = Model(_case(tmp_path, "off", karst=None), verbose=False, showlog=False)
    try:
        assert not m.karstOn and m.karstDiss is None
        m.runProcesses()
    finally:
        m.destroy()
    import h5py
    with h5py.File(tmp_path / "off_out" / "h5" / "gospl.1.p0.h5", "r") as f:
        assert "karstDiss" not in f


def test_karst_uniform_lowering(tmp_path, monkeypatch):
    """
    Protects the runoff law: with no focusing, uniform rain R and a pure
    carbonate bedrock, every land cell is lowered by exactly rate * R * t
    (and the sea floor not at all); the dissolved bulk and solid volumes
    match the lowering and the porosity of the lithology.
    """
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    m = Model(_case(tmp_path, "u"), verbose=False, showlog=False)
    try:
        h0 = m.hLocal.getArray().copy()
        m.runProcesses()
        dh = m.hLocal.getArray() - h0
        t = m.tEnd - m.tStart
        expect = 1.0e-4 * 1.0 * t
        # lakes do not dissolve; nothing else moves, so they keep their cells
        land = (h0 > m.sealevel + 10.0 * expect) & ~m._underLake(m.hLocal.getArray())
        sea = h0 < m.sealevel - 10.0 * expect
        assert land.sum() > 10 and sea.sum() > 10
        np.testing.assert_allclose(dh[land], -expect, rtol=1e-9)
        np.testing.assert_allclose(m.karstDiss[land], expect, rtol=1e-9)
        assert np.all(m.karstDiss[sea] == 0.0)
        own = m.inIDs == 1
        bulk = float((m.karstDiss * m.larea)[own].sum())
        assert abs(m._karstBulk - bulk) < 1e-9 * bulk
        assert abs(m._karstSolid.sum() - (1.0 - 0.3) * bulk) < 1e-9 * bulk
    finally:
        m.destroy()


def test_karst_selective_with_residue_cover(tmp_path, monkeypatch):
    """
    Protects the lithology selectivity: only the soluble piles are removed,
    the insoluble rest of a layer stays, and a residue thicker than `cover`
    stops the dissolution of the layers below.
    """
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    m = Model(_case(tmp_path, "s"), verbose=False, showlog=False)
    try:
        c = m.litName.index("carb")
        n = m.lpoints
        m.stratStep = 2
        # layer 2: 0.3 m pure carbonate; layer 1: 1 m of 60 % carbonate
        m.stratH[:, 2] = 0.3
        m.stratHL[c][:, 2] = 0.3
        m.stratH[:, 1] = 1.0
        m.stratHL[c][:, 1] = 0.6
        for l in range(m.litNb):
            m.phiL[l][:, 1:3] = 0.3
        H0 = m.stratH[:, :3].copy()
        rem, solid = m.dissolveStrat(np.full(n, 2.0), [c], cover=0.2)
        # layer 2 fully, layer 1's carbonate, then its 0.4 m residue stops it
        np.testing.assert_allclose(rem, 0.9)
        np.testing.assert_allclose(m.stratH[:, 2], 0.0)
        np.testing.assert_allclose(m.stratH[:, 1], 0.4)
        np.testing.assert_allclose(m.stratHL[c][:, 1], 0.0)
        np.testing.assert_allclose(m.stratH[:, 0], H0[:, 0])
        np.testing.assert_allclose(solid[:, 0], 0.9 * 0.7)
        # a thick cover lets it reach the bedrock
        m.stratH[:, 2] = 0.3
        m.stratHL[c][:, 2] = 0.3
        m.stratH[:, 1] = 1.0
        m.stratHL[c][:, 1] = 0.6
        rem, _ = m.dissolveStrat(np.full(n, 2.0), [c], cover=5.0)
        np.testing.assert_allclose(rem, 2.0)
        np.testing.assert_allclose(m.stratH[:, 0], H0[:, 0] - 1.1)
    finally:
        m.destroy()


def test_karst_partial_carbonate_bedrock(tmp_path, monkeypatch):
    """Protects: a bedrock of half carbonate loses its carbonate only; the
    coarse clastic pile of the layer is unchanged."""
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    m = Model(_case(tmp_path, "h", share=0.5), verbose=False, showlog=False)
    try:
        c = m.litName.index("carb")
        clas0 = m.stratH[:, 0] - sum(h[:, 0] for h in m.stratHL)
        m.runProcesses()
        k = m.stratStep
        clas = m.stratH[:, :k + 1].sum(axis=1) - sum(h[:, :k + 1].sum(axis=1)
                                                      for h in m.stratHL)
        np.testing.assert_allclose(clas, clas0, rtol=1e-12)
        assert m.karstDiss.max() > 0.0
        drop = m.stratHL[c][:, 0] - 0.5e6
        land = m.karstDiss > 0.0
        np.testing.assert_allclose(-drop[land], m.karstDiss[land], rtol=1e-6)
    finally:
        m.destroy()


def test_karst_focusing(tmp_path, monkeypatch):
    """Protects the focusing: the rate is rate * R * max(Q / (R A0), 1)^p,
    larger along the drainage than on the divides."""
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    m = Model(_case(tmp_path, "f", extra="    focus: 0.5\n    focus_area: 1.e6\n",
                    quiet=False), verbose=False, showlog=False)
    try:
        from gospl.flow.flowplex import FAMesh
        FAMesh.flowAccumulation(m)
        D = m._karstRate()
        Q = m.FAL.getArray()
        h = m.hLocal.getArray()
        land = (h > m.sealevel) & ~m._underLake(h)
        exp = 1.0e-4 * np.maximum(Q / 1.0e6, 1.0) ** 0.5
        np.testing.assert_allclose(D[land], exp[land], rtol=1e-12)
        assert D[land].max() > 3.0 * D[land].min()
    finally:
        m.destroy()


def test_karst_feeds_the_river_load(tmp_path, monkeypatch):
    """
    Protects the geochemistry link: the dissolved carbonate (solid volume
    times solute_per_m3) is added to the river routing source of the named
    species. The routing is linear in its source, so on the same state the
    ocean delivery with the karst minus without it is the routed karst load:
    positive and at most the dissolved amount (closed basins trap some).
    """
    from gospl.model import Model
    from gospl.flow.gwplex import GWMesh

    src = FIXTURES_DIR / "minimal_gw_river2.yml"
    if not src.exists():
        pytest.skip("minimal_gw_river2.yml not present")
    for f in ("mesh.npz", "soiltemp.npz"):
        shutil.copy(FIXTURES_DIR / f, tmp_path / f)
    monkeypatch.chdir(tmp_path)
    text = re.sub(r"dir: '[^']*'", "dir: 'g_out'", src.read_text())
    text = text.replace(
        "sea:\n", "strata:\n    dual: True\n    bedrock_coarse_frac: 1.0\n"
        "    lithologies:\n      - {name: carb, size: coarse, phi0: 0.3, "
        "bedrock_frac: 1.0}\n\nkarst:\n    rate: 1.e-4\n    lithologies: [carb]\n"
        "    species: carbonate\n    solute_per_m3: 2.\n\nsea:\n", 1)
    (tmp_path / "g.yml").write_text(text)
    m = Model("g.yml", verbose=False, showlog=False)
    rec = []
    orig = GWMesh._routeRiverSolute

    def spy(self):
        ks = self._karstSolute
        self._karstSolute = None
        orig(self)
        base = self.riverSoluteToOceanSp.copy()
        self._karstSolute = ks
        orig(self)
        own = self.inIDs == 1
        rec.append((base[self._karstSp], self.riverSoluteToOceanSp[self._karstSp],
                    float(ks[own].sum()), float(self.karstBudget[1:].sum())))
    monkeypatch.setattr(GWMesh, "_routeRiverSolute", spy)
    try:
        assert m._karstSp == list(m.gwGeoName).index("carbonate")
        m.runProcesses()
    finally:
        m.destroy()
    assert rec
    for base, kar, sol, solid in rec:
        assert abs(sol - 2.0 * solid / m.dt) < 1e-9 * sol
        assert sol > 0.0
        assert 0.0 < kar - base <= sol * (1 + 1e-6), (base, kar, sol)


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
    rows = np.column_stack([m.locIDs[own].astype(float), m.hLocal.getArray()[own],
                            m.karstDiss[own]] + [h[own, :top].sum(axis=1) for h in m.stratHL])
    allr = MPI.COMM_WORLD.gather(rows, root=0)
    if MPI.COMM_WORLD.rank == 0:
        a = np.vstack(allr)
        np.save(sys.argv[2], a[np.argsort(a[:, 0])])
finally:
    m.destroy()
'''


@pytest.mark.mpi
@pytest.mark.slow
def test_karst_partition_invariant(tmp_path):
    """Protects: the karst gives the same surface and dissolution at np=1
    and np=3 (focusing on, with river and hillslope erosion)."""
    yml = _case(tmp_path, "part", extra="    focus: 0.5\n", share=0.7, quiet=False)
    child = tmp_path / "child.py"
    child.write_text(_CHILD)
    out = {}
    for n in (1, 3):
        f = tmp_path / f"k{n}.npy"
        cmd = ["mpirun", "-n", str(n), sys.executable, str(child), yml, str(f)]
        r = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True,
                           env=mpi_child_env(), timeout=MPI_TIMEOUT)
        assert r.returncode == 0, r.stderr[-2000:]
        out[n] = np.load(f)
    a, b = out[1], out[3]
    assert np.array_equal(a[:, 0], b[:, 0])
    assert a[:, 2].max() > 0.0
    d = np.max(np.abs(a[:, 1:] - b[:, 1:]))
    assert d < 1.0e-5, f"np=1 vs np=3 differ by {d:.3e} m"
