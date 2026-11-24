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
# production favoured on margins and rims (open-water index, phase C3)
CARB_EXP = ("carbonate:\n    exposure_length: 3.e4\n    factories:\n"
            "      - {name: margin, growth: 500., extinction: 0.8, saturation: 60., "
            "framework: 0.5, sand: 0.7, response: [{control: exposure, "
            "x: [-50., 50.], f: [0.2, 1.]}]}\n"
            "      - {name: deep, growth: 100., extinction: 0.005, sand: 0.3}\n\n")


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
        (dict(block="carbonate:\n    factories:\n      - {name: a, growth: 1., "
                    "response: [{control: exposure, x: [0, 1], f: [1, 0]}]}\n\n"),
         "exposure_length"),
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


def test_carbonate_turbidity(tmp_path, monkeypatch):
    """The turbidity plume spreads the river fine load delivered to the sea:
    from a point source far from land it keeps the delivered volume (the
    finite-volume operator is conservative), peaks at the source and falls
    with distance; in a run it is zero on land and drives the response. A
    turbidity response without a length is refused."""
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    L = 4.0e5
    block = ("carbonate:\n    factories:\n      - {name: a, growth: 500., response: "
             "[{control: turbidity, length: %g, source: mouths, x: [0., 1.e-4], "
             "f: [1., 0.]}]}\n\n" % L)
    m = Model(_case(tmp_path, "t", block), verbose=False, showlog=False)
    try:
        assert m.carbFact[0]["turbidity"] == (L, "mouths")
        h = m.hLocal.getArray()
        xyz = m.lcoords
        land = h >= m.sealevel
        sea = (~land) & (m.lFill <= m.sealevel)
        # the sea cell farthest from land
        d = np.min(np.linalg.norm(xyz[sea][:, None, :] - xyz[land][None, ::4, :], axis=2),
                   axis=1)
        i = np.where(sea)[0][np.argmax(d)]
        assert d.max() > 6 * L
        v = np.zeros(m.lpoints)
        v[i] = 1.0e9                                     # m3/yr of fine load
        m.vSedFLocal.setArray(v)
        c = m._carbTurbidity(L, "mouths")
        assert np.all(c >= 0.0) and np.all(c[land] == 0.0)
        # conservative: only the tail reaching land (> 6 L away, ~exp(-6)
        # times the coastal cells' share) is absorbed by its identity rows
        np.testing.assert_allclose((c * m.larea).sum(), 1.0e9, rtol=1.0e-5)
        assert np.argmax(c) == i
        r = np.linalg.norm(xyz - xyz[i], axis=1)
        near, far = (r > 0) & (r < 1.5 * L), (r > 3 * L) & (r < 4.5 * L)
        assert c[near].mean() > 3.0 * c[far].mean() > 0.0
        st = m._carbStatic(h)
        assert np.allclose(st[0], np.interp(c, [0.0, 1.0e-4], [1.0, 0.0]))

        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()                                 # river load routed
        c = m._carbTurbidity(L, "mouths")
        assert np.all(c[m.hLocal.getArray() >= m.sealevel] == 0.0) and c.max() > 0.0
    finally:
        m.destroy()
    bad = block.replace("length: %g, " % L, "")
    with pytest.raises(ValueError, match="turbidity response needs"):
        Model(_case(tmp_path, "t2", bad), verbose=False, showlog=False)
    with pytest.raises(ValueError, match="turbidity source"):
        Model(_case(tmp_path, "t3", block.replace("source: mouths", "source: rivers")),
              verbose=False, showlog=False)


def test_carbonate_turbidity_transport(tmp_path, monkeypatch):
    """`source: transport` (default): the marine routing records the river
    mud carried through each cell; every deposited part of it passed through
    its node (throughflow >= mud deposit, both from the same routing), the
    river input counts at the nodes it enters, and the next step's plume uses
    it (zero on land). Without a transport turbidity response the
    routing records nothing."""
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    block = ("carbonate:\n    factories:\n      - {name: a, growth: 500., response: "
             "[{control: turbidity, length: 4.e5, x: [0., 1.e-3], f: [1., 0.]}]}\n\n")
    m = Model(_case(tmp_path, "tt", block), verbose=False, showlog=False)
    try:
        assert m.carbFact[0]["turbidity"] == (4.0e5, "transport") and m._turbTransport
        c0 = m._carbTurbidity(4.0e5)
        assert np.all(c0 == 0.0)                          # nothing routed yet
        dep = {}
        orig = type(m)._distOcean

        def spy(self, sedflux, provFlux=None):
            out = orig(self, sedflux, provFlux=provFlux)
            if getattr(self, "_turbCol", None) is not None:
                dep["mud"] = self._marDepProv[:, self._turbCol].copy()
                dep["in"] = np.maximum(provFlux[:, self._turbCol], 0.0)
            return out

        monkeypatch.setattr(type(m), "_distOcean", spy)
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()
        thr = m._turbThrough
        assert thr.sum() > 0.0 and "mud" in dep
        assert np.all(thr >= dep["mud"] - 1.0e-6 * max(dep["mud"].max(), 1.0))
        assert np.all(thr >= dep["in"])
        np.testing.assert_allclose(m._turbRate, thr / (m.larea * m.dt))
        c = m._carbTurbidity(4.0e5)
        assert c.max() > 0.0
        assert np.all(c[m.hLocal.getArray() >= m.sealevel] == 0.0)    # land rows
    finally:
        m.destroy()
    m = Model(_case(tmp_path, "tn", CARB), verbose=False, showlog=False)
    try:
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()
        assert getattr(m, "_turbThrough", None) is None
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
@pytest.mark.parametrize("pel", ["decay", "resuspend", "turbidity"])
def test_carbonate_partition_invariant(tmp_path, pel):
    """Carbonate production (with the open-water index, a parallel solve),
    framework, the routed loose sediment and the pelagic rain (clay decaying
    from the coast, another solve; settled or resuspended into the marine
    diffusion) give the same elevation and lithology piles at np=1 and
    np=3."""
    if shutil.which("mpirun") is None:
        pytest.skip("mpirun not on PATH")
    carb = CARB_EXP
    if pel == "turbidity":
        # the turbid plume of the river fine load (another parallel solve)
        carb = CARB_EXP.replace(
            "response: [{control: exposure,",
            "response: [{control: turbidity, length: 4.e5, x: [0., 1.e-3], "
            "f: [1., 0.2]}, {control: exposure,")
    yml = _case(tmp_path, "part",
                block=carb + (PEL_RESUSP if pel == "resuspend" else PEL_DECAY))
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


PEL = ("pelagic:\n  - start: 0.\n    carbonate: 50.\n    lysocline: 3500.\n"
       "    ccd: 4500.\n    clay: 20.\n\n")
# clay decaying from the coast (C2b), settled in place / resuspended
PEL_DECAY = PEL.replace("    clay: 20.\n", "    clay: 20.\n    clay_decay: 3.e5\n")
PEL_RESUSP = PEL_DECAY.replace("    clay_decay: 3.e5\n",
                               "    clay_decay: 3.e5\n    resuspend: True\n")


def test_pelagic_parse(tmp_path, monkeypatch):
    """pelagic: events are sorted and converted to m/yr, carb_mud is added
    to the lithology table, inconsistent events are refused."""
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    block = ("pelagic:\n  - start: 50.\n    clay: 2.\n"
             "  - start: 0.\n    carbonate: 5.\n    lysocline: 3000.\n    ccd: 4000.\n\n")
    m = Model(_case(tmp_path, "pp", block=block), verbose=False, showlog=False)
    try:
        assert m.pelOn and not m.carbOn
        assert [e["start"] for e in m.pelEvents] == [0.0, 50.0]
        assert abs(m.pelEvents[0]["carb"] - 5.0e-6) < 1e-18
        assert m.pelEvents[1]["carb"] == 0.0 and abs(m.pelEvents[1]["clay"] - 2e-6) < 1e-18
        assert m.litName == ["fine", "carb_mud"] and m.carbLit["carb_mud"] == 1
        m.tNow = 49.0
        assert m._pelEvent() == 0
        m.tNow = 50.0
        assert m._pelEvent() == 1
    finally:
        m.destroy()
    for block, msg in (
        ("pelagic:\n  - start: 0.\n    carbonate: 5.\n\n", "lysocline"),
        ("pelagic:\n  - start: 0.\n    carbonate: 5.\n    lysocline: 4000.\n"
         "    ccd: 3000.\n\n", "ccd >= lysocline"),
    ):
        with pytest.raises(ValueError, match=msg):
            Model(_case(tmp_path, "pb", block=block), verbose=False, showlog=False)
    with pytest.raises(ValueError, match="dual lithology"):
        Model(_case(tmp_path, "pd", block=PEL, dual=False), verbose=False, showlog=False)


def test_pelagic_settling_exact(tmp_path, monkeypatch):
    """
    Protects: the settled thickness. On a seafloor at three water depths
    (above the lysocline, half-way to the CCD, below the CCD) one step
    deposits ``rate * dt * (1, 0.5, 0)`` of carb_mud plus ``clay * dt`` of
    fine everywhere (20 m/Myr), exactly; a rate too small for the stratigraphic floor
    is carried over, not lost.
    """
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    m = Model(_case(tmp_path, "px", block=PEL), verbose=False, showlog=False)
    try:
        n = m.lpoints
        depth = np.array([1000.0, 4000.0, 5000.0])[np.arange(n) % 3]
        m.sealevel = 0.0
        h = -depth
        m.hLocal.setArray(h)
        m.dm.localToGlobal(m.hLocal, m.hGlobal)
        m.lFill = h.copy()
        m.tNow = 0.0
        top = m.stratStep
        c0 = m.stratHL[1][:, top].copy()
        f0 = m.stratHL[0][:, top].copy()
        m.pelagicRain()
        carb = m.stratHL[1][:, top] - c0
        clay = m.stratHL[0][:, top] - f0
        exp = 50.0e-6 * m.dt * np.array([1.0, 0.5, 0.0])[np.arange(n) % 3]
        assert np.allclose(carb, exp, rtol=1e-12, atol=1e-15)
        assert np.allclose(clay, 20.0e-6 * m.dt, rtol=1e-12)
        assert np.allclose(m.hLocal.getArray() - h, exp + 20.0e-6 * m.dt, rtol=1e-12)
        # a tiny rate: carried over until it reaches the 0.1 mm floor
        m.pelEvents[0]["carb"] = m._pelMaps[0]["carb"] = 1.0e-7
        m.pelEvents[0]["clay"] = m._pelMaps[0]["clay"] = 0.0
        before = m.stratHL[1][:, top].copy()
        nstep = 0
        while nstep < 200 and np.all(m.stratHL[1][::3, top] == before[::3]):
            m.pelagicRain()
            nstep += 1
        need = int(np.ceil(1.0e-4 / (1.0e-7 * m.dt)))
        assert nstep == need, (nstep, need)
        got = m.stratHL[1][::3, top] - before[::3]
        assert np.allclose(got, need * 1.0e-7 * m.dt, rtol=1e-9)
    finally:
        m.destroy()


@pytest.mark.slow
def test_pelagic_budget(tmp_path, monkeypatch):
    """Over a run: carbonate rain = settled + dissolved, the settled
    carbonate is the carb_mud of the stratigraphy (exactly, no carbonate
    block), and the output field exists."""
    from gospl.model import Model
    import h5py

    monkeypatch.chdir(tmp_path)
    m = Model(_case(tmp_path, "pg", block=PEL), verbose=False, showlog=False)
    try:
        m.runProcesses()
        own = m.inIDs == 1
        carry = (m._pelCarry[:, 0] * m.larea)[own].sum()
        assert m._pelRain > 0.0 and m._pelDissolved > 0.0 and m._pelClay > 0.0
        assert abs(m._pelRain - m._pelPreserved - m._pelDissolved) < 1e-9 * m._pelRain
        dep = dict(zip(m.litName, m._litDeposited))
        assert abs(dep["carb_mud"] + carry - m._pelPreserved) < 1e-9 * m._pelPreserved
    finally:
        m.destroy()
    with h5py.File(tmp_path / "pg_out" / "h5" / "gospl.1.p0.h5", "r") as f:
        assert "pelDep" in f


def test_carbonate_exposure_response(tmp_path, monkeypatch):
    """The open-water index is positive on highs and negative in lows of
    the seafloor, and its response table scales the factory accordingly."""
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    m = Model(_case(tmp_path, "e", block=CARB_EXP), verbose=False, showlog=False)
    try:
        h = m.hLocal.getArray().copy()
        st = m._carbStatic(h)
        E = m.carbExpo
        own = m.inIDs == 1
        # the highest and lowest local nodes sit above / below their smoothing
        assert E[own][np.argmax(h[own])] > 0.0 and E[own][np.argmin(h[own])] < 0.0
        exp = np.interp(E, [-50.0, 50.0], [0.2, 1.0])
        assert np.allclose(st[0], exp) and np.all(st[1] == 1.0)
    finally:
        m.destroy()


def test_pelagic_parse_c2b(tmp_path, monkeypatch):
    """clay_decay and resuspend are read per event; a non-positive decay
    length is refused."""
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    m = Model(_case(tmp_path, "q", block=PEL_RESUSP), verbose=False, showlog=False)
    try:
        e = m.pelEvents[0]
        assert e["decay"] == 3.0e5 and e["resuspend"] is True
        assert m._pelResuspend()
    finally:
        m.destroy()
    with pytest.raises(ValueError, match="clay_decay"):
        Model(_case(tmp_path, "qb", block=PEL.replace("    clay: 20.\n",
                                                      "    clay: 20.\n    clay_decay: 0.\n")),
              verbose=False, showlog=False)


@pytest.mark.slow
def test_pelagic_resuspension(tmp_path, monkeypatch):
    """
    Protects: the rain handed to the marine diffusion. Settled in place the
    carbonate mud is at most the rain of the run at a node, and none lies
    below the CCD; resuspended, the diffusion gathers it into the lows (its
    thickest cell several times the settled maximum), and the carb_mud of the
    stratigraphy is the rain handed over, to the marine floors and the open
    edges (1 %).
    """
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    out = {}
    for name, block in (("set", PEL_DECAY), ("res", PEL_RESUSP)):
        m = Model(_case(tmp_path, name, block=block), verbose=False, showlog=False)
        try:
            m.runProcesses()
            top = m.stratStep + 1
            cm = m.carbLit["carb_mud"]
            cmud = m.stratHL[cm][:, 1:top].sum(axis=1)
            deep = (m.sealevel - m.hLocal.getArray()) > 4600.0
            out[name] = (cmud.max(), cmud[deep].sum(), m._pelPreserved,
                         dict(zip(m.litName, m._litDeposited))["carb_mud"])
        finally:
            m.destroy()
    assert out["set"][1] == 0.0
    assert out["res"][0] > 5.0 * out["set"][0], (out["res"][0], out["set"][0])
    handed, laid = out["res"][2], out["res"][3]
    assert 0.99 * handed < laid <= handed * (1 + 1e-9), (laid, handed)


def test_sst_temperature(tmp_path, monkeypatch):
    """
    Protects: the sea-surface temperature (phase C5). Events are a step
    function, the curve is interpolated in time and added to the event map
    as an anomaly (or used alone); the ``temperature`` response scales a
    factory and the pelagic rain; temperature without an ``sst`` block is
    refused.
    """
    from gospl.model import Model

    monkeypatch.chdir(tmp_path)
    shutil.copy(FIXTURES_DIR / "mesh.npz", tmp_path / "mesh.npz")
    n = np.load(tmp_path / "mesh.npz")["v"].shape[0]
    np.savez(tmp_path / "sst0.npz", t=np.linspace(10.0, 30.0, n))
    np.savetxt(tmp_path / "anom.csv", np.array([[0.0, 0.0], [100.0, -4.0]]), delimiter=",")
    block = ("sst:\n    events:\n      - {start: 50., uniform: 20.}\n"
             "      - {start: 0., map: ['sst0', 't']}\n    curve: 'anom.csv'\n\n"
             "carbonate:\n    factories:\n      - {name: reef, growth: 500., "
             "response: [{control: temperature, x: [16., 18., 28., 32.], "
             "f: [0., 1., 1., 0.]}]}\n\n"
             "pelagic:\n  - start: 0.\n    carbonate: 50.\n    lysocline: 3500.\n"
             "    ccd: 4500.\n    temperature: {x: [10., 30.], f: [0., 1.]}\n\n")
    m = Model(_case(tmp_path, "sst", block=block), verbose=False, showlog=False)
    try:
        assert m.sstOn and [e["start"] for e in m.sstEvents] == [0.0, 50.0]
        t0 = np.asarray(m._sstMaps[0])
        m.tNow = 25.0                                   # map + (-1) anomaly
        assert np.allclose(m._sst(), t0 - 1.0)
        m.tNow = 75.0                                   # uniform 20 + (-3)
        assert np.allclose(m._sst(), 17.0)
        m.tNow = 200.0                                  # curve constant beyond its end
        assert np.allclose(m._sst(), 16.0)
        m.tNow = 0.0
        st = m._carbStatic(m.hLocal.getArray())
        assert np.allclose(st[0], np.interp(t0, [16., 18., 28., 32.], [0., 1., 1., 0.]))
        # pelagic rain scaled by the temperature table
        m.sealevel = 0.0
        h = np.full(m.lpoints, -1000.0)
        m.hLocal.setArray(h)
        m.dm.localToGlobal(m.hLocal, m.hGlobal)
        m.lFill = h.copy()
        top = m.stratStep
        c0 = m.stratHL[m.carbLit["carb_mud"]][:, top].copy()
        m.pelagicRain()
        got = m.stratHL[m.carbLit["carb_mud"]][:, top] - c0
        exp = 50.0e-6 * m.dt * np.interp(t0, [10.0, 30.0], [0.0, 1.0])
        lay = exp >= 1.0e-4
        assert np.allclose(got[lay], exp[lay], rtol=1e-12)
    finally:
        m.destroy()
    # the same factory without an sst block
    carb_only = "carbonate:" + block.split("carbonate:")[1].split("pelagic:")[0]
    with pytest.raises(ValueError, match="sst"):
        Model(_case(tmp_path, "nosst", block=carb_only), verbose=False, showlog=False)


def test_carbonate_provenance_marine_class(tmp_path, monkeypatch):
    """With provenance on, sediment made in the sea (carbonate framework and
    loose carbonate, pelagic carbonate and clay) gets one extra class after
    the user's: every layer stays exactly partitioned by class, the user
    classes keep their eroded/deposited balance (the carbonate used to be
    counted in them), and the marine class is deposited but never eroded."""
    from mpi4py import MPI
    from gospl.model import Model

    src = FIXTURES_DIR / "minimal_prov2.yml"
    if not src.exists() or not (FIXTURES_DIR / "mesh.npz").exists():
        pytest.skip("minimal_prov2.yml / mesh.npz not present")
    monkeypatch.chdir(tmp_path)
    shutil.copy(FIXTURES_DIR / "mesh.npz", tmp_path / "mesh.npz")
    v = np.load(tmp_path / "mesh.npz")["v"]
    np.savez(tmp_path / "prov_src.npz", rock=(v[:, 0] > v[:, 0].mean()).astype(np.int64))
    block = ("strata:\n    dual: True\n    coarse: {phi0: 0.49, z0: 3700.}\n"
             "    fine: {phi0: 0.63, z0: 1960., k_factor: 1.5}\n"
             "    bedrock_coarse_frac: 0.6\n\n"
             "carbonate:\n    factories:\n"
             "      - {name: reef, growth: 2000., extinction: 0.1, framework: 0.5, sand: 0.7}\n"
             "      - {name: deep, growth: 200., extinction: 0.005, sand: 0.3}\n\n"
             + PEL)
    text = re.sub(r"dir: '[^']*'", "dir: 'out'", src.read_text())
    (tmp_path / "p.yml").write_text(text.replace("output:", block + "output:", 1))
    m = Model("p.yml", verbose=False, showlog=False)
    try:
        assert m.provOn and m.provUserNb == 2 and m.provMarine == 2 and m.provNb == 3
        assert np.all(m.source_class < m.provUserNb)
        m.runProcesses()
        top = m.stratStep + 1
        H, P = m.stratH[:, :top], m.stratP[:, :top, :]
        assert np.abs(P.sum(axis=2) - H).sum() / max(H.sum(), 1.0) < 1.0e-9
        ero = MPI.COMM_WORLD.allreduce(m._provEroded, op=MPI.SUM)
        dep = MPI.COMM_WORLD.allreduce(m._provDeposited, op=MPI.SUM)
        assert ero[2] == 0.0 and dep[2] > 0.0
        for c in (0, 1):
            assert abs(dep[c] - ero[c]) < 1.0e-3 * ero[c]
        reef = m.stratHL[m.litName.index("reef")][:, :top]
        cell = (reef > 0.99 * H) & (H > 1.0e-3)
        if cell.any():
            np.testing.assert_allclose(P[:, :, 2][cell], H[cell], rtol=1.0e-9)
    finally:
        m.destroy()
