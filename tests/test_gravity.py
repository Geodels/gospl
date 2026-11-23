"""
Marine gravity currents (opt-in ``gravity:`` block), docs/DESIGN_GRAVITY_CURRENTS.md.

Guards: parsing and input checks; ``fraction: 0`` bitwise identical to no block;
provenance classes routed with the gravity pathway (phase G2);
the routed budget (supplied = deposited + exported); dual-lithology fine
conservation with the fine/coarse split between the two pathways; the routed
solve against an explicit cascade of the same deposition law; the flood
fraction against a Monte-Carlo estimate; the ``gravDep`` output.
"""

from __future__ import annotations

import glob
import os
import re
import shutil
import types

import numpy as np
import pytest

from _helpers import FIXTURES_DIR

pytest.importorskip(
    "gospl._fortran", reason="goSPL Fortran extension not built"
)

pytestmark = [pytest.mark.sediment]

GRAV_BLOCK = (
    "gravity:\n"
    "    runout: 3.0e5\n"
    "    runout_fine: 1.2e6\n"
    "    fraction: {frac}\n\n"
)


def _write_case(tmp_path, name, gravity=None, extra=""):
    """minimal_dual.yml copied to tmp_path, with an optional gravity block."""
    src = FIXTURES_DIR / "minimal_dual.yml"
    if not src.exists() or not (FIXTURES_DIR / "mesh.npz").exists():
        pytest.skip("minimal_dual.yml / mesh.npz not present")
    if not (tmp_path / "mesh.npz").exists():
        shutil.copy(FIXTURES_DIR / "mesh.npz", tmp_path / "mesh.npz")
    text = src.read_text()
    text = re.sub(r"dir: '[^']*'", f"dir: '{name}_out'", text)
    block = (gravity if gravity is not None else "") + extra
    if block:
        text = text.replace("output:", block + "output:", 1)
    (tmp_path / f"{name}.yml").write_text(text)
    return f"{name}.yml"


def _run(tmp_path, yml, attrs=()):
    from gospl.model import Model

    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        m = Model(yml, verbose=False, showlog=False)
        try:
            m.runProcesses()
            out = {a: getattr(m, a) for a in attrs}
            out["elev"] = m.hLocal.getArray().copy()
            out["flatModel"] = m.flatModel
        finally:
            m.destroy()
    finally:
        os.chdir(cwd)
    return out


def test_gravity_off_by_default(tmp_path):
    """No gravity block: the option is off and nothing is allocated."""
    from gospl.model import Model

    yml = _write_case(tmp_path, "off")
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        m = Model(yml, verbose=False, showlog=False)
        try:
            assert m.gravOn is False
            assert not hasattr(m, "gravDep")
        finally:
            m.destroy()
    finally:
        os.chdir(cwd)


@pytest.mark.parametrize(
    "block, msg",
    [
        ("gravity:\n    runout: -1.\n\n", "positive"),
        ("gravity:\n    fraction: 1.5\n\n", "between 0 and 1"),
        ("gravity:\n    flood: {rating_exp: 1.0}\n\n", "rating_exp"),
    ],
)
def test_gravity_input_checks(tmp_path, block, msg):
    """Invalid gravity inputs stop at parse time with a clear message."""
    from gospl.model import Model

    yml = _write_case(tmp_path, "bad", gravity=block)
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        with pytest.raises(ValueError, match=msg):
            Model(yml, verbose=False, showlog=False)
    finally:
        os.chdir(cwd)


@pytest.mark.slow
@pytest.mark.parametrize("dual", [False, True])
def test_gravity_provenance(tmp_path, dual):
    """Provenance classes ride the gravity pathway (phase G2): two classes from
    a map, half the load as gravity currents, single fraction and dual
    lithology. The per-class layer thicknesses partition the layer, the
    deposited class ratio tracks the eroded one, the gravity class deposits sum
    to the gravity deposit, and a node that received only a gravity deposit
    records the gravity composition (``_gravProvFraction``), not the
    domain-average fallback of ``_marineProvFraction``."""
    from mpi4py import MPI
    from gospl.model import Model

    src = FIXTURES_DIR / "minimal_prov2.yml"
    if not src.exists() or not (FIXTURES_DIR / "mesh.npz").exists():
        pytest.skip("minimal_prov2.yml / mesh.npz not present")
    shutil.copy(FIXTURES_DIR / "mesh.npz", tmp_path / "mesh.npz")
    v = np.load(tmp_path / "mesh.npz")["v"]
    np.savez(tmp_path / "prov_src.npz", rock=(v[:, 0] > v[:, 0].mean()).astype(np.int64))
    text = re.sub(r"dir: '[^']*'", "dir: 'out'", src.read_text())
    block = GRAV_BLOCK.format(frac=0.5)
    if dual:
        block += ("strata:\n    dual: True\n    coarse: {phi0: 0.49, z0: 3700.}\n"
                  "    fine: {phi0: 0.63, z0: 1960., k_factor: 1.5}\n"
                  "    bedrock_coarse_frac: 0.6\n\n")
    (tmp_path / "prov.yml").write_text(text.replace("output:", block + "output:", 1))
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        m = Model("prov.yml", verbose=False, showlog=False)
        try:
            if m.flatModel:
                pytest.skip("needs the closed sphere fixture")
            assert m.gravOn and m.provOn and m.stratLith == dual
            m.runProcesses()
            top = m.stratStep + 1
            H, P = m.stratH[:, :top], m.stratP[:, :top, :]
            assert np.abs(P.sum(axis=2) - H).sum() / max(H.sum(), 1.0) < 1.0e-6

            ero = MPI.COMM_WORLD.allreduce(m._provEroded, op=MPI.SUM)
            dep = MPI.COMM_WORLD.allreduce(m._provDeposited, op=MPI.SUM)
            assert ero.sum() > 0 and dep.sum() > 0
            assert abs(dep[1] / dep.sum() - ero[1] / ero.sum()) < 1.0e-2

            casc, gC, gF = m._gravThick
            gp = m._gravProvThick
            g = gC + gF
            hit = g > 1.0e-6 * g.max()
            assert hit.any()
            # separate solves (classes, total) agree to the KSP tolerance
            np.testing.assert_allclose(gp[hit].sum(axis=1), g[hit], rtol=0.0,
                                       atol=1.0e-6 * g.max())

            only = hit & (casc <= 0.0) & (m.inIDs == 1)
            assert only.any(), "no gravity-only node: blend check is vacuous"
            np.testing.assert_allclose(
                m.depoProvFrac[only], gp[only] / g[only, None], atol=1.0e-9
            )
        finally:
            m.destroy()
    finally:
        os.chdir(cwd)


@pytest.mark.slow
def test_gravity_zero_fraction_is_identical(tmp_path):
    """``fraction: 0`` diverts nothing, so elevations must equal the run
    without the block bitwise: the gravity path has no side effect."""
    ref = _run(tmp_path, _write_case(tmp_path, "ref"))
    zero = _run(tmp_path, _write_case(tmp_path, "zero", GRAV_BLOCK.format(frac=0.0)))
    np.testing.assert_array_equal(ref["elev"], zero["elev"])


@pytest.mark.slow
def test_gravity_budget_and_fine_conservation(tmp_path):
    """Half of the river load travels as gravity currents. The routed budget
    closes every step (supplied = deposited + exported, to the solver
    tolerance) and, with dual lithology, the fine volume eroded equals the fine
    volume deposited: the fine part of the gravity supply is removed from the
    cascade's fine budget and the node fine fraction blends both deposits."""
    from mpi4py import MPI

    out = _run(
        tmp_path, _write_case(tmp_path, "half", GRAV_BLOCK.format(frac=0.5)),
        attrs=("gravBudget", "gravDep", "_fineEroded", "_fineDeposited"),
    )
    if out["flatModel"]:
        pytest.skip("needs the closed sphere fixture")
    sup, dep, exp, ero = out["gravBudget"]
    assert sup > 0.0 and ero == 0.0           # no slope erosion configured
    assert abs(sup - dep - exp) / sup < 1.0e-6
    assert MPI.COMM_WORLD.allreduce(float(out["gravDep"].max()), op=MPI.MAX) > 0.0
    ero = MPI.COMM_WORLD.allreduce(out["_fineEroded"], op=MPI.SUM)
    fdep = MPI.COMM_WORLD.allreduce(out["_fineDeposited"], op=MPI.SUM)
    assert ero > 0.0
    assert abs(ero - fdep) / ero < 5.0e-3, f"fine {ero:.4e} vs {fdep:.4e}"


@pytest.mark.slow
def test_gravity_route_matches_explicit_cascade(tmp_path):
    """The routed solve (I - W^T diag(1 - Lambda)) q = s equals the explicit
    cascade of the same deposition law: deposit Lambda q at each node and pass
    (1 - Lambda) q to the receivers until nothing is left. Independent check of
    the operator assembly (column scaling, sinks, outlets)."""
    from gospl.model import Model

    yml = _write_case(tmp_path, "route", GRAV_BLOCK.format(frac=0.5))
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        m = Model(yml, verbose=False, showlog=False)
        try:
            m.runProcesses()          # builds the coast distance and pit state
            m.sinkIDs = m.lFill <= m.sealevel
            m._matOcean()
            m._gravNetwork(m.hLocal.getArray().copy())
            supply = np.zeros(m.lpoints)
            sea = np.where(m.sinkIDs & (m.inIDs == 1))[0]
            supply[sea[:: max(1, len(sea) // 20)]] = 1.0e9
            dep, exported = m._gravRoute(supply, m.gravRunout)

            lam = m._gravLambda(m.gravRunout)
            ref = np.zeros(m.lpoints)
            cur = supply.copy()
            for _ in range(10000):
                ref += lam * cur
                m.tmpL.setArray((1.0 - lam) * cur)
                m.dm.localToGlobal(m.tmpL, m.tmp)
                m._gW.mult(m.tmp, m.tmp1)
                m.dm.globalToLocal(m.tmp1, m.tmpL)
                cur = m.tmpL.getArray().copy()
                if cur.sum() < 1.0e-9 * supply.sum():
                    break
            m._gravRelease()
            m.dMat1.destroy()
            owned = m.inIDs == 1
            err = np.abs(dep - ref)[owned].sum() / supply[owned].sum()
            assert err < 1.0e-6, f"routed vs explicit deposit differ ({err:.2e})"
        finally:
            m.destroy()
    finally:
        os.chdir(cwd)


def test_gravity_flood_fraction_montecarlo():
    """The analytic flood fraction equals a Monte-Carlo estimate of the same
    model: log-normal daily discharge (mean 1, log-std sigma), rating curve
    Qs ~ Q^b, fraction of the load carried above the threshold concentration."""
    from gospl.sed.seaplex import SEAMesh

    stub = types.SimpleNamespace(
        gravRating=1.8, gravSigma=1.2, gravRhoS=2650.0, gravCcrit=40.0, lpoints=3
    )
    qs = np.array([1.0e6, 1.0e7, 0.0])          # m3/yr
    q = np.array([1.0e9, 1.0e9, 1.0e9])          # m3/yr -> Cbar 2.65, 26.5, 0 kg/m3
    frac = SEAMesh._gravFloodFraction(stub, qs, q)
    assert frac[2] == 0.0

    rng = np.random.default_rng(1)
    s, b = stub.gravSigma, stub.gravRating
    x = np.exp(rng.normal(-0.5 * s * s, s, 2_000_000))
    load = x ** b
    for k in range(2):
        cbar = stub.gravRhoS * qs[k] / q[k]
        conc = cbar * x ** (b - 1.0) / load.mean()
        mc = load[conc > stub.gravCcrit].sum() / load.sum()
        assert abs(frac[k] - mc) < 0.01, (k, frac[k], mc)
    assert frac[1] > frac[0] > 0.0


@pytest.mark.slow
def test_gravity_output_written(tmp_path):
    """``gravDep`` is written to every output step when the option is on."""
    import h5py

    _run(tmp_path, _write_case(tmp_path, "outg", GRAV_BLOCK.format(frac=0.5)))
    files = sorted(glob.glob(str(tmp_path / "outg_out" / "h5" / "gospl.*.p0.h5")))
    assert files
    with h5py.File(files[-1]) as h:
        assert "gravDep" in h
        assert h["gravDep"].shape[1] == 1
    with open(sorted(glob.glob(str(tmp_path / "outg_out" / "xmf" / "*.xmf")))[-1]) as f:
        assert "gravDep" in f.read()


@pytest.mark.slow
def test_gravity_restart_restores_deposit(tmp_path):
    """A run restarted from an output step recovers the cumulative gravity
    deposit written at that step (``outmesh.readData``)."""
    import h5py
    from gospl.model import Model

    _run(tmp_path, _write_case(tmp_path, "rst", GRAV_BLOCK.format(frac=0.5)))
    with h5py.File(tmp_path / "rst_out" / "h5" / "gospl.1.p0.h5") as h:
        saved = np.array(h["gravDep"])[:, 0]
    assert saved.max() > 0.0
    text = (tmp_path / "rst.yml").read_text()
    text = text.replace("    strat: 50\n", "    strat: 50\n    rstep: 1\n", 1)
    (tmp_path / "rst2.yml").write_text(text)
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        m = Model("rst2.yml", verbose=False, showlog=False)
        try:
            np.testing.assert_allclose(m.gravDep, saved, rtol=1e-6, atol=1e-6)
        finally:
            m.destroy()
    finally:
        os.chdir(cwd)


@pytest.mark.slow
def test_gravity_run_summary(tmp_path):
    """The ``--summary`` record of every step carries the gravity budget, which
    closes, and ``gospl-inspect`` reports no gravity issue; a tampered record
    is flagged."""
    import json
    from gospl.analyse.runinspect import inspect_summary
    from gospl.model import Model

    yml = _write_case(tmp_path, "summ", GRAV_BLOCK.format(frac=0.5))
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        m = Model(yml, verbose=False, showlog=False, summary="run.jsonl")
        try:
            m.runProcesses()
        finally:
            m.destroy()
    finally:
        os.chdir(cwd)
    recs = [json.loads(l) for l in open(tmp_path / "run.jsonl") if l.strip()]
    steps = [r for r in recs if "step" in r and "elev" in r]
    assert steps and all("gravity" in r for r in steps)
    g = steps[-1]["gravity"]
    assert g["supplied"] > 0.0
    assert abs(g["supplied"] - g["deposited"] - g["exported"]) / g["supplied"] < 1e-6
    rep = inspect_summary(str(tmp_path / "run.jsonl"))
    assert not [i for i in rep["issues"] if i["issue"] == "gravity_budget"]

    steps[-1]["gravity"]["deposited"] *= 0.9
    bad = tmp_path / "bad.jsonl"
    bad.write_text("".join(json.dumps(r) + "\n" for r in recs))
    rep = inspect_summary(str(bad))
    assert [i for i in rep["issues"] if i["issue"] == "gravity_budget"]


@pytest.mark.slow
def test_gravity_accommodation_cap(tmp_path):
    """The accommodation cap of ``_gravRoute``: when each node has room for a
    tenth of what it received uncapped, no non-sink node exceeds its room, the
    cap is active, the budget still closes (supplied = deposited + exported)
    and the capped volume moved on downstream instead of vanishing."""
    from gospl.model import Model

    yml = _write_case(tmp_path, "cap", GRAV_BLOCK.format(frac=0.5))
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        m = Model(yml, verbose=False, showlog=False)
        try:
            m.runProcesses()
            m.sinkIDs = m.lFill <= m.sealevel
            m._matOcean()
            m._gravNetwork(m.hLocal.getArray().copy())
            supply = np.zeros(m.lpoints)
            sea = np.where(m.sinkIDs & (m.inIDs == 1))[0]
            supply[sea[:: max(1, len(sea) // 20)]] = 1.0e9
            free, _ = m._gravRoute(supply, m.gravRunout)
            room = 0.1 * free
            dep, exp = m._gravRoute(supply, m.gravRunout, room=room)
            m._gravRelease()
            m.dMat1.destroy()
            owned = m.inIDs == 1
            ok = owned & ~m._gSink
            assert (dep[ok] <= room[ok] * (1.0 + 1.0e-6) + 1.0e-3).all()
            assert (free[ok] > room[ok] + 1.0).any(), "cap never active"
            sup = supply[owned].sum()
            assert abs(sup - dep[owned].sum() - exp[owned].sum()) / sup < 1.0e-6
            # the capped volume went elsewhere (sinks, other nodes or export)
            assert (dep[ok] < free[ok] - 1.0).any()
        finally:
            m.destroy()
    finally:
        os.chdir(cwd)


def test_gravity_raw_network_conserves(tmp_path):
    """The raw gravity network (``mfdrcvrs`` on the eps-filled bathymetry with
    6 receivers): every non-sink node passes ALL its flux on (outgoing weights
    sum to 1; mfdrcvrs pads unused slots with -1 and does not renormalise),
    a sink passes nothing, and with almost no deposition on the way the routed
    volume closes on sinks + export (an uphill receiver at a pit below the fill
    cut-off made a cycle and a singular solve)."""
    from gospl.model import Model

    yml = _write_case(tmp_path, "rawnet", GRAV_BLOCK.format(frac=0.5))
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        m = Model(yml, verbose=False, showlog=False)
        try:
            m.runProcesses()
            m.sinkIDs = m.lFill <= m.sealevel
            m._matOcean()
            m._gravNetwork(m.hLocal.getArray().copy())
            # outgoing weight of each node: column sums of W^T
            m.tmpL.set(1.0)
            m.dm.localToGlobal(m.tmpL, m.tmp)
            m._gW.multTranspose(m.tmp, m.tmp1)
            m.dm.globalToLocal(m.tmp1, m.tmpL)
            wout = m.tmpL.getArray().copy()
            owned = m.inIDs == 1
            inner = owned.copy()
            if m.flatModel:
                inner[m.outletIDs] = False
            ns = inner & ~m._gSink
            assert ns.any() and m._gSink[owned].any()
            assert np.abs(wout[ns] - 1.0).max() < 1.0e-12, "flux lost on the way"
            assert np.abs(wout[inner & m._gSink]).max() < 1.0e-12
            supply = np.zeros(m.lpoints)
            sea = np.where(m.sinkIDs & owned)[0]
            supply[sea[:: max(1, len(sea) // 20)]] = 1.0e9
            dep, exp = m._gravRoute(supply, 1.0e15)
            m._gravRelease()
            m.dMat1.destroy()
            sup = supply[owned].sum()
            err = abs(sup - dep[owned].sum() - exp[owned].sum()) / sup
            assert err < 1.0e-6, f"routed budget {err:.2e}"
            # nearly everything reached a sink or left the domain
            through = (dep[owned & m._gSink].sum() + exp[owned].sum()) / sup
            assert through > 0.5, f"only {through:.2f} reached sinks/export"
        finally:
            m.destroy()
    finally:
        os.chdir(cwd)


EROSION_BLOCK = (
    "gravity:\n"
    "    runout: 3.0e5\n"
    "    runout_fine: 1.2e6\n"
    "    fraction: 0.5\n"
    "    erosion_length: 5.0e4\n"
    "    erosion_slope: 0.005\n"
    "    max_erosion: 2.0\n\n"
)


@pytest.mark.slow
def test_gravity_slope_erosion(tmp_path):
    """Slope erosion by gravity currents (phase G3b): the currents pick up
    sediment, the budget closes with it (supplied + eroded = deposited +
    exported), the fine volume is conserved with the eroded fine routed on,
    and in the last step no node eroded more than ``max_erosion``, none above
    the bypass depth, and no node both eroded and received a deposit."""
    from mpi4py import MPI
    from gospl.model import Model
    from gospl.sed.seaplex import SEAMesh

    rec = {}
    orig = SEAMesh._gravErosion

    def spy(self, supply, hl):
        e, mask = orig(self, supply, hl)
        rec.update(e=e.copy(), mask=mask.copy(), hl=hl.copy())
        return e, mask

    yml = _write_case(tmp_path, "ero", gravity=EROSION_BLOCK)
    cwd = os.getcwd()
    os.chdir(tmp_path)
    SEAMesh._gravErosion = spy
    try:
        m = Model(yml, verbose=False, showlog=False)
        try:
            m.runProcesses()
            sup, dep, exp, ero = m.gravBudget
            assert ero > 0.0, "no slope erosion: test is vacuous"
            assert abs(sup + ero - dep - exp) / (sup + ero) < 1.0e-5
            fe = MPI.COMM_WORLD.allreduce(m._fineEroded, op=MPI.SUM)
            fd = MPI.COMM_WORLD.allreduce(m._fineDeposited, op=MPI.SUM)
            assert abs(fe - fd) / fe < 5.0e-3, f"fine {fe:.4e} vs {fd:.4e}"

            owned = m.inIDs == 1
            eTh = rec["e"] / m.larea
            assert (eTh[owned] <= m.gravMaxErosion * (1 + 1e-6)).all()
            shallow = (m.sealevel - rec["hl"]) <= m.gravBypassDepth
            assert (eTh[owned & shallow] == 0.0).all()
            casc, gC, gF = m._gravThick
            assert ((gC + gF)[rec["mask"] & owned] == 0.0).all()
        finally:
            m.destroy()
    finally:
        SEAMesh._gravErosion = orig
        os.chdir(cwd)


@pytest.mark.slow
def test_gravity_slope_erosion_provenance(tmp_path):
    """Eroded material carries its source classes: with slope erosion and two
    provenance classes, every layer stays exactly partitioned by class."""
    from gospl.model import Model

    src = FIXTURES_DIR / "minimal_prov2.yml"
    if not src.exists() or not (FIXTURES_DIR / "mesh.npz").exists():
        pytest.skip("minimal_prov2.yml / mesh.npz not present")
    shutil.copy(FIXTURES_DIR / "mesh.npz", tmp_path / "mesh.npz")
    v = np.load(tmp_path / "mesh.npz")["v"]
    np.savez(tmp_path / "prov_src.npz", rock=(v[:, 0] > v[:, 0].mean()).astype(np.int64))
    text = re.sub(r"dir: '[^']*'", "dir: 'out'", src.read_text())
    (tmp_path / "prov.yml").write_text(text.replace("output:", EROSION_BLOCK + "output:", 1))
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        m = Model("prov.yml", verbose=False, showlog=False)
        try:
            m.runProcesses()
            assert m.gravBudget[3] > 0.0, "no slope erosion: test is vacuous"
            top = m.stratStep + 1
            H, P = m.stratH[:, :top], m.stratP[:, :top, :]
            assert np.abs(P.sum(axis=2) - H).sum() / max(H.sum(), 1.0) < 1.0e-6
        finally:
            m.destroy()
    finally:
        os.chdir(cwd)


BELT_BLOCK = (
    "gravity:\n"
    "    runout: 3.0e5\n"
    "    runout_fine: 1.2e6\n"
    "    fraction: 0.5\n"
    "    network_smooth: 2\n"
    "    slope_smooth: 2\n"
    "    receivers: 2\n"
    "    flow_exp: 4.\n"
    "    receivers_fine: 6\n"
    "    flow_exp_fine: 1.\n"
    "    deposit_spread: {spread}\n"
    "    coarse_factor: {cf}\n\n"
)


@pytest.mark.slow
def test_gravity_channel_belt_options_conserve(tmp_path):
    """The channel-belt options (smoothed routing surface and slope, separate
    coarse and fine networks, lateral spreading of the deposit, coarse factor):
    the routed budget closes and the fine volume is conserved, as without
    them. Guards the network swap between the fractions and the spreading."""
    from mpi4py import MPI

    out = _run(
        tmp_path,
        _write_case(tmp_path, "belt", BELT_BLOCK.format(spread=4.0e5, cf=0.3)),
        attrs=("gravBudget", "_fineEroded", "_fineDeposited"),
    )
    if out["flatModel"]:
        pytest.skip("needs the closed sphere fixture")
    sup, dep, exp, ero = out["gravBudget"]
    assert sup > 0.0
    assert abs(sup - dep - exp) / sup < 1.0e-6
    fe = MPI.COMM_WORLD.allreduce(out["_fineEroded"], op=MPI.SUM)
    fd = MPI.COMM_WORLD.allreduce(out["_fineDeposited"], op=MPI.SUM)
    assert abs(fe - fd) / fe < 5.0e-3, f"fine {fe:.4e} vs {fd:.4e}"


def test_gravity_spread_conserves_and_keeps_uniform(tmp_path):
    """``_gravSpread`` conserves volume exactly, leaves a uniform thickness
    unchanged on the irregular mesh (the donor-share operator it replaced did
    not, and grew chains of dots at odd-sized cells), and smooths a single
    spike; ``_gravSpreadPasses`` turns the length into passes."""
    from gospl.model import Model

    yml = _write_case(tmp_path, "spread", BELT_BLOCK.format(spread=4.0e5, cf=1.0))
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        m = Model(yml, verbose=False, showlog=False)
        try:
            hl = np.full(m.lpoints, m.sealevel - 1000.0)  # every cell deep
            owned = m.inIDs == 1
            ok = np.ones(m.lpoints, bool)
            ok[m.outletIDs] = False
            uni = 3.0 * m.larea
            out = m._gravSpread(uni, 5, hl)
            th = out / m.larea
            assert np.abs(th[owned & ok] - 3.0).max() < 1.0e-10, "uniform layer changed"
            rng = np.random.default_rng(1)
            v = rng.random(m.lpoints) * m.larea
            v2 = m._gravSync(v)
            out = m._gravSpread(v2, 5, hl)
            assert abs(out[owned].sum() - v2[owned].sum()) / v2[owned].sum() < 1.0e-12
            spike = np.zeros(m.lpoints)
            i = np.where(owned & ok)[0][0]
            spike[i] = m.larea[i]
            out = m._gravSpread(m._gravSync(spike), 3, hl)
            assert out[i] < 0.9 * spike[i] and (out >= -1.0e-12).all()
            cell = np.sqrt(m.larea[owned]).mean()
            assert m._gravSpreadPasses() == int(np.ceil(2.0 * (4.0e5 / cell) ** 2))
        finally:
            m.destroy()
    finally:
        os.chdir(cwd)


def test_gravity_coarse_factor_scales_coarse_supply(tmp_path):
    """``coarse_factor`` diverts that share of the fraction for the coarse
    load only: on the first step (same state) the routed coarse supply scales
    by the factor and the fine supply is unchanged."""
    import gospl.sed.seaplex as sp
    from mpi4py import MPI

    seen = {}
    orig = sp.SEAMesh._gravRoute

    def spy(self, supply, runout, room=None, nodep=None):
        key = "c" if runout == self.gravRunout else "f"
        s0 = supply[:, 0] if np.ndim(supply) > 1 else supply
        if key not in seen:
            seen[key] = MPI.COMM_WORLD.allreduce(float(s0[self.inIDs == 1].sum()))
        return orig(self, supply, runout, room, nodep)

    res = {}
    try:
        sp.SEAMesh._gravRoute = spy
        for cf in (1.0, 0.25):
            seen.clear()
            _run(tmp_path, _write_case(tmp_path, f"cf{int(cf*100)}",
                                       BELT_BLOCK.format(spread=0.0, cf=cf)))
            res[cf] = dict(seen)
    finally:
        sp.SEAMesh._gravRoute = orig
    if "f" not in res[1.0]:
        pytest.skip("fixture without dual lithology")
    assert res[1.0]["c"] > 0.0
    assert abs(res[0.25]["c"] / res[1.0]["c"] - 0.25) < 1.0e-6
    assert abs(res[0.25]["f"] / res[1.0]["f"] - 1.0) < 1.0e-9


@pytest.mark.slow
def test_gravity_erosion_length_below_cell_is_capped(tmp_path):
    """An ``erosion_length`` shorter than the cell size gives a pick-up
    fraction above 1 per cell; uncapped, the flux grew as prod(1 + eps) along a
    path, every erosion solve failed and was discarded, and nothing eroded.
    With ``eps <= 1`` the erosion acts and the budget closes."""
    block = (
        "gravity:\n"
        "    runout: 3.0e6\n"
        "    runout_fine: 1.2e7\n"
        "    fraction: 0.5\n"
        "    erosion_length: 1.0e4\n"
        "    erosion_slope: 0.002\n"
        "    max_erosion: 2.0\n\n"
    )
    out = _run(tmp_path, _write_case(tmp_path, "short", block), attrs=("gravBudget",))
    if out["flatModel"]:
        pytest.skip("needs the closed sphere fixture")
    sup, dep, exp, ero = out["gravBudget"]
    assert ero > 0.0, "no erosion: the erosion solve failed"
    assert abs(sup + ero - dep - exp) / (sup + ero) < 1.0e-6
