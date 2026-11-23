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
    sup, dep, exp = out["gravBudget"]
    assert sup > 0.0
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
                m.dMat1.mult(m.tmp, m.tmp1)
                m.dm.globalToLocal(m.tmp1, m.tmpL)
                cur = m.tmpL.getArray().copy()
                if cur.sum() < 1.0e-9 * supply.sum():
                    break
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
            supply = np.zeros(m.lpoints)
            sea = np.where(m.sinkIDs & (m.inIDs == 1))[0]
            supply[sea[:: max(1, len(sea) // 20)]] = 1.0e9
            free, _ = m._gravRoute(supply, m.gravRunout)
            room = 0.1 * free
            dep, exp = m._gravRoute(supply, m.gravRunout, room=room)
            m.dMat1.destroy()
            owned = m.inIDs == 1
            ok = owned & ~m.is_sink_local
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
