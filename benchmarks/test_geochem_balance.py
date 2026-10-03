"""
Benchmarks for the Level-B groundwater geochemistry
(`gospl.flow.gwplex._updateSolute`, `_soluteAdvecCoeffs`).

Per step and per species, `_updateSolute` dissolves a source `D` (debiting the
source pool), solves the STEADY upwind advection-reaction problem

    [ upwind div(q .) + s + p ] c = D,      q = -T grad h

with `s = max(R - div q, 0)` the vertical discharge (seepage) sink, which
includes the recharge `R` at every node, and `p = precip_rate * Phi *
(1 - duriH/duriH_max)` the capillary-fringe precipitation sink, then books
`precipitated = p c A dt` and `exported = s c A dt`.

1. `test_geochem_budget_closure` (no closed form; an exact invariant): a 2-D
   dome with valleys, a real water table that seeps out in the valleys, a
   marine margin, erosion on, three species with different weatherability /
   precip_rate. Over 30 full model steps, every single step must satisfy, per
   species, `dissolved = precipitated + exported` (machine precision; the
   upwind face fluxes are antisymmetric and the solve is a direct LU), the
   pool debit must equal the cumulative dissolution, a zero-weatherability
   species must stay identically zero and a zero-precip species must export
   exactly what it dissolves.

2. `test_geochem_strip_analytic` (closed form): the Dupuit strip of
   `benchmarks/test_dupuit.py` (west drain, N/E/S walls, flat base, uniform
   recharge R). At steady state the discharge per unit width is
   `Q(u) = R u`, with `u = W - x` the distance from the east divide, and the
   discrete head balance gives `div q = R` at interior nodes, so the seepage
   sink `s` vanishes everywhere except at the drain column. With a uniform
   precipitation rate `p` (fringe favourability Phi = 1 via a very wide fringe,
   no crust feedback via a huge `max_thickness`) and a dissolution source `D0`
   confined to a band `u < Ls` next to the divide (per-vertex weatherability
   map), the 1-D steady balance is

       d(R u c)/du + p c = D(u)

   whose solution regular at the divide is

       c(u) = D0/(R+p)                              for u <= Ls
       c(u) = D0/(R+p) * (Ls/u)^(1 + p/R)           for u >  Ls.

   So `p = 0` gives the conservative dilution `c ~ 1/u` and `p = R` a `1/u^2`
   decay. The model is compared with BOTH this continuum solution (limited by
   first-order upwinding, O(dx/u)) and the exact solution of the same
   first-order upwind finite-volume recurrence on the strip, which it must
   reproduce to the head's distance from steady state. The exported fraction
   of the dissolved load (the drain column discharge `R W c(W)` over the band
   dissolution `D0 Ls`), the per-step budget and a zero-weatherability
   species are checked as well.
"""
import numpy as np
import pytest

scipy = pytest.importorskip("scipy")
from scipy.spatial import Delaunay  # noqa: E402

from mpi4py import MPI  # noqa: E402
from petsc4py import PETSc  # noqa: E402
from gospl.model import Model  # noqa: E402


def _rank0():
    return PETSc.COMM_WORLD.getRank() == 0


def _wrap_budget(m):
    """Record per-call increments of the geochem budget counters."""
    rec = []
    orig = m._updateSolute

    def wrapped():
        d0 = m.gwDissolved.copy()
        p0 = m.gwPrecip.copy()
        e0 = m.gwOceanFlux.copy()
        orig()
        rec.append((m.gwDissolved - d0, m.gwPrecip - p0, m.gwOceanFlux - e0))

    m._updateSolute = wrapped
    return rec


def _reduced(arr):
    return np.asarray(MPI.COMM_WORLD.allreduce(np.asarray(arr), op=MPI.SUM))


# =============================================================================
# 1. Budget closure on a 2-D dome with a marine margin
# =============================================================================
BX, BY, BDX = 61, 61, 250.0
N_STEPS_BUDGET = 30
TOL_CLOSE = 1.0e-10      # per-step |diss - precip - export| / diss


def _budget_mesh(path):
    xs = np.arange(BX) * BDX
    ys = np.arange(BY) * BDX
    X, Y = np.meshgrid(xs, ys)
    x, y = X.ravel(), Y.ravel()
    Lx = xs[-1]
    # A dome cut by radial valleys, tilted towards a sea along the west edge.
    xc, yc = 0.6 * Lx, 0.5 * Lx
    r = np.hypot(x - xc, y - yc)
    th = np.arctan2(y - yc, x - xc)
    z = (
        250.0 * np.exp(-((r / (0.35 * Lx)) ** 2))
        * (1.0 - 0.25 * np.cos(5.0 * th) ** 8)
        + 60.0 * (x / Lx) - 20.0
    )
    rng = np.random.default_rng(7)
    z += rng.normal(0.0, 0.5, z.size)
    v = np.column_stack([x, y, np.zeros_like(x)])
    cells = Delaunay(np.column_stack([x, y])).simplices.astype(np.int64)
    np.savez(path, v=v, c=cells, z=z.astype(np.float64))


def _budget_yaml(path):
    with open(path, "w") as f:
        f.write(
            "name: geochem budget closure benchmark\n"
            "domain:\n"
            "    npdata: ['budget_mesh','v','c','z']\n"
            "    flowdir: 2\n"
            "    bc: 'oooo'\n"
            "    seadepo: True\n"
            "    nodep: False\n"
            "time:\n"
            f"    start: 0.\n    end: {100.0 * N_STEPS_BUDGET}\n"
            f"    tout: {100.0 * N_STEPS_BUDGET}\n    dt: 100.\n"
            "spl:\n    K: 2.0e-6\n    d: 0.\n    m: 0.5\n"
            "diffusion:\n    hillslopeKa: 0.01\n    hillslopeKm: 0.1\n"
            "groundwater:\n"
            "    Ksat: 3650.\n"
            "    specific_yield: 0.1\n"
            "    aquifer_base: 40.\n"
            "    infiltration: 0.2\n"
            "    conserve_baseflow: True\n"
            "    duricrust:\n"
            "        max_thickness: 5.0\n"
            "        fringe_depth: 3.0\n"
            "        fringe_width: 4.0\n"
            "        armor_max: 0.5\n"
            "    geochem:\n"
            "        species:\n"
            "          - name: carbonate\n"
            "            weatherability: 1.0\n"
            "            precip_rate: 0.05\n"
            "            solid_volume: 0.01\n"
            "          - name: silica\n"
            "            weatherability: 0.3\n"
            "            precip_rate: 0.0\n"
            "          - name: inert\n"
            "            weatherability: 0.0\n"
            "            precip_rate: 0.5\n"
            "sea:\n    position: 0.\n"
            "climate:\n  - start: 0.\n    uniform: 1.0\n"
            "output:\n    dir: 'budget_out'\n"
        )


@pytest.mark.benchmark
@pytest.mark.slow
def test_geochem_budget_closure(tmp_path, monkeypatch):
    """Per-step, per-species dissolved = precipitated + exported (30 steps)."""
    monkeypatch.chdir(tmp_path)
    _budget_mesh(tmp_path / "budget_mesh.npz")
    _budget_yaml(tmp_path / "budget.yml")
    m = Model("budget.yml", verbose=False, showlog=False)
    try:
        assert m.gwGeochemOn and m.duriOn and m.gwNspecies == 3
        pool0 = m.gwSourcePool.copy()
        rec = _wrap_budget(m)
        m.runProcesses()
        assert len(rec) == N_STEPS_BUDGET, f"{len(rec)} solute updates"

        owned = m.inIDs == 1
        worst = np.zeros(3)
        for diss, prec, expo in rec:
            diss, prec, expo = _reduced(diss), _reduced(prec), _reduced(expo)
            for k in range(2):
                assert diss[k] > 0.0, "active species dissolved nothing"
                worst[k] = max(worst[k],
                               abs(diss[k] - prec[k] - expo[k]) / diss[k])
            # zero weatherability -> identically nothing
            assert diss[2] == 0.0 and prec[2] == 0.0 and expo[2] == 0.0
        cum_d = _reduced(m.gwDissolved)
        cum_p = _reduced(m.gwPrecip)
        cum_e = _reduced(m.gwOceanFlux)
        debit = _reduced(((pool0 - m.gwSourcePool)[owned]).sum(axis=0))
        pool_err = np.abs(debit - cum_d)[:2] / cum_d[:2]

        # Water table actually exists and seeps (non-trivial transport):
        wt = m.wtDepth[owned]
        sub = m._subaerialMask()[owned]
        nsat = MPI.COMM_WORLD.allreduce(int((wt[sub] < 0.5).sum()), op=MPI.SUM)
        nsub = MPI.COMM_WORLD.allreduce(int(sub.sum()), op=MPI.SUM)
        crust = MPI.COMM_WORLD.allreduce(
            float((m.duriHL.getArray()[owned] > 0).sum()), op=MPI.SUM)
        if _rank0():
            print(
                f"\n[geochem-budget] steps={len(rec)} subaerial={nsub} "
                f"seeping={nsat} crust nodes={crust:.0f}\n"
                f"  cumulative dissolved={cum_d}  precip={cum_p}  "
                f"export={cum_e}\n"
                f"  worst per-step closure: carbonate {worst[0]:.2e}, "
                f"silica {worst[1]:.2e}; pool-debit err {pool_err}"
            )
        assert np.all(worst[:2] < TOL_CLOSE), f"budget not closed: {worst}"
        assert np.all(pool_err < 1.0e-12), f"pool debit != dissolved: {pool_err}"
        # silica: no precipitation -> export == dissolution exactly
        assert cum_p[1] == 0.0
        assert abs(cum_e[1] - cum_d[1]) / cum_d[1] < TOL_CLOSE
        # carbonate actually precipitates (the sink is exercised)
        assert cum_p[0] > 1.0e-3 * cum_d[0], "no precipitation exercised"
        assert 0 < nsat < nsub, "need both seeping and unsaturated land"
        assert np.isfinite(m.gwSolute).all() and (m.gwSolute >= 0).all()
        assert (m.gwSolute[:, 2] == 0.0).all()
    finally:
        m.destroy()


# =============================================================================
# 2. Closed-form strip: Dupuit aquifer + band source + linear precipitation
# =============================================================================
NX, NY = 121, 9
DX = 100.0
W = (NX - 1) * DX                # 12 km
Z_DRAIN = 50.0
SLOPE = 0.02
KSAT = 3650.0
RAIN = 1.0
INFIL = 0.1
R = INFIL * RAIN                 # 0.1 m/yr
L0 = 2000.0                      # source band: nodes with W - x <= L0
LS = L0 + 0.5 * DX               # band edge = cell face
D0 = RAIN                        # proxy weathering supply = net rain (exp 1)
SPECIES = (                      # (name, weatherability map?, precip_rate)
    ("pR", True, R),             # p = R   -> c ~ u^-2 downstream of the band
    ("p0", True, 0.0),           # p = 0   -> conservative dilution c ~ 1/u
    ("off", False, 0.5),         # weatherability 0 -> c == 0
)
N_RELAX = 300
# Measured (serial): max dev vs discrete ~2e-5; RMSE vs continuum ~7e-3 (the
# first-order upwind error); exported fraction vs discrete ~1e-5.
TOL_DISCRETE = 1.0e-3            # max rel. dev from the discrete recurrence
TOL_CONT = 2.0e-2                # RMSE / c_band vs the continuum solution


def _strip_mesh(path):
    xs = np.arange(NX) * DX
    ys = np.arange(NY) * DX
    X, Y = np.meshgrid(xs, ys)
    x, y = X.ravel(), Y.ravel()
    v = np.column_stack([x, y, np.zeros_like(x)])
    cells = Delaunay(np.column_stack([x, y])).simplices.astype(np.int64)
    z = (Z_DRAIN + SLOPE * x).astype(np.float64)
    wsrc = ((W - x) <= L0 + 1.0e-6).astype(np.float64)
    np.savez(path, v=v, c=cells, z=z, abase=z.copy(), wsrc=wsrc)


def _strip_yaml(path):
    sp = ""
    for name, mapped, prate in SPECIES:
        wab = "['strip_mesh','wsrc']" if mapped else "0.0"
        sp += (
            f"          - name: {name}\n"
            f"            weatherability: {wab}\n"
            f"            precip_rate: {prate}\n"
            "            solid_volume: 1.0e-12\n"
        )
    with open(path, "w") as f:
        f.write(
            "name: geochem strip analytic benchmark\n"
            "domain:\n"
            "    npdata: ['strip_mesh','v','c','z']\n"
            "    flowdir: 2\n"
            "    bc: 'wwwf'\n"
            "    seadepo: False\n"
            "    nodep: True\n"
            "time:\n"
            "    start: 0.\n    end: 10.\n    tout: 10.\n    dt: 10.\n"
            "spl:\n    K: 1.0e-20\n    d: 0.\n    m: 0.5\n"
            "diffusion:\n    hillslopeKa: 0.\n    hillslopeKm: 0.\n"
            "groundwater:\n"
            f"    Ksat: {KSAT}\n"
            "    specific_yield: 0.1\n"
            "    aquifer_base: ['strip_mesh','abase']\n"
            f"    infiltration: {INFIL}\n"
            "    conserve_baseflow: False\n"
            "    picard_its: 3\n"
            "    seepage_passes: 4\n"
            "    duricrust:\n"
            # Phi = exp(-((wt-d0)/w)^2) ~ 1 everywhere; no crust feedback.
            "        fringe_depth: 0.\n"
            "        fringe_width: 1.0e7\n"
            "        max_thickness: 1.0e12\n"
            "        armor_max: 0.\n"
            "        break_rate: 0.\n"
            "        decay_rate: 0.\n"
            "    geochem:\n"
            "        species:\n" + sp +
            "sea:\n    position: -1000.\n"
            "climate:\n  - start: 0.\n"
            f"    uniform: {RAIN}\n"
            "output:\n    dir: 'strip_out'\n"
        )


def _continuum(u, p):
    c = np.full_like(u, D0 / (R + p))
    far = u > LS
    c[far] *= (LS / u[far]) ** (1.0 + p / R)
    return c


def _discrete(p):
    """Exact solution of the 1-D first-order upwind FV recurrence on the strip.

    Column j sits at u_j = j*DX from the divide (j = 0 .. NX-1, the drain is
    j = NX-1). Cell widths are DX except the half cells at both ends. At
    steady state the face discharge (per unit width) at u is R*u.
    """
    u = np.arange(NX) * DX
    src = (u <= L0 + 1.0e-6) * D0
    c = np.zeros(NX)
    c[0] = src[0] / (R + p)                    # divide half cell
    for j in range(1, NX - 1):
        qin, qout = R * (u[j] - 0.5 * DX), R * (u[j] + 0.5 * DX)
        c[j] = (qin * c[j - 1] + src[j] * DX) / (qout + p * DX)
    j = NX - 1                                  # drain half cell: sink R + Qin/A
    qin, a = R * (u[j] - 0.5 * DX), 0.5 * DX
    c[j] = (qin * c[j - 1] + src[j] * a) / (R * a + qin + p * a)
    return u, c


def _figure(curves, outdir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    for name, (u, cm, cc) in curves.items():
        ax.plot(u / 1e3, cm, "o", ms=3, label=f"goSPL {name}")
        ax.plot(u / 1e3, cc, "k-", lw=1)
    ax.set_yscale("log")
    ax.set_xlabel("distance from divide u (km)")
    ax.set_ylabel("solute concentration")
    ax.set_title("strip aquifer: band source + linear precipitation")
    ax.legend(frameon=False)
    fig.tight_layout()
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / "geochem_strip.png", dpi=120)
    plt.close(fig)


@pytest.mark.benchmark
@pytest.mark.slow
def test_geochem_strip_analytic(tmp_path, monkeypatch, request):
    """Steady solute along a Dupuit strip vs the closed-form 1-D solution."""
    from pathlib import Path

    if PETSc.COMM_WORLD.getSize() > 1:
        pytest.skip("1-D column comparison is written for a serial run")
    monkeypatch.chdir(tmp_path)
    _strip_mesh(tmp_path / "strip_mesh.npz")
    _strip_yaml(tmp_path / "strip.yml")
    m = Model("strip.yml", verbose=False, showlog=False)
    try:
        assert m.gwGeochemOn and m.duriOn and m.flatModel
        rec = _wrap_budget(m)
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()
        for _ in range(N_RELAX):
            m.updateGroundwater()

        x = m.lcoords[:, 0]
        u_node = W - x
        col = np.rint(u_node / DX).astype(int)
        z = m.hLocal.getArray()
        h = m.headL.getArray()
        assert len(m.seaID) == 0
        assert (h[x > 0.5 * DX] < z[x > 0.5 * DX]).all(), "table hit surface"

        curves = {}
        msgs = []
        for k, (name, mapped, prate) in enumerate(SPECIES):
            c = m.gwSolute[:, k]
            if not mapped:
                assert (c == 0.0).all(), f"{name}: zero-source species not zero"
                continue
            # 2-D strip must be 1-D: identical across the rows of a column.
            cm = np.array([c[col == j].mean() for j in range(NX)])
            spread = max(np.ptp(c[col == j]) for j in range(NX)) / cm.max()
            ud, cd = _discrete(prate)
            uc = ud
            cc = _continuum(uc, prate)
            dev_d = float(np.max(np.abs(cm - cd) / cd))
            rmse_c = float(np.sqrt(np.mean((cm - cc) ** 2)) / (D0 / (R + prate)))
            # export at the drain column vs the budget
            curves[name] = (uc, cm, cc)
            msgs.append(
                f"  {name}: p/R={prate / R:.1f} c_band={cm[0]:.4f} "
                f"(exact {D0 / (R + prate):.4f})  c_drain={cm[-1]:.4e} "
                f"(continuum {cc[-1]:.4e}, discrete {cd[-1]:.4e})  "
                f"max dev vs discrete={dev_d:.2e}  RMSE vs continuum="
                f"{rmse_c:.2e}  row spread={spread:.1e}"
            )
            assert spread < 1.0e-6, f"{name}: solution not 1-D ({spread:.1e})"
            assert dev_d < TOL_DISCRETE, (
                f"{name}: off the discrete upwind solution by {dev_d:.2e}")
            assert rmse_c < TOL_CONT, (
                f"{name}: off the continuum solution, RMSE {rmse_c:.2e}")

        # Budget (last step): closure, and p=0 exports exactly what dissolves.
        diss, prec, expo = rec[-1]
        close = np.abs(diss - prec - expo)[:2] / diss[:2]
        # Exact global partition (continuum): precipitated fraction of species
        # pR = 1 - export/dissolved, export = Q(W) c(W) = R W c(W) per width.
        cW = _continuum(np.array([W]), R)[0]
        frac_exp_cont = R * W * cW / (D0 * LS)
        frac_exp_exact = R * W * _discrete(R)[1][-1] / (D0 * LS)
        frac_exp = expo[0] / diss[0]
        msgs.append(
            f"  budget closure {close}; pR exported fraction {frac_exp:.5f} "
            f"(discrete {frac_exp_exact:.5f}, continuum {frac_exp_cont:.5f})"
        )
        if _rank0():
            print("\n[geochem-strip] R=%.2f m/yr  Ls=%.0f m  W=%.0f m\n"
                  % (R, LS, W) + "\n".join(msgs))
            _figure(curves, Path(request.config.rootpath) / "results"
                    / "geochem_balance")
        assert np.all(close < TOL_CLOSE)
        assert prec[1] == 0.0 and abs(expo[1] - diss[1]) / diss[1] < TOL_CLOSE
        assert abs(frac_exp - frac_exp_exact) / frac_exp_exact < TOL_DISCRETE
        assert abs(frac_exp - frac_exp_cont) / frac_exp_cont < 0.1
    finally:
        m.destroy()
