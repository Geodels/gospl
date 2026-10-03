"""
Analytical benchmark: marine sediment diffusion as the 2-D heat kernel.

Validates goSPL's marine deposit diffusion (`gospl.sed.hillslope._diffuseOcean`,
which calls `_diffuseImplicit` (adaptive PETSc TS, `marineSolver: ts`) or
`_diffuseImplicitPicard` (lagged-diffusivity backward Euler,
`marineSolver: picard`)) against the exact solution of the linear diffusion
equation for an instantaneous point source spread to a Gaussian:

    g(r, t) = M / (4 pi D (t + t0)) * exp(-r^2 / (4 D (t + t0)))

with `M` the deposit volume, `D` the diffusivity and `t0` the virtual age of
the initial Gaussian (initial variance per axis `2 D t0`).

Effective diffusivity in the code. Both solvers advance
`dh/dt = div(Cd grad h)` on the ABSOLUTE surface `h = bed + deposit`, with

    Cd = minDiff + nlK * (1 - exp(-dexp * deposit)),   deposit >= 0.1 m
    Cd = minDiff,                                       deposit <  0.1 m

(`hillslope._evalFunctionMardDiff` / `_diffuseImplicitPicard`; `minDiff = 1e-4`,
`dexp = 0.05` are hard-coded in `hillSLP.__init__`, `nlK` = YAML
`diffusion: nonlinKm`). Cd is therefore a function of the local deposit
THICKNESS, not of water depth. To make it a constant we put the Gaussian on a
uniform background deposit `B = 300 m` (so `exp(-dexp*deposit) <= 3e-7`): then
`D = nlK + minDiff` to ~1e-6 relative everywhere, the uniform part carries no
flux, the bed is flat (no bed-driven flux), the `ndepo >= 0` clamp never fires
and the `_diffuseOcean` volume rescale is the identity. What remains is linear
diffusion of the Gaussian, which is compared with the heat kernel above.

The whole mesh is marine (flat bed at -500 m, sea level 0), so the marine mask
covers the domain and the only boundary is the mesh edge, which the Gaussian
never reaches (domain half-width ~4.5 final standard deviations).

Note on the TS tolerances: `_diffuseImplicit` uses `atol=5e-3, rtol=1e-4` on
the ABSOLUTE surface, so its local error budget scales with |bed + deposit|
(here ~200 m -> ~2 cm per step on a 20 m Gaussian). The bed is placed at -500 m
rather than in the abyss for that reason; the physics is depth-independent.

Checks, for both solvers, after N_STEPS model steps of `dt`:
  - volume of the Gaussian conserved (relative < 1e-8);
  - peak decay: peak(t) * (t + t0) constant (the 1/(t+t0) law);
  - effective D from the growth of the second moment <r^2> = 4 D (t+t0);
  - RMSE of the profile vs the exact kernel, relative to the exact peak.
"""
import time
from pathlib import Path

import numpy as np
import pytest

scipy = pytest.importorskip("scipy")
from scipy.spatial import Delaunay  # noqa: E402

from mpi4py import MPI  # noqa: E402
from petsc4py import PETSc  # noqa: E402
from gospl.model import Model  # noqa: E402

# ---- Geometry / parameters (SI, years) --------------------------------------
NX = NY = 101            # 101 x 101 structured grid
DX = 1000.0              # 1 km spacing -> 100 km square
L = (NX - 1) * DX
XC = YC = 0.5 * L        # Gaussian centre
BED = -500.0             # flat seafloor (m), sea level 0
BACKGROUND = 300.0       # uniform background deposit (m) -> Cd ~ constant
AMP0 = 20.0              # initial Gaussian peak (m)
NLK = 1.0e4              # diffusion: nonlinKm (m^2/yr)
MINDIFF = 1.0e-4         # hard-coded hillSLP.minDiff
D_EFF = NLK + MINDIFF    # effective constant diffusivity
SIGMA0 = 5000.0          # initial std dev per axis (m)
T0 = SIGMA0 ** 2 / (2.0 * D_EFF)   # virtual age of the initial Gaussian (yr)
DT = 1000.0              # model step (yr)
N_STEPS = 5              # -> sigma grows 5 km -> ~11.2 km
PICARD_SUB = 50          # backward-Euler sub-steps per dt (Picard)

# Tolerances. Measured (this configuration; identical at np=1, 2 and 4):
#   picard: D err 1.3e-4, peak-law dev 6.4e-3, RMSE/peak 5.1e-4
#   ts:     D err 3.3e-5, peak-law dev 1.2e-3, RMSE/peak 1.5e-4
# History: until 2026-10 the cached marine TS (`_diffuseImplicit`) used
# `ksp_type preonly` + `pc gasm`, i.e. ONE preconditioner application per
# Rosenbrock stage instead of a stage solve. It under-diffused (D err 4.9e-2
# here; D_eff/D = 0.973 / 0.868 for single steps of 1e3 / 5e3 yr) and depended
# on the partition (D err 0.26 at np=2). The stage system is now Krylov-solved
# (gmres + gasm, rtol 1e-8), and both solvers share the same tolerances.
TOL_VOLUME = 1.0e-8
TOL_PEAK_LAW = {"ts": 1.0e-2, "picard": 1.0e-2}  # max dev of peak*(t+t0)
TOL_D = {"ts": 2.0e-3, "picard": 2.0e-3}         # moment-derived D, relative
TOL_RMSE = {"ts": 2.0e-3, "picard": 2.0e-3}      # RMSE / exact peak


def _write_mesh(path, bed=None):
    xs = np.arange(NX) * DX
    ys = np.arange(NY) * DX
    X, Y = np.meshgrid(xs, ys)
    x, y = X.ravel(), Y.ravel()
    v = np.column_stack([x, y, np.zeros_like(x)])
    cells = Delaunay(np.column_stack([x, y])).simplices.astype(np.int64)
    z = np.full(x.shape, BED if bed is None else bed, dtype=np.float64)
    np.savez(path, v=v, c=cells, z=z)


def _write_yaml(path, solver):
    with open(path, "w") as f:
        f.write(
            "name: marine diffusion heat-kernel benchmark\n"
            "domain:\n"
            "    npdata: ['marine_mesh','v','c','z']\n"
            "    flowdir: 2\n"
            "    bc: 'oooo'\n"
            "    seadepo: True\n"
            "time:\n"
            f"    start: 0.\n    end: {DT * N_STEPS}\n"
            f"    tout: {DT * N_STEPS}\n    dt: {DT}\n"
            "spl:\n    K: 1.0e-20\n    d: 0.\n    m: 0.5\n"
            "diffusion:\n"
            "    hillslopeKa: 0.\n    hillslopeKm: 0.\n"
            f"    nonlinKm: {NLK}\n"
            f"    marineSolver: '{solver}'\n"
            f"    picardSub: {PICARD_SUB}\n"
            "    picardIts: 2\n"
            "sea:\n    position: 0.\n"
            "climate:\n  - start: 0.\n    uniform: 1.\n"
            f"output:\n    dir: 'marine_out_{solver}'\n"
        )


def _gaussian(r2, t):
    """Exact heat kernel thickness at squared radius r2, elapsed time t."""
    s = 4.0 * D_EFF * (t + T0)
    mass = AMP0 * 4.0 * np.pi * D_EFF * T0
    return mass / (np.pi * s) * np.exp(-r2 / s)


def _gsum(arr):
    return MPI.COMM_WORLD.allreduce(float(arr), op=MPI.SUM)


def _gmax(arr):
    return MPI.COMM_WORLD.allreduce(float(arr), op=MPI.MAX)


def _run_solver(tmp_path, solver):
    """Diffuse the Gaussian N_STEPS steps; return per-step diagnostics."""
    _write_yaml(tmp_path / f"marine_{solver}.yml", solver)
    m = Model(f"marine_{solver}.yml", verbose=False, showlog=False)
    try:
        assert m.flatModel
        assert m.marineSolver == solver
        bed = m.hLocal.getArray().copy()
        m.sealevel = 0.0
        m.seaID = np.where(bed <= m.sealevel)[0]
        owned = m.inIDs == 1
        A = m.larea
        x, y = m.lcoords[:, 0], m.lcoords[:, 1]
        r2 = (x - XC) ** 2 + (y - YC) ** 2

        g = _gaussian(r2, 0.0)
        mass0 = _gsum(np.sum((g * A)[owned]))
        out = {"t": [], "peak": [], "mass": [], "r2m": [], "rmse": [],
               "maxerr": [], "profile": None}
        t = 0.0
        wall = time.perf_counter()
        for _ in range(N_STEPS):
            m._diffuseOcean(BACKGROUND + g)
            m.dm.globalToLocal(m.tmp, m.tmpL)
            g = m.tmpL.getArray().copy() - BACKGROUND
            t += m.dt
            ex = _gaussian(r2, t)
            err = (g - ex)[owned]
            n = MPI.COMM_WORLD.allreduce(int(owned.sum()), op=MPI.SUM)
            mass = _gsum(np.sum((g * A)[owned]))
            out["t"].append(t)
            out["peak"].append(_gmax(g[owned].max()))
            out["mass"].append(mass)
            out["r2m"].append(_gsum(np.sum((g * A * r2)[owned])) / mass)
            out["rmse"].append(np.sqrt(_gsum(np.sum(err ** 2)) / n)
                               / ex[owned].max())
            out["maxerr"].append(_gmax(np.abs(err).max()) / ex[owned].max())
        out["wall"] = time.perf_counter() - wall
        out["mass0"] = mass0
        out["final"] = g.copy()
        out["owned"] = owned
        out["ts_steps"] = (m._ts_marine.getStepNumber()
                           if getattr(m, "_ts_marine", None) is not None else 0)
        # Centre-line profile (y == YC) for the figure, serial only.
        line = owned & (np.abs(y - YC) < 0.5 * DX)
        order = np.argsort(x[line])
        out["profile"] = (x[line][order] - XC, g[line][order],
                          _gaussian(r2, t)[line][order],
                          _gaussian(r2, 0.0)[line][order])
        return out
    finally:
        m.destroy()


def _figure(results, outdir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for solver, res in results.items():
        xr, g, ex, g0 = res["profile"]
        ax[0].plot(xr / 1e3, g, lw=1.5, label=f"goSPL ({solver})")
        ax[1].plot(res["t"], np.array(res["rmse"]) * 100, "o-", label=solver)
    ax[0].plot(xr / 1e3, ex, "k--", lw=1, label="exact")
    ax[0].plot(xr / 1e3, g0, color="0.6", lw=1, label="initial")
    ax[0].set_xlabel("distance from centre (km)")
    ax[0].set_ylabel("deposit above background (m)")
    ax[0].legend(frameon=False)
    ax[1].set_xlabel("elapsed time (yr)")
    ax[1].set_ylabel("RMSE / exact peak (%)")
    ax[1].legend(frameon=False)
    fig.tight_layout()
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / "marine_diffusion.png", dpi=120)
    plt.close(fig)


@pytest.mark.benchmark
@pytest.mark.slow
def test_marine_diffusion_heat_kernel(tmp_path, monkeypatch, request):
    """Marine deposit diffusion vs the 2-D heat kernel (TS and Picard)."""
    monkeypatch.chdir(tmp_path)
    _write_mesh(tmp_path / "marine_mesh.npz")

    results = {s: _run_solver(tmp_path, s) for s in ("ts", "picard")}
    rank0 = PETSc.COMM_WORLD.getRank() == 0
    if rank0:
        _figure(results,
                Path(request.config.rootpath) / "results" / "marine_diffusion")

    for solver, res in results.items():
        t = np.array(res["t"])
        peak = np.array(res["peak"])
        mass = np.array(res["mass"])
        r2m = np.array(res["r2m"])
        dvol = float(np.max(np.abs(mass - res["mass0"]) / res["mass0"]))
        law = peak * (t + T0)
        law0 = AMP0 * T0
        dlaw = float(np.max(np.abs(law / law0 - 1.0)))
        # <r^2> = 4 D (t + t0): slope over the run gives D (initial moment
        # taken from the discrete initial field via the exact 4 D t0).
        d_fit = float(np.polyfit(t, r2m, 1)[0] / 4.0)
        d_err = abs(d_fit - D_EFF) / D_EFF
        rmse = float(res["rmse"][-1])
        if rank0:
            print(
                f"\n[marine-diffusion:{solver}] D={D_EFF:.4g} t0={T0:.0f} yr  "
                f"vol err={dvol:.2e}  peak*(t+t0) dev={dlaw:.2e}  "
                f"D_fit={d_fit:.5g} (err {d_err:.2e})  "
                f"RMSE/peak={rmse:.2e}  max|err|/peak={res['maxerr'][-1]:.2e}  "
                f"wall={res['wall']:.1f} s"
            )
        assert dvol < TOL_VOLUME, f"{solver}: volume not conserved ({dvol:.2e})"
        assert dlaw < TOL_PEAK_LAW[solver], (
            f"{solver}: peak does not decay as 1/(t+t0) (max dev {dlaw:.2e})"
        )
        assert d_err < TOL_D[solver], (
            f"{solver}: moment-derived D={d_fit:.5g} vs {D_EFF:.5g} "
            f"(rel err {d_err:.2e})"
        )
        assert rmse < TOL_RMSE[solver], (
            f"{solver}: profile RMSE/peak {rmse:.2e} > {TOL_RMSE[solver]:.0e}"
        )


@pytest.mark.benchmark
@pytest.mark.slow
def test_marine_diffusion_depth_invariance(tmp_path, monkeypatch):
    """The same deposit on a 500 m and a 5000 m deep seafloor must diffuse
    identically: the TS error control acts on the deposit, not on the absolute
    surface `bed + deposit` (until 2026-10 a scalar rtol on |bed + deposit|
    loosened the tolerance tenfold at 5000 m: 1.6e-4 of the peak apart)."""
    finals = {}
    for bed in (-500.0, -5000.0):
        sub = tmp_path / f"bed{int(-bed)}"
        sub.mkdir()
        monkeypatch.chdir(sub)
        _write_mesh(sub / "marine_mesh.npz", bed=bed)
        res = _run_solver(sub, "ts")
        finals[bed] = (res["final"], res["owned"])
    (a, own), (b, _) = finals[-500.0], finals[-5000.0]
    dev = _gmax(np.abs(a - b)[own].max()) / _gmax(a[own].max())
    if PETSc.COMM_WORLD.getRank() == 0:
        print(f"\n[marine-diffusion depth invariance] max|d|/peak = {dev:.2e}")
    assert dev < 1.0e-9, f"deposit depends on water depth ({dev:.2e} of peak)"
