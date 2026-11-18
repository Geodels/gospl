"""
Analytical benchmark: steady-state ridge under NON-LINEAR hillslope diffusion.

Validates the three non-linear creep laws of `hillslope._hillSlopeNL` (Fortran
`hillslp_nl`) against their exact 1-D steady states.

A strip of length L, uplifted at a uniform rate U, is held at z=0 on its east
and west edges (fixed, 'f') and closed to the north and south (walls, 'w'). At
steady state the flux at distance x' from the divide carries all the uplift
upslope of it, q(x') = U x', so the slope S solves D S F(S) = U x', with F the
law's flux factor q/(D S), and z(x') is the integral of S from x' to L/2:

    linear              F = 1                    z = U/(2D) (L^2/4 - x'^2)
    Taylor series (N=2) F = 1 + (S/Sc)^2          Ganti et al. (2012); Barnhart
                                                  et al. (2019), Eq. 14
    Andrews-Bucknam     F = 1 / (1 - (S/Sc)^2)   Roering et al. (1999)
    power law (p=2)     F = S^(p-1)              Wang et al. (2024)

The non-linear profiles are integrated numerically from the exact S(x').

The mesh is a hexagonal lattice: its Voronoi cells are regular (the FV stencil
is consistent) and two of its three edge directions are at 60 degrees to the
ridge gradient. The linear law, which does not use the non-linear kernel, is run
on the same mesh as the discretisation reference: each non-linear law must match
its own exact profile as well as the linear law matches its own.

Two faults this guards against (both fixed 2026-11):
- the kernel used the slope DERIVATIVE dq/dS as the flux factor ((2j+1) series
  coefficients, p*S^(p-1), (1+f^2)/(1-f^2)^2), so no option solved its law;
- the factor was evaluated on the edge slope |dz/L|, which is S*cos(theta) on an
  edge at angle theta to the gradient, so the ridge stood too high on any mesh
  whose edges are not all aligned with the gradient, at any resolution.
"""
import numpy as np
import pytest

scipy = pytest.importorskip("scipy")
from scipy.spatial import Delaunay  # noqa: E402
from scipy.integrate import cumulative_trapezoid  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

from mpi4py import MPI  # noqa: E402
from gospl.model import Model  # noqa: E402

L = 1000.0       # strip length (m), divide at x = L/2
DX = 12.5        # lattice spacing (m)
NROW = 9         # lattice rows across the strip
U = 5.0e-4       # uplift (m/yr)
D = 0.5          # diffusivity (m^2/yr)
SC = 0.6         # critical slope: S/Sc reaches ~0.6, a strong non-linearity
TEND = 3.0e6     # >> L^2/(pi^2 D) ~ 2e5 yr: steady state
DT = 5.0e3       # implicit: the steady state does not depend on dt
# A non-linear law may be at most this much LESS accurate than the linear law
# on the same mesh, in % of its analytic divide (linear itself: about +2.4 %).
TOL_DIVIDE_PCT = 1.0
TOL_RMSE_PCT = 1.0

# name -> (extra YAML lines in the diffusion block, flux factor F(S), upper
# bracket for S: Andrews-Bucknam diverges at Sc, the others are unbounded)
LAWS = {
    "linear": ("", lambda S: 1.0, 50.0),
    "series": (f"    hillslopeSc: {SC}\n    hillslopeNb: 2\n",
               lambda S: 1.0 + (S / SC) ** 2, 50.0),
    "andrews_bucknam": (f"    hillslopeSc: {SC}\n    hillslopeNb: 0\n",
                        lambda S: 1.0 / (1.0 - (S / SC) ** 2), 0.999 * SC),
    "power2": ("    hillslopenl: 2.0\n", lambda S: S, 50.0),
}


def _analytic(F, smax, xp):
    """Steady elevation at distance xp from the divide for flux factor F:
    solve D S F(S) = U x' for S, then integrate S from x' to the boundary."""
    s = np.linspace(0.0, L / 2, 4001)
    S = np.array([0.0 if q == 0 else
                  brentq(lambda v: D * v * F(v) - q, 0.0, smax)
                  for q in U * s])
    zcum = cumulative_trapezoid(S, s, initial=0.0)
    return np.interp(xp, s, zcum[-1] - zcum)


def _write_mesh(path):
    """Hexagonal lattice strip; odd rows are offset by DX/2 and get an extra
    node on each fixed edge so both edges are straight lines x=0 and x=L."""
    dy = DX * np.sqrt(3.0) / 2.0
    pts = []
    for j in range(NROW):
        if j % 2 == 0:
            xs = np.arange(0.0, L + 0.5 * DX, DX)
        else:
            xs = np.concatenate([[0.0], np.arange(0.5 * DX, L, DX), [L]])
        pts += [(x, j * dy) for x in xs]
    p = np.array(pts)
    x, y = p[:, 0], p[:, 1]
    t = np.where((x == 0.0) | (x == L), 0.0, U)
    v = np.column_stack([x, y, np.zeros_like(x)])
    cells = Delaunay(p).simplices.astype(np.int64)
    np.savez(path, v=v, c=cells, z=np.zeros_like(x), t=t)


def _write_yaml(path, extra):
    with open(path, "w") as f:
        f.write(
            "name: non-linear hillslope benchmark\n"
            "domain:\n"
            "    npdata: ['ridge_mesh','v','c','z']\n"
            "    flowdir: 1\n"
            "    bc: 'wfwf'\n"
            "    seadepo: False\n"
            "    nodep: True\n"
            "time:\n"
            f"    start: 0.\n    end: {TEND}\n    tout: {TEND / 10}\n    dt: {DT}\n"
            "spl:\n    K: 0.\n    d: 0.\n    m: 0.5\n"
            f"diffusion:\n    hillslopeKa: {D}\n    hillslopeKm: 0.\n{extra}"
            "sea:\n    position: -1.\n"
            "climate:\n  - start: 0.\n    uniform: 0.\n"
            "tectonics:\n  - start: 0.\n    end: 1.e8\n"
            "    upsub: ['ridge_mesh','t']\n"
            "output:\n    dir: 'ridge_out'\n"
        )


def _run(workdir, name):
    """Run one law to steady state; return (divide error %, RMSE %, divide m)
    over the OWNED nodes of every rank (rank-count independent)."""
    import os

    extra, F, smax = LAWS[name]
    workdir.mkdir(parents=True, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(workdir)  # the YAML paths are relative
    _write_mesh(workdir / "ridge_mesh.npz")
    _write_yaml(workdir / "ridge.yml", extra)
    m = Model("ridge.yml", verbose=False, showlog=False)
    try:
        nonlinear = not (m.K_nl == 1.0 and m.K_sc == 0.0)
        assert nonlinear == (name != "linear"), "YAML did not select the law"
        m.runProcesses()
        own = m.inIDs == 1
        x = m.lcoords[own, 0].copy()
        z = m.hLocal.getArray()[own].copy()
    finally:
        m.destroy()
        os.chdir(cwd)
    comm = MPI.COMM_WORLD
    x = np.concatenate(comm.allgather(x))
    z = np.concatenate(comm.allgather(z))
    xp = np.abs(x - L / 2)
    za = _analytic(F, smax, xp)
    zdiv = _analytic(F, smax, np.array([0.0]))[0]
    top = xp < 0.75 * DX
    err_div = 100.0 * (z[top] - za[top]).mean() / zdiv
    rmse = 100.0 * np.sqrt(np.mean((z - za) ** 2)) / zdiv
    return err_div, rmse, zdiv


@pytest.fixture(scope="module")
def linear_reference(tmp_path_factory):
    """The linear law on the same mesh: the discretisation error to match."""
    return _run(tmp_path_factory.mktemp("hill_nl") / "linear", "linear")


@pytest.mark.benchmark
@pytest.mark.slow
@pytest.mark.parametrize("name", ["series", "andrews_bucknam", "power2"])
def test_nonlinear_hillslope_ridge(tmp_path, linear_reference, name):
    """Steady ridge of each non-linear law vs its exact profile, as accurate
    as the linear law is against its own on the same mesh."""
    lin_div, lin_rmse, _ = linear_reference
    err_div, rmse, zdiv = _run(tmp_path / name, name)
    print(f"[{name}] analytic divide {zdiv:.2f} m: divide error {err_div:+.2f} % "
          f"(linear {lin_div:+.2f} %), RMSE {rmse:.2f} % (linear {lin_rmse:.2f} %)")
    assert abs(err_div) < abs(lin_div) + TOL_DIVIDE_PCT, (
        f"{name}: divide error {err_div:+.2f} % vs linear {lin_div:+.2f} % "
        f"(tol {TOL_DIVIDE_PCT} %)")
    assert rmse - lin_rmse < TOL_RMSE_PCT, (
        f"{name}: RMSE {rmse:.2f} % vs linear {lin_rmse:.2f} % (tol {TOL_RMSE_PCT} %)")
