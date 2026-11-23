"""
Analytical benchmark: gravity-current deposit over a run-out length.

Validates the routed deposition of the marine gravity currents
(`gospl.sed.seaplex._gravRoute`, docs/DESIGN_GRAVITY_CURRENTS.md 2.1-2.2)
against its exact solution. On a uniform seafloor slope, a supply Q0 entering a
cross-slope line of cells loses a fraction Lambda = 1 - exp(-l/L) of its flux in
every cell it crosses (l the cell size, L the run-out length), so the deposit
summed across the slope in the k-th column below the source is

    D_k = Q0 (1 - e^(-dx/L)) e^(-k dx/L)

i.e. it decays as e^(-x/L) whatever the mesh spacing.

Setup: a square lattice (cells of exactly dx x dx, so l = dx), a plane dipping
east (all of it deeper than twice the bypass depth, so the depth taper is 1),
walls north, south and west, an open east edge that exports what is left. On a
square lattice every receiver (straight or diagonal) is exactly one column
further down-slope, so the cross-slope sums follow D_k exactly for any number of
receivers. The supply is injected in a band of the first interior column,
far enough from the side walls that the flow never reaches them (the wall
cells are half cells, so their size l = sqrt(area) is dx / sqrt(2): they
deposit less, and a narrow strip made the decay length depend on dx). The
triangulation's diagonals all point north-east, so with several receivers the
flux drifts north (about a third of a cell per column): the band starts in the
south so that it stays clear of the north wall. The
operator is driven directly (`_gravNetwork`, `_gravRoute`), with no deposit
spreading. Two resolutions with the same L check that the decay length does not
depend on the spacing.
"""
import numpy as np
import pytest

pytest.importorskip("scipy")
from mpi4py import MPI  # noqa: E402
from petsc4py import PETSc  # noqa: E402

from gospl.model import Model  # noqa: E402

LX, LY = 60.0e3, 40.0e3    # domain (m): the deposit decays over x
Y_SRC, BAND = 8.0e3, 2.0e3  # supply in |y - Y_SRC| <= BAND, far from the walls
RUNOUT = 1.0e4             # run-out length L (m)
SLOPE = 0.01               # seafloor gradient, dipping east
Z0 = -300.0                # depth at x = 0 (m), deeper than 2 x bypass depth
Q0 = 1.0e9                 # supply (m3 over the step)
TOL = 1.0e-6               # relative error on every column


def _write_case(path, dx, receivers):
    nx, ny = int(round(LX / dx)) + 1, int(round(LY / dx)) + 1
    X, Y = np.meshgrid(np.arange(nx) * dx, np.arange(ny) * dx)
    x, y = X.ravel(), Y.ravel()
    cells = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            a = j * nx + i
            cells += [[a, a + 1, a + nx + 1], [a, a + nx + 1, a + nx]]
    v = np.column_stack([x, y, np.zeros_like(x)])
    z = Z0 - SLOPE * x
    np.savez(path / "grav_mesh.npz", v=v, c=np.array(cells, dtype=np.int64), z=z)
    (path / "grav.yml").write_text(
        "name: gravity run-out benchmark\n"
        "domain:\n"
        "    npdata: ['grav_mesh','v','c','z']\n"
        "    flowdir: 1\n"
        "    bc: 'woww'\n"              # N wall, E open (export), S and W walls
        "    seadepo: True\n"
        "time:\n    start: 0.\n    end: 10.\n    tout: 10.\n    dt: 10.\n"
        "spl:\n    K: 1.0e-20\n    d: 0.\n    m: 0.5\n"
        "diffusion:\n    hillslopeKa: 0.\n    hillslopeKm: 0.\n"
        "sea:\n    position: 0.\n"
        "climate:\n  - start: 0.\n    uniform: 1.\n"
        "gravity:\n"
        f"    runout: {RUNOUT}\n"
        "    fraction: 0.5\n"
        f"    receivers: {receivers}\n"
        "    deposit_spread: 0.\n"
        "output:\n    dir: 'grav_out'\n"
    )
    return nx


def _columns(m, vals, dx, nx):
    """Sum a local array over each x column (owned nodes, all ranks)."""
    own = m.inIDs == 1
    k = np.rint(m.lcoords[own, 0] / dx).astype(int)
    out = np.bincount(k, weights=vals[own], minlength=nx)
    return MPI.COMM_WORLD.allreduce(out, op=MPI.SUM)


@pytest.mark.benchmark
@pytest.mark.slow
@pytest.mark.parametrize("dx,receivers", [(1000.0, 1), (500.0, 1), (500.0, 6)])
def test_gravity_runout_exponential(tmp_path, monkeypatch, dx, receivers):
    """Deposit on a uniform slope = Q0 (1 - e^(-dx/L)) e^(-k dx/L), exactly,
    at two spacings and with one or six receivers."""
    monkeypatch.chdir(tmp_path)
    nx = _write_case(tmp_path, dx, receivers)
    m = Model("grav.yml", verbose=False, showlog=False)
    try:
        assert m.gravOn and m.flatModel
        hl = m.hLocal.getArray().copy()
        m.sinkIDs = np.ones(m.lpoints, dtype=bool)
        m._gravNetwork(hl)
        x = m.lcoords[:, 0]
        y = m.lcoords[:, 1]
        src = (np.rint(x / dx) == 1) & (np.abs(y - Y_SRC) <= BAND) & (m.inIDs == 1)
        nsrc = MPI.COMM_WORLD.allreduce(int(src.sum()), op=MPI.SUM)
        supply = np.where(src, Q0 / nsrc, 0.0)
        dep, exp = m._gravRoute(supply, m.gravRunout)
        m._gravRelease()

        D = _columns(m, dep, dx, nx)
        E = _columns(m, exp, dx, nx)
        k = np.arange(nx) - 1                     # columns below the source
        lam = 1.0 - np.exp(-dx / RUNOUT)
        ref = np.where(k >= 0, Q0 * lam * np.exp(-np.maximum(k, 0) * dx / RUNOUT), 0.0)
        ref[-1] = 0.0                             # the open edge exports instead
        inner = slice(1, nx - 1)
        err = np.abs(D[inner] - ref[inner]) / Q0
        budget = abs(D.sum() + E.sum() - Q0) / Q0
        # decay length fitted on the deposit (log-linear), against L
        kk = np.arange(1, nx - 1)
        Lfit = -dx / np.polyfit(kk, np.log(D[inner]), 1)[0]
        if PETSc.COMM_WORLD.getRank() == 0:
            print(f"\n[gravity runout] dx={dx:.0f} m receivers={receivers}: "
                  f"max column error {err.max():.2e} Q0, budget {budget:.1e}, "
                  f"fitted L {Lfit:.1f} m (L = {RUNOUT:.0f} m)")
        assert budget < 1.0e-8, f"supply not conserved ({budget:.1e})"
        assert err.max() < TOL, f"deposit off the exponential by {err.max():.2e} Q0"
        assert abs(Lfit / RUNOUT - 1.0) < 1.0e-6
    finally:
        m.destroy()
