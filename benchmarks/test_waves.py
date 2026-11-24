"""
Analytical benchmarks: wave energy and shadows (``waves:``,
docs/DESIGN_CARBONATES.md 2.10, phase C7).

On a square lattice the first-order upwind operator along a mesh axis is
exact in 1-D: with a constant attenuation kappa, (d . grad + kappa) E = 0 and
E = 1 on the inflow edge give E_i = (1 + kappa dx)^(-i) at i cells from the
edge, which tends to exp(-kappa x) at first order in kappa dx. An island
blocks the waves: the cells straight downwind of it get no energy (the
operator has no lateral spreading along an axis), the open water beside it is
unchanged. A climate of two directions gives the weighted sum of the two
fields. The half cells of the walls and of the outflow edge follow a
different discrete balance and are left out of the exact comparisons.
"""
import numpy as np
import pytest

pytest.importorskip("scipy")
from gospl.model import Model  # noqa: E402

LX, LY = 40.0e3, 20.0e3
DECAY = 1.0e4                 # deep-water attenuation length (m)


def _case(path, dx, climate, island=False):
    nx, ny = int(round(LX / dx)) + 1, int(round(LY / dx)) + 1
    X, Y = np.meshgrid(np.arange(nx) * dx, np.arange(ny) * dx)
    x, y = X.ravel(), Y.ravel()
    cells = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            a = j * nx + i
            cells += [[a, a + 1, a + nx + 1], [a, a + nx + 1, a + nx]]
    z = np.full_like(x, -500.0)                    # deep everywhere: kappa = 1/DECAY
    if island:
        z[(np.abs(x - 10.0e3) <= 1.0e3) & (np.abs(y - 10.0e3) <= 2.0e3)] = 50.0
    np.savez(path / "wave.npz", v=np.column_stack([x, y, np.zeros_like(x)]),
             c=np.array(cells, dtype=np.int64), z=z)
    clim = "".join("      - {from: %r, weight: %r}\n" % c for c in climate)
    (path / "wave.yml").write_text(
        "name: wave energy\n"
        "domain:\n    npdata: ['wave','v','c','z']\n    flowdir: 1\n    bc: 'wwww'\n"
        "    seadepo: True\n"
        "time:\n    start: 0.\n    end: 10.\n    tout: 10.\n    dt: 10.\n"
        "spl:\n    K: 1.0e-20\n    d: 0.\n    m: 0.5\n"
        "diffusion:\n    hillslopeKa: 0.\n    hillslopeKm: 0.\n"
        "sea:\n    position: 0.\n"
        "climate:\n  - start: 0.\n    uniform: 1.\n"
        "waves:\n    decay: %r\n    shoal_depth: 20.\n    shoal_length: 2.e3\n"
        "    climate:\n%s"
        "output:\n    dir: 'wave_out'\n" % (DECAY, clim))
    return nx


@pytest.mark.benchmark
@pytest.mark.slow
@pytest.mark.parametrize("dx", [500.0, 250.0])
def test_wave_strip_decay(tmp_path, monkeypatch, dx):
    """Waves from the west on a deep strip: E = (1 + dx/L)^(-i) exactly,
    exp(-x/L) to first order."""
    from mpi4py import MPI

    monkeypatch.chdir(tmp_path)
    _case(tmp_path, dx, [(270.0, 1.0)])
    m = Model("wave.yml", verbose=False, showlog=False)
    try:
        m.sealevel = 0.0
        E = m.waveEnergy()
        x, y = m.lcoords[:, 0], m.lcoords[:, 1]
        i = np.rint(x / dx)
        # away from the side walls and the outflow edge, whose half cells
        # follow a different discrete balance
        own = (m.inIDs == 1) & (y > 2 * dx) & (y < LY - 2 * dx) & (x < LX - 0.5 * dx)
        disc = np.max(np.abs(E[own] - (1.0 + dx / DECAY) ** (-i[own])))
        cont = np.max(np.abs(E[own] - np.exp(-x[own] / DECAY)))
        disc = MPI.COMM_WORLD.allreduce(disc, op=MPI.MAX)
        cont = MPI.COMM_WORLD.allreduce(cont, op=MPI.MAX)
        print("\n[waves] dx=%g: vs (1+dx/L)^-i %.2e, vs exp(-x/L) %.2e" % (dx, disc, cont))
        assert disc < 1.0e-8
        assert cont < 0.5 * dx / DECAY
    finally:
        m.destroy()


@pytest.mark.benchmark
@pytest.mark.slow
def test_wave_island_shadow(tmp_path, monkeypatch):
    """An island in a westerly swell: nothing behind it, the open water
    beside it as without the island."""
    from mpi4py import MPI

    dx = 250.0
    monkeypatch.chdir(tmp_path)
    _case(tmp_path, dx, [(270.0, 1.0)], island=True)
    m = Model("wave.yml", verbose=False, showlog=False)
    try:
        m.sealevel = 0.0
        E = m.waveEnergy()
        x, y = m.lcoords[:, 0], m.lcoords[:, 1]
        own = m.inIDs == 1
        lee = own & (x > 11.5e3) & (np.abs(y - 10.0e3) <= 1.5e3)
        side = own & (x > 11.5e3) & (np.abs(y - 10.0e3) >= 3.5e3) & (y > 2 * dx) \
            & (y < LY - 2 * dx) & (x < LX - 0.5 * dx)
        i = np.rint(x / dx)
        lee_max = MPI.COMM_WORLD.allreduce(float(E[lee].max()) if lee.any() else 0.0,
                                           op=MPI.MAX)
        side_err = MPI.COMM_WORLD.allreduce(
            float(np.max(np.abs(E[side] - (1.0 + dx / DECAY) ** (-i[side]))))
            if side.any() else 0.0, op=MPI.MAX)
        print("\n[waves] island: max E in the lee %.2e, open water error %.2e"
              % (lee_max, side_err))
        assert lee_max < 1.0e-10
        assert side_err < 1.0e-8
    finally:
        m.destroy()


@pytest.mark.benchmark
@pytest.mark.slow
def test_wave_climate_weights(tmp_path, monkeypatch):
    """Two directions: the field is the weighted sum of the single-direction
    fields."""
    monkeypatch.chdir(tmp_path)
    fields = {}
    for name, clim in (("w", [(270.0, 1.0)]), ("s", [(180.0, 1.0)]),
                       ("ws", [(270.0, 0.7), (180.0, 0.3)])):
        _case(tmp_path, 500.0, clim)
        m = Model("wave.yml", verbose=False, showlog=False)
        try:
            m.sealevel = 0.0
            fields[name] = m.waveEnergy().copy()
        finally:
            m.destroy()
    assert np.allclose(fields["ws"], 0.7 * fields["w"] + 0.3 * fields["s"], atol=1e-10)


def _fetch_case(path, dx, fetch, island):
    """As _case, with a fetch (growth) term."""
    _case(path, dx, [(270.0, 1.0)], island=island)
    yml = (path / "wave.yml").read_text().replace(
        "    shoal_depth: 20.\n", "    shoal_depth: 20.\n    fetch: %r\n" % fetch, 1)
    (path / "wave.yml").write_text(yml)


@pytest.mark.benchmark
@pytest.mark.slow
def test_wave_fetch_growth(tmp_path, monkeypatch):
    """With a fetch L_f the upwind scheme along the waves is the recurrence
    E_i = (E_(i-1) + g dx) / (1 + (g + k) dx), g = 1/L_f, k = 1/decay: from
    the inflow edge (E_0 = 1) towards E_inf = g / (g + k), and behind an
    island from 0 again (fetch-limited regrowth)."""
    from mpi4py import MPI

    dx, fetch = 250.0, 5.0e3
    g, k = 1.0 / fetch, 1.0 / DECAY
    monkeypatch.chdir(tmp_path)
    _fetch_case(tmp_path, dx, fetch, island=True)
    m = Model("wave.yml", verbose=False, showlog=False)
    try:
        m.sealevel = 0.0
        E = m.waveEnergy()
        x, y = m.lcoords[:, 0], m.lcoords[:, 1]
        i = np.rint(x / dx).astype(int)
        nmax = int(np.rint(LX / dx)) + 1
        # open water row: from the inflow edge
        rec = np.empty(nmax)
        rec[0] = 1.0
        for n in range(1, nmax):
            rec[n] = (rec[n - 1] + g * dx) / (1.0 + (g + k) * dx)
        own = m.inIDs == 1
        side = own & (np.abs(y - 10.0e3) >= 3.5e3) & (y > 2 * dx) & (y < LY - 2 * dx) \
            & (x < LX - 0.5 * dx)
        e1 = float(np.max(np.abs(E[side] - rec[i[side]]))) if side.any() else 0.0
        # lee of the island: from 0 at the last island column (x = 11 km)
        i0 = int(np.rint(11.0e3 / dx))
        lee_rec = np.zeros(nmax)
        for n in range(i0 + 1, nmax):
            lee_rec[n] = (lee_rec[n - 1] + g * dx) / (1.0 + (g + k) * dx)
        lee = own & (x > 11.0e3 + 0.5 * dx) & (np.abs(y - 10.0e3) <= 1.5e3) \
            & (x < LX - 0.5 * dx)
        e2 = float(np.max(np.abs(E[lee] - lee_rec[i[lee]]))) if lee.any() else 0.0
        e1 = MPI.COMM_WORLD.allreduce(e1, op=MPI.MAX)
        e2 = MPI.COMM_WORLD.allreduce(e2, op=MPI.MAX)
        print("\n[waves fetch] open water %.2e, lee regrowth %.2e (E_inf %.4f)"
              % (e1, e2, g / (g + k)))
        assert e1 < 1.0e-10 and e2 < 1.0e-10
    finally:
        m.destroy()
