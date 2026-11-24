"""
Analytical benchmarks: benthic carbonate production (``carbonate:``,
docs/DESIGN_CARBONATES.md, phase C1).

One factory, all of it framework (built in place), on a seafloor whose water
depth ``w`` is the same everywhere. Production follows Bosscher & Schlager
(1992), ``P(w) = g tanh(I0 exp(-k w) / Ik)``, and the water depth evolves as

    dw/dt = sigma - P(w)

under a relative sea-level rise ``sigma``. Three exact properties:

- keep-up or give-up: ``w*`` with ``P(w*) = sigma``, i.e.
  ``w* = ln(I0 / (Ik atanh(sigma/g))) / k``, is an unstable equilibrium
  (production grows as the water shallows): a platform shallower than ``w*``
  reaches sea level and keeps up there, a deeper one drowns;
- drowning: for ``sigma > g`` the depth grows without bound, at a rate that
  tends to ``sigma``;
- the integration over a step does not depend on the time step: with constant
  sea level the depth after a time ``T`` matches the ODE solution (solved here
  with a tight tolerance) for ``dt = 1 kyr`` and ``dt = 10 kyr``.

``carbProduction`` is driven directly (sea level and the filled surface set by
the test); the river and marine processes are not run.
"""
import os
import re
import shutil

import numpy as np
import pytest

pytest.importorskip("scipy")
from scipy.integrate import solve_ivp  # noqa: E402

from gospl.model import Model  # noqa: E402

FIX = os.path.join(os.path.dirname(__file__), "..", "tests", "fixtures")
G, K, IK, I0 = 500.0, 0.1, 60.0, 400.0          # m/Myr, 1/m, W/m2, W/m2


def _rate(w):
    return G * 1e-6 * np.tanh(I0 * np.exp(-K * np.maximum(w, 0.0)) / IK)   # m/yr


def _model(tmp_path, dt):
    src = os.path.join(FIX, "minimal_dual.yml")
    if not os.path.exists(src):
        pytest.skip("minimal_dual.yml not present")
    shutil.copy(os.path.join(FIX, "mesh.npz"), tmp_path / "mesh.npz")
    text = re.sub(r"dir: '[^']*'", "dir: 'carb_out'", open(src).read())
    text = text.replace("output:", (
        "carbonate:\n    insolation: %r\n    factories:\n"
        "      - {name: reef, growth: %r, extinction: %r, saturation: %r, "
        "framework: 1.0, sand: 1.0}\n\noutput:" % (I0, G, K, IK)), 1)
    (tmp_path / "carb.yml").write_text(text)
    m = Model("carb.yml", verbose=False, showlog=False)
    m.dt = dt
    return m


def _evolve(m, w0, sigma, tend):
    """Water depth after ``tend`` years from ``w0`` (uniform), relative sea
    level rising at ``sigma`` (m/yr); returns the depths after every step."""
    h = m.hLocal.getArray()
    h[:] = -w0
    m.hLocal.setArray(h)
    m.dm.localToGlobal(m.hLocal, m.hGlobal)
    m.sealevel = 0.0
    out = []
    t = 0.0
    while t < tend - 1e-9:
        m.sealevel += sigma * m.dt
        m.lFill = m.hLocal.getArray().copy()
        m.carbProduction()
        t += m.dt
        out.append(m.sealevel - m.hLocal.getArray().copy())
    return np.array(out)


@pytest.mark.benchmark
@pytest.mark.slow
@pytest.mark.parametrize("side", [-1.0, 1.0])
def test_carbonate_keep_up_threshold(tmp_path, monkeypatch, side):
    """Production grows as the water shallows, so w* (P(w*) = sigma) is the
    give-up threshold: half a metre shallower the platform outgrows the rise,
    reaches sea level and keeps up there (growth limited to sigma); half a
    metre deeper it falls behind and drowns."""
    sigma = 200.0e-6
    wstar = np.log(I0 / (IK * np.arctanh(sigma / (G * 1e-6)))) / K
    monkeypatch.chdir(tmp_path)
    m = _model(tmp_path, 1000.0)
    try:
        w = _evolve(m, wstar + 0.5 * side, sigma, 5.0e5)
        print("\n[carbonate threshold] w* = %.3f m, start %+.1f m: depth %.3f m"
              % (wstar, 0.5 * side, w[-1].mean()))
        if side < 0:
            assert np.all(w[-1] < 1e-6), "platform above w* did not keep up"
            # keeping up: the depth stays at sea level over the last steps
            assert np.max(np.abs(w[-50:])) < 1e-6
        else:
            assert np.all(w[-1] > wstar + 10.0), "platform below w* did not drown"
    finally:
        m.destroy()


@pytest.mark.benchmark
@pytest.mark.slow
def test_carbonate_drowning(tmp_path, monkeypatch):
    """A relative sea-level rise faster than the maximum growth rate drowns
    the platform: the depth keeps increasing at a rate tending to sigma."""
    sigma = 600.0e-6
    monkeypatch.chdir(tmp_path)
    m = _model(tmp_path, 1000.0)
    try:
        w = _evolve(m, 5.0, sigma, 4.0e5)
        rate = (w[-1] - w[-11]).mean() / (10 * m.dt)
        print("\n[carbonate drowning] depth %.1f m, deepening %.3g m/yr (sigma %.3g)"
              % (w[-1].mean(), rate, sigma))
        assert np.all(np.diff(w.mean(axis=1)) > 0.0)
        assert abs(rate / sigma - 1.0) < 0.02
    finally:
        m.destroy()


@pytest.mark.benchmark
@pytest.mark.slow
@pytest.mark.parametrize("dt", [1000.0, 10000.0])
def test_carbonate_step_independent(tmp_path, monkeypatch, dt):
    """Constant sea level: the depth after 100 kyr matches the ODE solution
    dw/dt = -P(w) for 1 and 10 kyr steps."""
    w0, tend = 40.0, 1.0e5
    ref = solve_ivp(lambda t, w: -_rate(w), (0.0, tend), [w0], rtol=1e-10,
                    atol=1e-10).y[0, -1]
    monkeypatch.chdir(tmp_path)
    m = _model(tmp_path, dt)
    try:
        w = _evolve(m, w0, 0.0, tend)
        print("\n[carbonate dt=%g] depth %.4f m, ODE %.4f m" % (dt, w[-1].mean(), ref))
        assert np.allclose(w[-1], ref, rtol=1e-4)
    finally:
        m.destroy()
