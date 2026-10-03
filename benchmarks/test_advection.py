"""
Analytical benchmark: horizontal advection of a passive surface field.

Validates goSPL's finite-volume horizontal advection
(``gospl.mesher.tectonics._varAdvector``, YAML ``domain: advect:``) against the
exact solution of the linear advection equation with a uniform velocity ``v``:

    dh/dt + v . grad(h) = 0    =>    h(x, t) = h0(x - v t)

i.e. the initial profile is translated rigidly by ``v t`` with no change of
shape or volume. Three schemes are exercised:

* ``upwind`` (``advscheme=1``): first-order implicit upwind. Unconditionally
  stable but diffusive. Its modified equation adds an along-flow numerical
  diffusivity ``D_num ~ U dx / 2 + U^2 dt / 2 = (U dx / 2)(1 + C)`` (spatial
  upwinding plus backward Euler), with ``C = U dt / dx`` the Courant number. A
  Gaussian of width ``s0`` therefore widens along the flow to
  ``s^2 = s0^2 + 2 D_num t`` and, being stretched in one direction only, its
  peak decays to ``s0 / s``. The benchmark measures ``D_num`` from the second
  moment of the advected field and checks it against that estimate.
* ``iioe1`` (``advscheme=2``): Inflow-Implicit / Outflow-Explicit (Mikula and
  Ohlberger 2014), formally second order for smooth solutions: almost no
  diffusion, small dispersive undershoots allowed.
* ``iioe2`` (``advscheme=3``): IIOE plus the Scheme-2 anti-overshoot correction,
  applied only where the Scheme-1 result leaves the local neighbourhood range.
  It must stay bounded (no new extrema). The *no-overshoot* path is exercised
  separately with a uniform field: there the global excess is exactly zero and
  the correction is skipped. That path zeroed every field before the 2026-06-29
  fix (see AGENTS.md > Fixed > "IIOE2 advection").

Setup. A flat hexagonal-lattice mesh (rows offset by dx/2, so the Voronoi cells
are regular hexagons and the Delaunay triangulation is non-degenerate) built in
``tmp_path`` at runtime. A Gaussian bump ``A exp(-r^2 / (2 s0^2))`` on a zero
background is advected by a uniform velocity read from a ``tectonics: hdisp``
map (m/yr). The advection kernel is driven directly by calling
``getTectonics()`` once per step (first call loads the velocity, every call
advects by ``dt``), so no erosion, deposition or flow routing perturbs the field.
The bump is kept more than 5 widths from every edge so the Dirichlet edges and
the ``fitedges`` edge reset never see it.
"""
from pathlib import Path

import numpy as np
import pytest

scipy = pytest.importorskip("scipy")
matplotlib = pytest.importorskip("matplotlib")
from scipy.spatial import Delaunay  # noqa: E402

from gospl.model import Model  # noqa: E402

# ---- Problem definition (SI, years) ----------------------------------------
DX = 250.0                 # lattice spacing (m)
LX, LY = 60.0e3, 40.0e3    # domain (m)
SIG0 = 3.0e3               # Gaussian half-width (m)
AMP = 100.0                # bump amplitude (m)
SPEED = 0.1                # |v| (m/yr) = 100 km/Myr
TRAVEL = 20.0e3            # distance travelled (m)
X0 = 15.0e3                # start x of the bump centre (m)

SCHEMES = {"upwind": 1, "iioe1": 2, "iioe2": 3}


def _hex_lattice():
    dy = DX * np.sqrt(3.0) / 2.0
    xs = np.arange(0.0, LX + 0.5 * DX, DX)
    nrow = int(LY / dy) + 1
    x = np.concatenate([xs + (0.5 * DX if j % 2 else 0.0) for j in range(nrow)])
    y = np.concatenate([np.full(xs.size, j * dy) for j in range(nrow)])
    return x, y


def _setup(tmp, scheme, angle_deg, courant, field="bump"):
    """Write mesh + YAML; return (yaml, x0, y0, ux, uy, dt, nstep)."""
    x, y = _hex_lattice()
    a = np.radians(angle_deg)
    ux, uy = SPEED * np.cos(a), SPEED * np.sin(a)
    y0 = 0.5 * (y.max() - TRAVEL * np.sin(a))
    if field == "bump":
        z = AMP * np.exp(-0.5 * ((x - X0) ** 2 + (y - y0) ** 2) / SIG0 ** 2)
    else:
        z = np.full(x.size, AMP)
    cells = Delaunay(np.column_stack([x, y])).simplices.astype(np.int64)
    uv = np.zeros((x.size, 3))
    uv[:, 0], uv[:, 1] = ux, uy
    np.savez(tmp / "adv_mesh.npz", v=np.column_stack([x, y, 0.0 * x]),
             c=cells, z=z, uv=uv)

    dt = courant * DX / SPEED
    nstep = int(round(TRAVEL / (SPEED * dt)))
    tend = nstep * dt
    (tmp / "adv.yml").write_text(
        "name: horizontal advection benchmark\n"
        "domain:\n"
        "    npdata: ['adv_mesh','v','c','z']\n"
        "    flowdir: 5\n"
        "    bc: 'oooo'\n"
        f"    advect: '{scheme}'\n"
        "time:\n"
        f"    start: 0.\n    end: {tend}\n    tout: {tend}\n    dt: {dt}\n"
        "spl:\n    K: 0.\n    d: 0.\n    m: 0.5\n"
        "sea:\n    position: -1000.\n"
        "climate:\n  - start: 0.\n    uniform: 1.\n"
        "tectonics:\n"
        f"    - start: 0.\n      end: {tend}\n      hdisp: ['adv_mesh','uv']\n"
        "output:\n    dir: 'adv_out'\n    makedir: True\n"
    )
    return "adv.yml", X0, y0, ux, uy, dt, nstep


def _moments(h, x, y, area, mask):
    """Volume, centroid and weights of the positive part of h."""
    w = np.where(mask, np.maximum(h, 0.0) * area, 0.0)
    vol = w.sum()
    xc, yc = (w * x).sum() / vol, (w * y).sum() / vol
    return vol, xc, yc, w


def _run(tmp, scheme, angle_deg, courant, field="bump", count_excess=False):
    yml, x0, y0, ux, uy, dt, nstep = _setup(tmp, scheme, angle_deg, courant,
                                            field)
    m = Model(yml, verbose=False, showlog=False)
    try:
        assert m.flatModel and m.advscheme == SCHEMES[scheme]
        excess = {"zero": 0, "nonzero": 0}
        if count_excess:
            # Record whether each Scheme-2 call takes the no-overshoot path,
            # with the same overshoot measure the kernel uses; the original
            # method is then called unchanged.
            orig = m._advectorIIOE2

            def wrapped(gvec, vL, newv, vmin, vmax, nbOut):
                d = np.maximum(newv - vmax, 0.0) + np.maximum(vmin - newv, 0.0)
                excess["zero" if d.max() == 0.0 else "nonzero"] += 1
                return orig(gvec, vL, newv, vmin, vmax, nbOut)

            m._advectorIIOE2 = wrapped

        h0 = m.hLocal.getArray().copy()
        for _ in range(nstep):
            m.getTectonics()
        h1 = m.hLocal.getArray().copy()
        out = dict(
            h0=h0, h1=h1, x=m.lcoords[:, 0].copy(), y=m.lcoords[:, 1].copy(),
            area=m.larea.copy(), own=(m.inIDs == 1), x0=x0, y0=y0,
            ux=ux, uy=uy, dt=dt, nstep=nstep, t=nstep * dt, excess=excess,
        )
    finally:
        m.destroy()
    return out


def _metrics(r):
    x, y, area, own = r["x"], r["y"], r["area"], r["own"]
    xe = r["x0"] + r["ux"] * r["t"]
    ye = r["y0"] + r["uy"] * r["t"]
    exact = AMP * np.exp(-0.5 * ((x - xe) ** 2 + (y - ye) ** 2) / SIG0 ** 2)
    _, xc, yc, w = _moments(r["h1"], x, y, area, own)
    # Signed volume conservation (no clipping at 0, so undershoots count).
    m0 = np.sum(r["h0"][own] * area[own])
    m1 = np.sum(r["h1"][own] * area[own])
    # Second moments along / across the velocity direction.
    ang = np.arctan2(r["uy"], r["ux"])
    s = (x - xc) * np.cos(ang) + (y - yc) * np.sin(ang)
    n = -(x - xc) * np.sin(ang) + (y - yc) * np.cos(ang)
    var_s = (w * s ** 2).sum() / w.sum()
    var_n = (w * n ** 2).sum() / w.sum()
    ipk = np.argmax(np.where(own, r["h1"], -np.inf))
    err = r["h1"][own] - exact[own]
    return dict(
        mass_rel=(m1 - m0) / m0,
        centroid_err=float(np.hypot(xc - xe, yc - ye)),
        peak_err=float(np.hypot(x[ipk] - xe, y[ipk] - ye)),
        peak_ratio=float(r["h1"][own].max() / AMP),
        hmin=float(r["h1"][own].min() / AMP),
        rmse=float(np.sqrt(np.mean(err ** 2)) / AMP),
        l2=float(np.sqrt(np.sum(area[own] * err ** 2)
                         / np.sum(area[own] * exact[own] ** 2))),
        D_along=float((var_s - SIG0 ** 2) / (2.0 * r["t"])),
        D_across=float((var_n - SIG0 ** 2) / (2.0 * r["t"])),
        exact=exact,
    )


def _figure(request, name, r, met, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x, y, own = r["x"], r["y"], r["own"]
    ye = r["y0"] + r["uy"] * r["t"]
    band = own & (np.abs(y - ye) < 0.6 * DX)
    order = np.argsort(x[band])
    fig, ax = plt.subplots(figsize=(7, 3.2))
    ax.plot(x[band][order] / 1e3, met["exact"][band][order], "k-", lw=2,
            label="exact (translated)")
    ax.plot(x[band][order] / 1e3, r["h1"][band][order], "C1--", lw=1.5,
            label="goSPL")
    ax.plot(x[band][order] / 1e3, r["h0"][band][order], color="0.6", lw=1,
            label="initial (row through end centre)")
    ax.set_xlabel("x (km)")
    ax.set_ylabel("h (m)")
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    dst = Path(request.config.rootpath) / "results" / "advection"
    dst.mkdir(parents=True, exist_ok=True)
    fig.savefig(dst / f"{name}.png", dpi=120)
    plt.close(fig)


# Courant number C = U dt / dx = 0.5 for the three main cases.
COURANT = 0.5

# Per-scheme acceptance thresholds. Volume is relative; distances in metres;
# peak/min/rmse relative to the amplitude.
TOL = {
    # first-order: strongly diffusive but bounded and conservative
    "upwind": dict(mass=1e-4, centroid=0.5 * DX, peak=(0.60, 0.80),
                   hmin=-1e-10, rmse=0.05),
    # second order: tiny undershoots, almost no amplitude loss
    "iioe1": dict(mass=1e-4, centroid=0.5 * DX, peak=(0.98, 1.0 + 1e-9),
                  hmin=-1e-4, rmse=0.005),
    # limited second order: bounded, some extremum clipping
    "iioe2": dict(mass=2e-3, centroid=1.0 * DX, peak=(0.93, 1.0 + 1e-9),
                  hmin=-1e-10, rmse=0.006),
}


@pytest.mark.benchmark
@pytest.mark.slow
@pytest.mark.parametrize("scheme", ["upwind", "iioe1", "iioe2"])
def test_advection_gaussian_translation(scheme, tmp_path, monkeypatch, request):
    """A Gaussian bump advected 20 km by a uniform velocity vs the exact
    translated profile, per scheme."""
    monkeypatch.chdir(tmp_path)
    r = _run(tmp_path, scheme, 0.0, COURANT)
    met = _metrics(r)
    tol = TOL[scheme]
    d_theory = 0.5 * SPEED * DX * (1.0 + COURANT)
    print(
        f"\n[advection:{scheme}] steps={r['nstep']} dt={r['dt']:.0f} yr "
        f"C={COURANT}  mass={met['mass_rel']:+.2e}  "
        f"centroid err={met['centroid_err']:.1f} m  "
        f"peak err={met['peak_err']:.0f} m  peak={met['peak_ratio']:.4f}  "
        f"min={met['hmin']:+.2e}  RMSE={100 * met['rmse']:.3f}%A  "
        f"L2={100 * met['l2']:.2f}%  D_along={met['D_along']:.2f} "
        f"D_across={met['D_across']:.2f} m2/yr (upwind est. {d_theory:.2f})"
    )
    _figure(request, f"advection_{scheme}", r, met,
            f"{scheme}: Gaussian advected {TRAVEL / 1e3:.0f} km, C={COURANT}")

    assert np.isfinite(r["h1"]).all()
    assert abs(met["mass_rel"]) < tol["mass"], "volume not conserved"
    assert met["centroid_err"] < tol["centroid"], "bump displaced by != v t"
    assert met["peak_err"] <= 2.0 * DX, "peak not at x0 + v t"
    lo, hi = tol["peak"]
    assert lo <= met["peak_ratio"] <= hi, (
        f"peak attenuation {met['peak_ratio']:.4f} outside [{lo}, {hi}]")
    assert met["hmin"] >= tol["hmin"], f"undershoot {met['hmin']:.2e}"
    assert met["rmse"] < tol["rmse"], f"RMSE {100 * met['rmse']:.3f}%A"

    if scheme == "upwind":
        # The widening is the first-order numerical diffusion: compare the
        # measured along-flow diffusivity with the modified-equation estimate
        # (U dx / 2)(1 + C). On a hexagonal lattice the effective spacing is
        # not exactly dx, hence the wide band.
        ratio = met["D_along"] / d_theory
        assert 0.7 < ratio < 1.6, f"upwind D_num/estimate = {ratio:.2f}"
        # The predicted peak decay of a 1-D-stretched Gaussian.
        s_along = np.sqrt(SIG0 ** 2 + 2.0 * met["D_along"] * r["t"])
        s_across = np.sqrt(SIG0 ** 2 + 2.0 * max(met["D_across"], 0.0) * r["t"])
        pred = SIG0 ** 2 / (s_along * s_across)
        assert abs(met["peak_ratio"] - pred) < 0.05, (
            f"peak {met['peak_ratio']:.3f} vs diffusion estimate {pred:.3f}")
    else:
        # Second-order schemes: negligible numerical diffusion.
        assert met["D_along"] < 0.1 * d_theory, "IIOE too diffusive"
    if scheme == "iioe2":
        # The anti-overshoot correction must not create new maxima.
        assert r["h1"][r["own"]].max() <= AMP * (1.0 + 1e-9)


@pytest.mark.benchmark
@pytest.mark.slow
def test_advection_oblique_iioe1(tmp_path, monkeypatch):
    """Oblique (30 deg) velocity: translation must not depend on the lattice
    orientation."""
    monkeypatch.chdir(tmp_path)
    r = _run(tmp_path, "iioe1", 30.0, COURANT)
    met = _metrics(r)
    print(
        f"\n[advection:iioe1 30deg] mass={met['mass_rel']:+.2e}  "
        f"centroid err={met['centroid_err']:.1f} m  peak={met['peak_ratio']:.4f}"
        f"  RMSE={100 * met['rmse']:.3f}%A  D_along={met['D_along']:.2f} "
        f"D_across={met['D_across']:.2f}"
    )
    assert abs(met["mass_rel"]) < 1e-4
    assert met["centroid_err"] < 0.5 * DX
    assert met["peak_ratio"] > 0.98
    assert met["rmse"] < 0.005


@pytest.mark.benchmark
@pytest.mark.slow
def test_advection_iioe2_no_overshoot_field(tmp_path, monkeypatch):
    """A uniform field has no overshoot (neighbourhood min == max), so every
    IIOE2 step takes the excess == 0 path. The exact solution is the field
    itself. Before the 2026-06-29 fix this path returned zeros."""
    monkeypatch.chdir(tmp_path)
    r = _run(tmp_path, "iioe2", 30.0, COURANT, field="uniform",
             count_excess=True)
    own = r["own"]
    err = np.abs(r["h1"][own] - AMP).max()
    print(f"\n[advection:iioe2 uniform] excess calls {r['excess']}  "
          f"max |h - h0| = {err:.2e} m")
    # The no-overshoot path was actually exercised, on every call.
    assert r["excess"]["zero"] > 0 and r["excess"]["nonzero"] == 0
    assert np.isfinite(r["h1"]).all()
    assert err < 1e-6 * AMP, "uniform field not preserved by iioe2"
