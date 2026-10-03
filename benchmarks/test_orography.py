"""
Analytical benchmark: mesh-native orographic precipitation (linear theory).

goSPL's orographic rain (``gospl.tools.addprocess.cptOrography``) is the
Smith and Barstad (2004) linear model with the stratified mountain-wave term
dropped, written as two steady advection-relaxation equations solved on the
unstructured mesh (first-order upwind FV, see ``_buildOroMat``/``_oroSolve``):

    (v . grad + 1/tau_c) q_c = Cw v . grad(h_s)       (cloud water)
    (v . grad + 1/tau_f) q_s = q_c / tau_c            (hydrometeors)
    P_oro = q_s / tau_f                               (kg m^-2 s^-1 = mm/s)

with ``h_s = max(h, sea level)`` (no forcing from submarine bathymetry) and
``Cw = ref_density * moist_lapse_rate / env_lapse_rate`` (``oro_cw``).

For a 1-D ridge h(x) (uniform across the wind) and a uniform wind ``u`` along x
the pair is linear with constant coefficients, so in Fourier space

    P_oro_hat(k) = Cw i k u h_hat(k) / ((1 + i k u tau_c)(1 + i k u tau_f))

which is evaluated here with a dense, heavily padded FFT (the input is a
compactly supported clamped ridge, so periodic wrap-around is negligible). The
sign of ``u`` carries the wind direction; ``wind_dir`` is the bearing the wind
comes FROM, so ``wind_dir=270`` blows toward +x and ``wind_dir=90`` toward -x.

goSPL then converts and clips (mirrored exactly in the reference):

    P [mm/h]   = max(3600 P_oro + precip_base, precip_min)
    rain [m/yr] = 0.366 * rainfall_frequency * P

The comparison is done on the full profile, clipping included, AND separately
on the raw orographic term ``q_s / tau_f`` (``_oroQs``), which isolates the
discretisation error from the post-processing.

Mesh: flat hexagonal lattice (rows offset by dx/2) built at runtime in
``tmp_path``. Only the central rows are compared: every domain edge is a zero
Dirichlet row for q_c and q_s (``advectBorders``), and on a hexagonal lattice the
faces at +-60 deg to the wind let the zero lateral (wind-parallel) edges bleed
into the first few km of rows next to them (measured: ~36% RMSE on the edge
row, decaying to the 0.4% interior level ~3 km in). Upstream of the coast the
zero inflow edge is exact (no forcing there), and the outflow edge does not
feed back under upwinding, so only the wind-parallel edges matter.
"""
from pathlib import Path

import numpy as np
import pytest

scipy = pytest.importorskip("scipy")
matplotlib = pytest.importorskip("matplotlib")
from scipy.spatial import Delaunay  # noqa: E402

from gospl.model import Model  # noqa: E402

# ---- Problem definition (SI) -----------------------------------------------
LX = 150.0e3          # domain length along the wind (m)
NROW = 41             # lattice rows across the wind
X0 = 75.0e3           # ridge crest (m)
SIG = 8.0e3           # ridge half-width (m)
H0 = 1800.0           # ridge amplitude above the abyssal base (m)
ZBASE = -300.0        # far-field bathymetry (m): only the ridge top is land
SEA = 0.0             # sea level (m)
WIND = 10.0           # wind speed (m/s)
TAU_C = 500.0         # conversion time (s)  -> L_c = u tau_c = 5 km
TAU_F = 500.0         # fall-out time (s)
P_BASE = 2.0          # background precipitation (mm/h)
P_MIN = 0.01          # floor (mm/h)
FREQ = 1.0            # rainfall_frequency


def _ridge(x):
    return ZBASE + H0 * np.exp(-0.5 * ((x - X0) / SIG) ** 2)


def _hex_lattice(dx):
    dy = dx * np.sqrt(3.0) / 2.0
    xs = np.arange(0.0, LX + 0.5 * dx, dx)
    x = np.concatenate([xs + (0.5 * dx if j % 2 else 0.0)
                        for j in range(NROW)])
    y = np.concatenate([np.full(xs.size, j * dy) for j in range(NROW)])
    return x, y, dy


def _write_case(tmp, dx, wind_dir):
    x, y, dy = _hex_lattice(dx)
    cells = Delaunay(np.column_stack([x, y])).simplices.astype(np.int64)
    np.savez(tmp / "oro_mesh.npz", v=np.column_stack([x, y, 0.0 * x]),
             c=cells, z=_ridge(x))
    (tmp / "oro.yml").write_text(
        "name: orographic precipitation benchmark\n"
        "domain:\n"
        "    npdata: ['oro_mesh','v','c','z']\n"
        "    flowdir: 5\n"
        "    bc: 'oooo'\n"
        "time:\n    start: 0.\n    end: 10.\n    tout: 10.\n    dt: 10.\n"
        "spl:\n    K: 0.\n    d: 0.\n    m: 0.5\n"
        f"sea:\n    position: {SEA}\n"
        "orography:\n"
        f"    wind_speed: {WIND}\n"
        f"    wind_dir: {wind_dir}\n"
        f"    conv_time: {TAU_C}\n"
        f"    fall_time: {TAU_F}\n"
        f"    precip_base: {P_BASE}\n"
        f"    precip_min: {P_MIN}\n"
        f"    rainfall_frequency: {FREQ}\n"
        "output:\n    dir: 'oro_out'\n    makedir: True\n"
    )
    return dy


def exact_orographic(x, u, cw, n=2 ** 17, pad=8.0):
    """Exact P_oro(x) (mm/h) for the clamped 1-D ridge and signed wind u."""
    length = pad * LX
    xf = (np.arange(n) * length / n) - 0.5 * (length - LX)
    hs = np.maximum(_ridge(xf), SEA)
    hs -= hs[0]                     # constant offset carries no forcing
    k = 2.0 * np.pi * np.fft.fftfreq(n, length / n)
    ph = cw * 1j * k * u * np.fft.fft(hs) / (
        (1.0 + 1j * k * u * TAU_C) * (1.0 + 1j * k * u * TAU_F))
    pf = np.real(np.fft.ifft(ph)) * 3600.0
    return np.interp(x, xf, pf)


def exact_rain(p_oro):
    """goSPL post-processing: background, floor clip, mm/h -> m/yr."""
    return np.maximum(p_oro + P_BASE, P_MIN) * 0.366 * FREQ


def _run(tmp, dx, wind_dir):
    dy = _write_case(tmp, dx, wind_dir)
    m = Model("oro.yml", verbose=False, showlog=False)
    try:
        assert m.oroOn and m.flatModel
        m.cptOrography()
        m.dm.globalToLocal(m._oroQs, m.tmpL)
        p_oro = m.tmpL.getArray().copy() / m.oro_fall_time * 3600.0
        u, v = m._windVector()
        out = dict(
            x=m.lcoords[:, 0].copy(), y=m.lcoords[:, 1].copy(),
            rain=np.asarray(m.rainVal).copy(), p_oro=p_oro,
            own=(m.inIDs == 1), cw=m.oro_cw, u=u, v=v, dy=dy,
            border=np.isin(np.arange(m.lpoints), m.advectBorders),
        )
    finally:
        m.destroy()
    return out


def _compare(r, dx):
    x, y = r["x"], r["y"]
    ymid = 0.5 * (y.max() + y.min())
    band = (r["own"] & ~r["border"] & (np.abs(y - ymid) < 2.6 * r["dy"])
            & (x > dx) & (x < LX - dx))
    ref_oro = exact_orographic(x, r["u"], r["cw"])
    ref_rain = exact_rain(ref_oro)
    pk = np.abs(ref_oro[band]).max()
    e_oro = r["p_oro"][band] - ref_oro[band]
    e_rain = r["rain"][band] - ref_rain[band]
    xb = x[band]
    clipped_ref = band & (ref_oro + P_BASE <= P_MIN)
    clipped = band & (r["p_oro"] + P_BASE <= P_MIN)
    land = r["own"] & (_ridge(x) > SEA)
    # windward ocean: upwind of every land node
    if r["u"] > 0:
        windward = band & (x < x[land].min() - 2 * dx)
    else:
        windward = band & (x > x[land].max() + 2 * dx)
    return dict(
        band=band, ref_oro=ref_oro, ref_rain=ref_rain,
        rmse_oro=float(np.sqrt(np.mean(e_oro ** 2)) / pk),
        max_oro=float(np.abs(e_oro).max() / pk),
        rmse_rain=float(np.sqrt(np.mean(e_rain ** 2))
                        / ref_rain[band].max()),
        xpk_num=float(xb[np.argmax(r["p_oro"][band])]),
        xpk_ref=float(xb[np.argmax(ref_oro[band])]),
        pk_ratio=float(r["p_oro"][band].max() / ref_oro[band].max()),
        min_ratio=float(r["p_oro"][band].min() / ref_oro[band].min()),
        n_clipped=int(clipped.sum()),
        clip_overlap=float((clipped & clipped_ref).sum()
                           / max((clipped | clipped_ref).sum(), 1)),
        clip_err=float(np.abs(r["rain"][clipped]
                              - P_MIN * 0.366 * FREQ).max())
        if clipped.any() else 0.0,
        windward_err=float(np.abs(r["rain"][windward]
                                  - P_BASE * 0.366 * FREQ).max()),
        n_windward=int(windward.sum()),
    )


def _figure(request, name, r, c, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    b = c["band"]
    o = np.argsort(r["x"][b])
    xs = r["x"][b][o] / 1e3
    fig, ax = plt.subplots(2, 1, figsize=(7, 5), sharex=True)
    ax[0].plot(xs, np.maximum(_ridge(r["x"][b][o]), SEA), "k-")
    ax[0].plot(xs, _ridge(r["x"][b][o]), color="0.6", lw=0.8)
    ax[0].set_ylabel("z (m)")
    ax[1].plot(xs, c["ref_rain"][b][o], "k-", lw=2, label="exact")
    ax[1].plot(xs, r["rain"][b][o], "C0--", lw=1.5, label="goSPL")
    ax[1].set_ylabel("rain (m/yr)")
    ax[1].set_xlabel("x (km)")
    ax[1].legend(fontsize=7, frameon=False)
    ax[0].set_title(title, fontsize=9)
    fig.tight_layout()
    dst = Path(request.config.rootpath) / "results" / "orography"
    dst.mkdir(parents=True, exist_ok=True)
    fig.savefig(dst / f"{name}.png", dpi=120)
    plt.close(fig)


DX = 250.0


@pytest.mark.benchmark
@pytest.mark.slow
@pytest.mark.parametrize("wind_dir", [270.0, 90.0])
def test_orography_ridge_exact(wind_dir, tmp_path, monkeypatch, request):
    """Orographic rain across a Gaussian coastal ridge vs the exact linear
    solution (clipping and background included)."""
    monkeypatch.chdir(tmp_path)
    r = _run(tmp_path, DX, wind_dir)
    c = _compare(r, DX)
    print(
        f"\n[orography wind_dir={wind_dir:.0f}] u={r['u']:+.1f} m/s  "
        f"P_oro RMSE={100 * c['rmse_oro']:.2f}% max={100 * c['max_oro']:.2f}%"
        f" of peak  rain RMSE={100 * c['rmse_rain']:.2f}%  "
        f"peak x num/ref={c['xpk_num'] / 1e3:.2f}/{c['xpk_ref'] / 1e3:.2f} km"
        f"  peak ratio={c['pk_ratio']:.4f}  lee-min ratio={c['min_ratio']:.4f}"
        f"  clipped={c['n_clipped']} (overlap {c['clip_overlap']:.3f})"
        f" windward={c['n_windward']}"
    )
    _figure(request, f"orography_{int(wind_dir)}", r, c,
            f"Gaussian ridge, wind_dir={wind_dir:.0f}, dx={DX:.0f} m")

    own = r["own"]
    assert np.isfinite(r["rain"][own]).all() and (r["rain"][own] > 0).all()
    # Wind direction convention: from the West blows +x, from the East -x.
    assert np.sign(r["u"]) == (1.0 if wind_dir == 270.0 else -1.0)
    # Raw orographic term vs exact (first-order upwind, dx/L_c = 0.05).
    assert c["rmse_oro"] < 0.015, f"P_oro RMSE {100 * c['rmse_oro']:.2f}%"
    assert c["max_oro"] < 0.05
    # Peak location and magnitude, lee (rain-shadow) minimum.
    assert abs(c["xpk_num"] - c["xpk_ref"]) <= 2.0 * DX
    assert abs(c["pk_ratio"] - 1.0) < 0.03
    assert abs(c["min_ratio"] - 1.0) < 0.05
    # Full rainfall field (background + clip + unit conversion).
    assert c["rmse_rain"] < 0.01, f"rain RMSE {100 * c['rmse_rain']:.2f}%"
    # The lee clip is exercised, lands on the floor exactly, and covers the
    # same nodes as the reference clip (to within the edge of the dry zone).
    assert c["n_clipped"] > 0 and c["clip_err"] < 1e-12
    assert c["clip_overlap"] > 0.95
    # Sea-level clamp: upwind of the coast (submarine slope) there is no
    # orographic forcing at all, so the rain is exactly the background.
    assert c["n_windward"] > 0 and c["windward_err"] < 1e-9


@pytest.mark.benchmark
@pytest.mark.slow
def test_orography_first_order_convergence(tmp_path, monkeypatch):
    """Halving dx should roughly halve the error (first-order upwind)."""
    errs = {}
    for dx in (500.0, 250.0):
        sub = tmp_path / f"dx{int(dx)}"
        sub.mkdir()
        monkeypatch.chdir(sub)
        r = _run(sub, dx, 270.0)
        errs[dx] = _compare(r, dx)["rmse_oro"]
    order = np.log2(errs[500.0] / errs[250.0])
    print(f"\n[orography convergence] RMSE dx=500: {100 * errs[500.0]:.3f}%  "
          f"dx=250: {100 * errs[250.0]:.3f}%  observed order {order:.2f}")
    assert 0.7 < order < 1.5, f"unexpected convergence order {order:.2f}"
