"""
Analytical benchmarks -- flexural isostasy (thin elastic plate / shell).

Validates both goSPL flexure solvers against exact solutions of the thin-plate
equation on an elastic (inviscid-fluid) foundation,

    D lap^2 w + drho g w = q,      D = E Te^3 / (12 (1 - nu^2)),

with q = rho_s g dz the load of a surface thickness change dz and drho = rho_a
(goSPL takes the moat infill density as 0). goSPL returns the deflection with
the sign convention NEGATIVE = subsidence under deposition.

Every test builds its own mesh + YAML in ``tmp_path`` (no committed binary
meshes) and drives the model's own entry point ``applyFlexure`` with a
prescribed load: ``hOldFlex`` is set to ``h - dz`` so the load seen by the
solver is exactly ``dz``, and the response is read from ``localFlex``. This
exercises the full production path (load assembly, BCs, solve, scatter) without
running erosion, so the comparison is free of any other process.

1. Flat model, ``flexure: method: fem`` (parallel FV biharmonic,
   ``addprocess._buildFlexFEM`` / ``_cmptFlexFEM``)
   --------------------------------------------------------------------------
   The natural boundary of the FV operator ``Lm diag(D) Lm`` is the
   ``0Slope0Shear`` condition (w' = 0 and w''' = 0 on every side). Two loads
   have an EXACT solution under exactly that condition on a bounded box, so the
   boundary is part of the test instead of something to stay away from:

   (a) Cosine modes ``q = q0 cos(kx x) cos(ky y)`` with ``kx = n pi/Lx``,
       ``ky = m pi/Ly``: every derivative of odd order vanishes on the edges,
       so the cosine is the exact solution on the box,
           w = q0 cos(kx x) cos(ky y) / (D k^4 + drho g),   k^2 = kx^2 + ky^2.
       On a square-lattice Delaunay mesh the Voronoi cells are squares, so the
       FV Laplacian is the 5-point stencil (diagonal edges have zero dual
       length) and the cosine is ALSO an exact eigenvector of the discrete
       operator with eigenvalue ``lam = (4/h^2)[sin^2(kx h/2) + sin^2(ky h/2)]``.
       That gives two checks: goSPL must reproduce the discrete solution to
       round-off (assembly/BC/solve are exact), and the continuum solution with
       an O((k h)^2) truncation error that falls 4x per halving of h.

   (b) A line load of V N/m along x = xc. The infinite-plate solution
           w(s) = w0 exp(-|s|/a) (cos(|s|/a) + sin(|s|/a)),
           a = (4D/(drho g))^(1/4),   w0 = V a^3 / (8 D)
       is the classic Hetenyi (1946) / Turcotte & Schubert line-load result. On the
       box [0, L] a 0Slope0Shear edge is a mirror, so the bounded-domain
       solution is the image sum  sum_j w(x - xc - j L)  -- exact, no "far
       from the boundary" approximation needed. This checks the flexural
       wavelength, the forebulge and the peak amplitude together.

   ``0Displacement0Slope`` (clamped) has no comparably clean closed form for
   these loads and is covered by the regression suite instead.

2. Global model, ``flexure: method: global`` (spherical harmonics on rank 0,
   ``_buildDHGrid`` / ``_cmptFlexGlobal``)
   --------------------------------------------------------------------------
   For a load proportional to a single spherical harmonic of degree l, the
   operator is diagonal in degree: with L = l(l+1),
       w_l = q_l / (drho g + D P_l),   P_l = (L^2 - 4L) / R^4,
   (degrees 0 and 1 removed), which is goSPL's documented thin-shell bending
   operator ``D lap(lap + 4/R^2)``. Two checks:

   (c) Spectral core on the Driscoll-Healy grid: the load is evaluated
       ANALYTICALLY (scipy ``lpmv``, independent of pyshtools) on the DH grid,
       which samples a band-limited field exactly, so ``_cmptFlexGlobal`` must
       return the degree-l response to round-off for every tested (l, m). This
       pins degree indexing, the l = 0/1 removal, the DH/DH2 grids and units.
       It also checks the flat-plate limit: for l >> 1 the shell operator must
       tend to D k^4 with k^2 = L / R^2 (independent physics, not the code's
       own formula).
   (d) Full mesh pipeline (mesh -> DH inverse-distance -> SH solve -> bilinear
       back to the mesh) through ``applyFlexure`` on a Fibonacci sphere. Here
       the error is the interpolation error, second order in grid/mesh spacing.

   Note on the physics (not asserted): goSPL's shell term is bending-only.
   The full thin-shell response (e.g. Turcotte et al. 1981; Wieczorek 2007)
   has a bending term D L^2 (L - 4) / (R^4 (L - 1 + nu)) -- differing from
   goSPL's by the factor L / (L - 1 + nu), i.e. 12% at l = 2, 4% at l = 4,
   <1% for l >= 9 -- plus a membrane term E Te (L - 2) / (R^2 (L - 1 + nu)),
   which is a few % of drho g for Earth-like Te. Both only matter at the very
   lowest degrees of a stiff lithosphere; the benchmark prints the size of
   the difference for the tested degrees.
"""
import shutil
from pathlib import Path

import numpy as np
import pytest

scipy = pytest.importorskip("scipy")
matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.spatial import ConvexHull, Delaunay  # noqa: E402
from scipy.special import lpmv  # noqa: E402

from petsc4py import PETSc  # noqa: E402
from gospl.model import Model  # noqa: E402

# ---- Shared lithosphere parameters (SI) ------------------------------------
RHO_S = 2300.0           # load (crust/sediment) density
RHO_A = 3300.0           # asthenosphere density; drho = rho_a (no infill)
YOUNG = 65.0e9           # Young's modulus
NU = 0.25                # Poisson ratio

# Flat model: Te = 20 km -> flexural parameter a ~ 49 km.
TE_FLAT = 20.0e3
L_BOX = 400.0e3          # square box for the cosine modes
# Global model: Te = 300 km on Earth -> D k^4 ~ drho g at l ~ 24, so the tested
# degrees span the isostatic -> rigid transition.
TE_GLOB = 300.0e3
R_EARTH = 6378137.0
N_SPHERE = 40000         # Fibonacci-sphere nodes (~113 km spacing)
RES_DEG = 0.5            # DH grid resolution -> lmax = 179

# ---- Tolerances ------------------------------------------------------------
# (a) discrete identity: the solve is a direct LU, so round-off only.
TOL_DISCRETE = 1.0e-8
# (a) continuum: 5-point truncation ~ (k h)^2 / 12 on the D k^4 term. Asserted
#     on the finer 101^2 mesh (h = 4 km), where the stiffest mode (6,0)
#     (k h ~ 0.19, D k^4 ~ 11 drho g) is ~0.5% off; the 51^2 mesh only enters the
#     second-order convergence check (error ratio >= 3 per halving of h).
TOL_COS_MAX_PCT = 1.0
# (b) line load: h = 2 km against a = 49 km (a/h ~ 24) -> ~0.05% measured; the
#     bounds (0.5% RMSE, 1% max) leave room for the discretised delta load.
TOL_LINE_RMSE_PCT = 0.5
TOL_LINE_MAX_PCT = 1.0
# (c) spectral core: the load is exactly band-limited on the DH grid.
TOL_SPECTRAL = 1.0e-9
# (d) mesh pipeline: interpolation error ~ (l * spacing / R)^2, ~0.3% RMSE /
#     ~2% max (at the pole of a zonal l = 24 mode) on a 40k-node sphere.
TOL_GLOB_RMSE_PCT = 1.0
TOL_GLOB_MAX_PCT = 4.0

FLAT_YAML = """name: flexure analytical benchmark (flat fem)
domain:
    npdata: ['{mesh}','v','c','z']
    flowdir: 2
    bc: 'ffff'
    nodep: True
time:
    start: 0.
    end: 10.
    tout: 10.
    dt: 10.
spl:
    K: 1.0e-20
    d: 0.
    m: 0.5
flexure:
    method: 'fem'
    thick: {te}
    rhoc: {rhos}
    rhoa: {rhoa}
    young: {young}
    nu: {nu}
    bcN: '0Slope0Shear'
    bcS: '0Slope0Shear'
    bcE: '0Slope0Shear'
    bcW: '0Slope0Shear'
sea:
    position: -1000.
climate:
  - start: 0.
    uniform: 1.
output:
    dir: '{mesh}_out'
"""

GLOBAL_YAML = """name: flexure analytical benchmark (global)
domain:
    npdata: ['{mesh}','v','c','z']
    flowdir: 2
    nodep: True
    radius: {radius}
time:
    start: 0.
    end: 10.
    tout: 10.
    dt: 10.
spl:
    K: 1.0e-20
    d: 0.
    m: 0.5
flexure:
    method: 'global'
    thick: {te}
    rhoc: {rhos}
    rhoa: {rhoa}
    young: {young}
    nu: {nu}
    res_deg: {res}
sea:
    position: -10000.
climate:
  - start: 0.
    uniform: 1.
output:
    dir: '{mesh}_out'
"""


def _rigidity(te):
    return YOUNG * te ** 3 / (12.0 * (1.0 - NU ** 2))


def _errors(w, wa):
    """(RMSE %, max %) of w - wa relative to max |wa|."""
    scale = float(np.abs(wa).max())
    e = w - wa
    return (100.0 * float(np.sqrt(np.mean(e ** 2))) / scale,
            100.0 * float(np.abs(e).max()) / scale)


def _grid_model(tmp_path, name, nx, ny, h):
    """Square-lattice Delaunay mesh (Voronoi = squares) + flat 'fem' model."""
    xs = np.arange(nx) * h
    ys = np.arange(ny) * h
    X, Y = np.meshgrid(xs, ys)
    x, y = X.ravel(), Y.ravel()
    v = np.column_stack([x, y, np.zeros_like(x)])
    cells = Delaunay(np.column_stack([x, y])).simplices.astype(np.int64)
    np.savez(tmp_path / f"{name}.npz", v=v, c=cells, z=np.full(x.size, 500.0))
    (tmp_path / f"{name}.yml").write_text(FLAT_YAML.format(
        mesh=name, te=TE_FLAT, rhos=RHO_S, rhoa=RHO_A, young=YOUNG, nu=NU))
    return Model(f"{name}.yml", verbose=False, showlog=False)


def _apply_load(m, dz):
    """Flex the model under the thickness load dz via the production entry
    point; return the deflection and restore the elevation."""
    hl = m.hLocal.getArray().copy()
    m.hOldFlex.setArray(hl - dz)            # load seen by applyFlexure == dz
    m.localFlex[:] = 0.0
    m.applyFlexure()
    w = m.localFlex.copy()
    m.hLocal.setArray(hl)
    m.dm.localToGlobal(m.hLocal, m.hGlobal)
    return w


def _gather(m, *arrays):
    """Owned values of local arrays, concatenated over all ranks (identical on
    every rank), so error metrics / profiles are global at any np."""
    from mpi4py import MPI
    own = m.inIDs == 1
    out = []
    for a in arrays:
        parts = MPI.COMM_WORLD.allgather(np.ascontiguousarray(a[own]))
        out.append(np.concatenate(parts))
    return out if len(out) > 1 else out[0]


def _save(fig, tmp_path, request, name):
    fig.savefig(tmp_path / name, dpi=110, bbox_inches="tight")
    plt.close(fig)
    dst = Path(request.config.rootpath) / "results" / "flexure"
    dst.mkdir(parents=True, exist_ok=True)
    shutil.copy2(tmp_path / name, dst / name)


# ---------------------------------------------------------------------------
# 1a. Flat fem: cosine modes (exact on the 0Slope0Shear box)
# ---------------------------------------------------------------------------
@pytest.mark.benchmark
@pytest.mark.slow
def test_flexure_fem_cosine_modes(tmp_path, monkeypatch):
    """FV biharmonic vs exact cosine-mode deflection; 2nd-order convergence."""
    monkeypatch.chdir(tmp_path)
    D = _rigidity(TE_FLAT)
    kf = RHO_A * 9.81
    modes = [(1, 0), (2, 0), (4, 0), (4, 2), (3, 3), (6, 0)]
    amp = 100.0
    max_err = {}
    for nx in (51, 101):
        h = L_BOX / (nx - 1)
        m = _grid_model(tmp_path, f"cos{nx}", nx, nx, h)
        try:
            assert m.flatModel and m.flexOn and m.flex_method == "fem"
            assert m.gravity == 9.81
            x, y = m.lcoords[:, 0], m.lcoords[:, 1]
            print(f"\n[flex-fem cos] {nx}x{nx}, h={h / 1e3:.0f} km, "
                  f"a={(4 * D / kf) ** 0.25 / 1e3:.1f} km")
            for n, k in modes:
                kx, ky = n * np.pi / L_BOX, k * np.pi / L_BOX
                dz = amp * np.cos(kx * x) * np.cos(ky * y)
                w = _apply_load(m, dz)
                q = RHO_S * 9.81 * dz
                wa = -q / (D * (kx ** 2 + ky ** 2) ** 2 + kf)
                lam = 4.0 * (np.sin(kx * h / 2) ** 2 + np.sin(ky * h / 2) ** 2) / h ** 2
                wd = -q / (D * lam ** 2 + kf)
                w, wa, wd = _gather(m, w, wa, wd)
                rel_disc = float(np.abs(w - wd).max() / np.abs(wd).max())
                rmse, emax = _errors(w, wa)
                max_err[(nx, n, k)] = emax
                print(f"  mode ({n},{k}): Dk^4/drho g={D * (kx**2 + ky**2)**2 / kf:6.3f}"
                      f"  RMSE={rmse:.4f}%  max={emax:.4f}%  |w-w_discrete|={rel_disc:.1e}")
                assert rel_disc < TOL_DISCRETE, (
                    f"mode ({n},{k}): not the exact discrete solution "
                    f"(rel {rel_disc:.2e}) -- assembly/BC/solve error")
                assert nx < 101 or emax < TOL_COS_MAX_PCT, (
                    f"mode ({n},{k}) at {nx}^2: max error {emax:.3f}% "
                    f"> {TOL_COS_MAX_PCT}%")
        finally:
            m.destroy()
    # Second-order convergence on the stiffest modes (error falls ~4x per
    # halving of h; accept >= 3x).
    for n, k in [(4, 0), (6, 0), (4, 2)]:
        ratio = max_err[(51, n, k)] / max_err[(101, n, k)]
        print(f"  convergence ({n},{k}): 51->101 error ratio {ratio:.2f}")
        assert ratio > 3.0, f"mode ({n},{k}) not second-order: ratio {ratio:.2f}"


# ---------------------------------------------------------------------------
# 1b. Flat fem: line load (Hetenyi), exact via mirror images
# ---------------------------------------------------------------------------
@pytest.mark.benchmark
@pytest.mark.slow
def test_flexure_fem_line_load(tmp_path, monkeypatch, request):
    """FV biharmonic vs Hetenyi line-load profile (peak, wavelength, bulge)."""
    monkeypatch.chdir(tmp_path)
    nx, ny, h = 601, 11, 2000.0             # 1200 km strip, a/h ~ 24
    m = _grid_model(tmp_path, "line", nx, ny, h)
    try:
        D = _rigidity(TE_FLAT)
        kf = RHO_A * m.gravity
        a = (4.0 * D / kf) ** 0.25
        L = (nx - 1) * h
        xc = 0.5 * L
        x = m.lcoords[:, 0]
        dz0 = 1000.0
        dz = np.where(np.abs(x - xc) < 0.25 * h, dz0, 0.0)  # one node column
        V = dz0 * h * RHO_S * m.gravity    # N/m: load per unit length of line
        w = _apply_load(m, dz)

        def w_inf(s):
            s = np.abs(s)
            return V * a ** 3 / (8.0 * D) * np.exp(-s / a) * (np.cos(s / a) + np.sin(s / a))

        # 0Slope0Shear edges at x=0, L are mirrors -> image sum is exact.
        wa = -sum(w_inf(x - xc - j * L) for j in range(-4, 5))
        x, w, wa = _gather(m, x, w, wa)
        rmse, emax = _errors(w, wa)
        bulge, bulge_a = float(w.max()), float(wa.max())
        # first zero crossing (analytic: s = 3 pi a / 4)
        xs_ = np.sort(np.unique(np.round(x, 6)))
        prof = np.array([w[np.isclose(x, xv)].mean() for xv in xs_])
        right = xs_ >= xc
        sgn = np.sign(prof[right])
        i0 = int(np.argmax(sgn > 0))
        x0 = xs_[right][i0 - 1] - prof[right][i0 - 1] * h / (prof[right][i0] - prof[right][i0 - 1])
        s0, s0_a = x0 - xc, 0.75 * np.pi * a
        print(f"\n[flex-fem line] a={a / 1e3:.1f} km  w0={abs(wa).max():.3f} m  "
              f"peak goSPL={w.min():.4f} ana={wa.min():.4f}  bulge goSPL={bulge:.4f} "
              f"ana={bulge_a:.4f}  zero-crossing {s0 / 1e3:.2f} vs {s0_a / 1e3:.2f} km  "
              f"RMSE={rmse:.4f}%  max={emax:.4f}%")

        fig, ax = plt.subplots(figsize=(7, 3.5))
        sel = np.abs(xs_ - xc) < 6 * a
        ax.plot((xs_[sel] - xc) / 1e3, prof[sel], "o", ms=2.5, label="goSPL fem")
        sx = np.linspace(-6 * a, 6 * a, 1000)
        ax.plot(sx / 1e3, -sum(w_inf(sx - j * L) for j in range(-4, 5)), "k-", lw=1,
                label="Hetenyi (exact)")
        ax.set_xlabel("distance from line load (km)")
        ax.set_ylabel("deflection (m)")
        ax.set_title(f"Line load, Te={TE_FLAT / 1e3:.0f} km  (RMSE {rmse:.3f}%)")
        ax.legend()
        if PETSc.COMM_WORLD.getRank() == 0:
            _save(fig, tmp_path, request, "flexure_fem_line_load.png")
        else:
            plt.close(fig)

        assert rmse < TOL_LINE_RMSE_PCT, f"line-load RMSE {rmse:.3f}% > {TOL_LINE_RMSE_PCT}%"
        assert emax < TOL_LINE_MAX_PCT, f"line-load max error {emax:.3f}% > {TOL_LINE_MAX_PCT}%"
        assert abs(bulge - bulge_a) < 0.02 * bulge_a, "forebulge amplitude off by >2%"
        assert abs(s0 - s0_a) < 0.02 * s0_a, "flexural zero crossing off by >2%"
    finally:
        m.destroy()


# ---------------------------------------------------------------------------
# 2. Global (spherical harmonic) flexure
# ---------------------------------------------------------------------------
def _sphere_model(tmp_path):
    """Fibonacci sphere at Earth radius + 'global' flexure model."""
    n = N_SPHERE
    i = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    th = np.pi * (1.0 + 5.0 ** 0.5) * i
    v = R_EARTH * np.column_stack([np.cos(th) * np.sin(phi),
                                   np.sin(th) * np.sin(phi), np.cos(phi)])
    cells = ConvexHull(v).simplices.astype(np.int64)   # spherical Delaunay
    np.savez(tmp_path / "sphere.npz", v=v, c=cells, z=np.full(n, 500.0))
    (tmp_path / "sphere.yml").write_text(GLOBAL_YAML.format(
        mesh="sphere", radius=R_EARTH, te=TE_GLOB, rhos=RHO_S, rhoa=RHO_A,
        young=YOUNG, nu=NU, res=RES_DEG))
    return Model("sphere.yml", verbose=False, showlog=False)


def _ylm(l, mo, lat, lon):
    """Real (unnormalised) spherical harmonic P_l^m(sin lat) cos(m lon), max 1."""
    y = lpmv(mo, l, np.sin(lat)) * np.cos(mo * lon)
    return y / np.abs(y).max()


GLOBAL_MODES = [(2, 0), (4, 0), (8, 0), (8, 3), (12, 0), (16, 0), (16, 5),
                (24, 0), (24, 12)]


@pytest.mark.benchmark
@pytest.mark.slow
def test_flexure_global_spherical_harmonics(tmp_path, monkeypatch, request):
    """SH thin-shell flexure: exact per-degree response (spectral core) and
    the full mesh -> DH -> mesh pipeline within interpolation tolerance."""
    pytest.importorskip("pyshtools")
    from mpi4py import MPI
    monkeypatch.chdir(tmp_path)
    m = _sphere_model(tmp_path)
    try:
        assert not m.flatModel and m.flexOn and m.flex_method == "global"
        D = _rigidity(TE_GLOB)
        kf = RHO_A * m.gravity
        R = m.radius

        def resp(l):                       # goSPL's documented shell operator
            L = l * (l + 1)
            return 1.0 / (kf + D * (L * L - 4 * L) / R ** 4)

        # ---- (c) spectral core on the DH grid (rank 0 owns the DH grid) ----
        if MPI.COMM_WORLD.Get_rank() == 0:
            lon2, lat2 = np.meshgrid(np.deg2rad(m.dh_lon), np.deg2rad(m.dh_lat))
            print(f"\n[flex-global] DH {m.dh_N}x{2 * m.dh_N} (lmax {m.dh_N // 2 - 1}),"
                  f" Te={TE_GLOB / 1e3:.0f} km, mesh {m.mpoints} nodes")
            for l, mo in GLOBAL_MODES + [(60, 0), (120, 7)]:
                dz = 100.0 * _ylm(l, mo, lat2, lon2)
                w = m._cmptFlexGlobal(dz, float(m.flex_eet))
                wa = -RHO_S * m.gravity * dz * resp(l)
                rel = float(np.abs(w - wa).max() / np.abs(wa).max())
                print(f"  core (l={l:3d}, m={mo:2d}): |w-w_exact|/|w|={rel:.1e}")
                assert rel < TOL_SPECTRAL, (
                    f"spectral core: degree {l} order {mo} response off by {rel:.2e}")
            # Degree 0 and 1 must be removed (no mean / centre-of-mass motion).
            for l in (0, 1):
                w = m._cmptFlexGlobal(100.0 * _ylm(l, 0, lat2, lon2), float(m.flex_eet))
                assert np.abs(w).max() < 1.0e-9, f"degree {l} not removed"
        # Flat-plate limit (independent physics): k^2 = L / R^2 for l >> 1.
        for l in (60, 120):
            L = l * (l + 1)
            plate = 1.0 / (kf + D * (L / R ** 2) ** 2)
            dev = abs(resp(l) - plate) / plate
            assert dev < 4.0 / L + 1e-12, f"l={l}: shell not -> plate limit ({dev:.2e})"

        # ---- (d) full pipeline through applyFlexure on the mesh -----------
        xyz = m.lcoords
        r = np.linalg.norm(xyz, axis=1)
        lat, lon = np.arcsin(xyz[:, 2] / r), np.arctan2(xyz[:, 1], xyz[:, 0])
        rows, ls, ratio = [], [], []
        for l, mo in GLOBAL_MODES:
            dz = 100.0 * _ylm(l, mo, lat, lon)
            w = _apply_load(m, dz)
            wa = -RHO_S * m.gravity * dz * resp(l)
            w, wa = _gather(m, w, wa)
            rmse, emax = _errors(w, wa)
            amp = float(np.dot(w, wa) / np.dot(wa, wa))
            L = l * (l + 1)
            # size of the (unasserted) bending-term difference vs full shell theory
            full = 1.0 / (kf + D * L * L * (L - 4) / (R ** 4 * (L - 1 + NU)))
            shell_diff = 100.0 * (resp(l) - full) / full
            flexfac = D * (L * L - 4 * L) / R ** 4 / kf
            print(f"  mesh (l={l:2d}, m={mo:2d}): D P_l/drho g={flexfac:6.3f}  "
                  f"RMSE={rmse:.3f}%  max={emax:.3f}%  amplitude ratio={amp:.4f}  "
                  f"[goSPL vs full-shell bending: {shell_diff:+.2f}%]")
            rows.append((l, mo, rmse, emax))
            ls.append(l)
            ratio.append(amp * resp(l))
            assert rmse < TOL_GLOB_RMSE_PCT, (
                f"(l={l}, m={mo}) RMSE {rmse:.3f}% > {TOL_GLOB_RMSE_PCT}%")
            assert emax < TOL_GLOB_MAX_PCT, (
                f"(l={l}, m={mo}) max error {emax:.3f}% > {TOL_GLOB_MAX_PCT}%")

        if MPI.COMM_WORLD.Get_rank() == 0:
            fig, ax = plt.subplots(figsize=(6, 3.8))
            lg = np.arange(2, 60)
            ax.plot(lg, [resp(l) * kf for l in lg], "k-", lw=1,
                    label="exact: drho g / (drho g + D P_l)")
            ax.plot(lg, [kf / (kf + D * (l * (l + 1) / R ** 2) ** 2) for l in lg],
                    "k:", lw=1, label="flat-plate limit")
            ax.plot(ls, np.array(ratio) * kf, "o", label="goSPL (mesh pipeline)")
            ax.set_xlabel("spherical-harmonic degree l")
            ax.set_ylabel("isostatic response  w / w_Airy")
            ax.set_title(f"Global flexure, Te={TE_GLOB / 1e3:.0f} km")
            ax.legend()
            _save(fig, tmp_path, request, "flexure_global_response.png")
    finally:
        m.destroy()
