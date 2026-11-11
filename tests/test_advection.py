"""
Horizontal advection schemes (serial).

Protects: AGENTS.md > Fixed (IIOE2 zeroed a no-overshoot field).

Split out of the former tests/test_regression.py. Run this group alone with
`pytest -m advection`; see tests/README.md for the marker list.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from _helpers import FIXTURES_DIR

# Skip the whole module (rather than erroring at collection) when the
# goSPL runtime stack is not installed.
pytest.importorskip(
    "gospl.tools.inputparser",
    reason="goSPL runtime deps (petsc4py / ruamel.yaml / scipy) not installed",
)

pytestmark = [pytest.mark.advection]


def test_advection_iioe2_smooth_field_not_zeroed(tmp_path):
    """
    Protects: IIOE2 (advect='iioe2', advscheme=3) must NOT zero a field when the
    Scheme-1 result has no overshoot. `tectonics._advectorIIOE2` summed the
    overshoot into `self.tmp` — which still held the Scheme-1 advected field, the
    result to keep when there is no overshoot. When the global excess was 0 the
    anti-diffusion branch was skipped, leaving the all-zero `diff` in `tmp`, so
    `_advectField` copied zeros back into the field: the elevation flattened to
    0 on the first advection step (reported on the stratigraphic_record
    `input-advect.yml` example). `test_advection_parallel` advects a sharp RIDGE,
    which always overshoots (excess>0) and so never hit this path. A spatially
    **uniform** field is used here: its neighbourhood min==max, so the Scheme-1
    advection produces zero overshoot (excess==0) deterministically — advection
    of a uniform field must return that field, but the bug zeroed it.
    """
    import shutil

    fixtures_dir = FIXTURES_DIR
    yml_src = fixtures_dir / "flat_advect.yml"
    npz_src = fixtures_dir / "flat_advect.npz"
    if not (yml_src.exists() and npz_src.exists()):
        pytest.skip("flat_advect fixtures not present.")
    try:
        from gospl.model import Model
    except Exception as exc:  # pragma: no cover - environment-dependent
        pytest.skip(f"Cannot import gospl.model: {exc!r}")

    shutil.copy(yml_src, tmp_path / "flat_advect.yml")
    shutil.copy(npz_src, tmp_path / "flat_advect.npz")
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        model = Model("flat_advect.yml", verbose=False, showlog=False)
        assert model.advscheme == 3, "fixture must use the iioe2 scheme"

        # Uniform field: neighbourhood min == max everywhere, so the IIOE2
        # anti-diffusion correction finds zero overshoot (excess == 0) — the
        # exact no-overshoot path the bug zeroed. Advection must preserve a
        # uniform field (SPL K=0, no diffusion/flexure block, so the step is
        # pure advection of a flat field above the -1000 m sea level).
        const = 50.0
        model.hLocal.setArray(np.full(model.lpoints, const))
        model.dm.localToGlobal(model.hLocal, model.hGlobal)
        model.dm.globalToLocal(model.hGlobal, model.hLocal)

        model.runProcesses()
        z1 = model.hLocal.getArray().copy()
        model.destroy()
    finally:
        os.chdir(cwd)

    assert np.isfinite(z1).all(), "advected elevation is not finite"
    # The uniform field must survive — the bug collapsed it to 0 everywhere.
    assert z1.min() > 25.0 and z1.max() > 25.0, (
        "iioe2 zeroed a uniform field (excess==0 bug): elev range "
        f"[{z1.min():.2f}, {z1.max():.2f}], expected ~{const:.0f}"
    )
    np.testing.assert_allclose(z1.mean(), const, atol=1.0)


def test_plate_advection_moves_stratigraphy(tmp_path):
    """
    Protects: plate-interpolation advection (`advect: interp` + `hdisp`) with
    stratigraphy on (`tectonics._advectStrati` -> `STRAMesh.stratalRecord`).

    Silent failure prevented: the path crashed at the first plate update
    (`AttributeError: 'Model' object has no attribute 'idNbghs'`, unassigned
    since the July 2024 refactor, 5ea92c3) and, past that, indexed the 1-D
    `tec_sumw` as 3-D. It also moved only stratH/stratZ/phiS, leaving stratK
    (and the dual/provenance/duricrust layers) behind (fixed 2026-11-11;
    docs/dev/FIXED_BUGS.md).

    Invariant: under a shear plate motion, right after the plate update every
    per-layer field equals an independent numpy inverse-distance
    interpolation of its pre-update values with the same `tec_*` tables
    (to 1e-12 relative; exact on macOS), and the run completes.
    """
    import shutil

    from gospl.model import Model

    for f in ("minimal_strat.yml", "mesh.npz"):
        shutil.copy(FIXTURES_DIR / f, tmp_path / f)
    with np.load(tmp_path / "mesh.npz") as d:
        v = d["v"]
    vel = np.zeros_like(v)
    vel[:, 0] = (np.ptp(v[:, 0]) / 600.0) * (v[:, 1] - v[:, 1].mean()) / np.ptp(v[:, 1])
    np.savez(tmp_path / "shear.npz", xyz=vel)
    yml = (tmp_path / "minimal_strat.yml").read_text()
    yml = yml.replace("domain:\n", "domain:\n    advect: 'interp'\n", 1)
    yml += ("\ntectonics:\n  - start: 0.\n    end: 60.\n"
            "    hdisp: ['shear', 'xyz']\n  - start: 60.\n    end: 100.\n")
    (tmp_path / "plate_strat.yml").write_text(yml)

    cwd = os.getcwd()
    os.chdir(tmp_path)
    model = None
    checked = []
    try:
        model = Model("plate_strat.yml", verbose=False)
        orig = model._advectStrati

        def checked_advect():
            before = {k: a.copy() for k, a in model._stratFields().items()}
            ids, w, on = model.tec_IDs, model.tec_weights, model.tec_onIDs
            orig()
            after = model._stratFields()
            for k, a in before.items():
                exp = np.einsum("np,np...->n...", w, a[ids]) / w.sum(1)[:, None]
                if len(on):
                    exp[on] = a[ids[on, 0]]
                np.testing.assert_allclose(after[k], exp, rtol=1e-12, atol=1e-9,
                                           err_msg=f"{k} not advected as expected")
                checked.append(k)

        model._advectStrati = checked_advect
        model.runProcesses()
        assert model.tNow == model.tEnd
    finally:
        if model is not None:
            model.destroy()
        os.chdir(cwd)
    assert {"stratH", "stratZ", "phiS", "stratK"} <= set(checked), checked


def test_tectonic_event_without_maps_runs(tmp_path):
    """
    Protects: `Tectonics.getTectonics` with an event that carries no map
    (only `start`/`end`, e.g. to stop the motion), which the parser accepts
    (`_defineTectonics` fills "empty").

    Silent failure prevented: the method ended with `del mdata`, a name only
    bound when the event loads a map, so such an event crashed the run with
    `UnboundLocalError: cannot access local variable 'mdata'` (fixed
    2026-11-11; docs/dev/FIXED_BUGS.md).

    Invariant: `minimal.yml` with a map-free tectonics event runs to the end
    and leaves no horizontal/vertical forcing set.
    """
    import shutil

    from gospl.model import Model

    for f in ("minimal.yml", "mesh.npz"):
        shutil.copy(FIXTURES_DIR / f, tmp_path / f)
    yml = (tmp_path / "minimal.yml").read_text()
    yml += "\ntectonics:\n  - start: 0.\n    end: 100.\n"
    (tmp_path / "notec.yml").write_text(yml)

    cwd = os.getcwd()
    os.chdir(tmp_path)
    model = None
    try:
        model = Model("notec.yml", verbose=False)
        model.runProcesses()
        assert model.tNow == model.tEnd
        assert model.hdisp is None and model.upsub is None
    finally:
        if model is not None:
            model.destroy()
        os.chdir(cwd)


@pytest.mark.parametrize("scheme", ["upwind", "iioe1", "iioe2"])
def test_scheme_advection_moves_stratigraphy(tmp_path, scheme):
    """
    Protects: the advection schemes move the stratigraphic record
    (`tectonics._advectStratiScheme`: semi-Lagrangian, departure point
    ``x - v*dt``, then `stratalRecord`).

    Silent failure prevented: `_varAdvector` advected the surface fields only,
    so the layers stayed where they were deposited while the surface above
    them moved (fixed 2026-11-11; docs/dev/FIXED_BUGS.md).

    Invariant: on a regular grid with a velocity of exactly one cell per
    step, every departure point is a mesh node, so after each advection step
    every per-layer field at interior node (i, j) equals its pre-step value
    at (i-1, j) exactly; the upstream column (departure outside the mesh) is
    excluded.
    """
    import shutil

    from gospl.model import Model

    n, dx, dt = 30, 1000.0, 1000.0
    x, y = np.meshgrid(np.arange(n) * dx, np.arange(n) * dx)
    v = np.c_[x.ravel(), y.ravel(), np.zeros(n * n)]
    i, j = np.meshgrid(np.arange(n - 1), np.arange(n - 1))
    a = (j * n + i).ravel()
    cells = np.r_[np.c_[a, a + 1, a + n + 1], np.c_[a, a + n + 1, a + n]].astype(np.int32)
    z = 100.0 + 0.01 * v[:, 0] + 0.003 * v[:, 1] ** 1.2       # varies along x and y
    vel = np.zeros_like(v)
    vel[:, 0] = dx / dt                                       # exactly one cell per step
    np.savez(tmp_path / "grid.npz", v=v, c=cells, z=z, uv=vel)
    (tmp_path / "adv.yml").write_text(f"""name: strata advection
domain:
    npdata: ['grid','v','c','z']
    flowdir: 5
    bc: 'oooo'
    advect: '{scheme}'
time:
    start: 0.
    end: 4000.
    tout: 4000.
    dt: {dt}
    strat: {dt}
spl:
    K: 0.
    d: 0.
    m: 0.5
sea:
    position: -1000.
climate:
  - start: 0.
    uniform: 1.
tectonics:
    - start: 0.
      end: 4000.
      hdisp: ['grid','uv']
output:
    dir: 'out'
""")

    cwd = os.getcwd()
    os.chdir(tmp_path)
    model = None
    checks = []
    try:
        model = Model("adv.yml", verbose=False)
        orig = model._advectStratiScheme
        col = np.rint(model.lcoords[:, 0] / dx).astype(int)
        row = np.rint(model.lcoords[:, 1] / dx).astype(int)
        key = {(c, r): k for k, (c, r) in enumerate(zip(col, row))}
        inner = np.where(col >= 1)[0]
        src = np.array([key[(col[k] - 1, row[k])] for k in inner])

        def checked():
            before = {k: a.copy() for k, a in model._stratFields().items()}
            orig()
            after = model._stratFields()
            for k, b in before.items():
                assert np.array_equal(after[k][inner], b[src]), f"{k} not shifted by one cell"
            checks.append(model.stratStep)

        model._advectStratiScheme = checked
        model.runProcesses()
    finally:
        if model is not None:
            model.destroy()
        os.chdir(cwd)
    assert checks and max(checks) >= 2, checks


def _cyclic_strata_model(tmp_path, scheme, speed):
    """minimal cyclic-cylinder model with strata recorded every step."""
    import shutil

    shutil.copy(FIXTURES_DIR / "cyclic_cyl_advect.npz", tmp_path / "cyl.npz")
    with np.load(tmp_path / "cyl.npz") as d:
        data = {k: d[k] for k in d.files}
    data["uv"] = np.zeros_like(data["uv"])
    data["uv"][:, 0] = speed
    np.savez(tmp_path / "cyl.npz", **data)
    yml = (FIXTURES_DIR / "cyclic_cyl_advect.yml").read_text()
    yml = yml.replace("cyclic_cyl_advect','v'", "cyl','v'").replace(
        "['cyclic_cyl_advect','uv']", "['cyl','uv']")
    yml = yml.replace("advect: 'iioe2'", f"advect: '{scheme}'")
    yml = yml.replace("    dt: 2.5e3\n", "    dt: 2.5e3\n    strat: 2.5e3\n")
    assert "strat: 2.5e3" in yml and "['cyl','uv']" in yml
    (tmp_path / "cyl.yml").write_text(yml)
    return "cyl.yml"


@pytest.mark.parametrize("scheme", ["upwind", "iioe2"])
def test_cyclic_scheme_advection_moves_stratigraphy_across_seam(tmp_path, scheme):
    """
    Protects: on a cyclic (cylinder) model the semi-Lagrangian stratigraphy
    move (`_advectStratiScheme`) uses the cylinder-TANGENT node velocity
    (`_cylinderVelocity`), so layers wrap across the periodic seam instead of
    being pushed radially off the surface or stopped at the seam.

    Invariant: the cylinder fixture has 48 angular columns 250 m apart; at
    0.1 m/yr and dt=2500 yr the record moves one column per step. A marker
    set in the seam column ends up, in every row, in the next column across
    the seam (>= 99% of its value: the departure point is computed exactly on
    the cylinder, `_cylinderDisplace`), and the seam column is emptied.
    """
    from gospl.model import Model

    yml = _cyclic_strata_model(tmp_path, scheme, 0.1)
    cwd = os.getcwd()
    os.chdir(tmp_path)
    model = None
    seen = []
    try:
        model = Model(yml, verbose=False)
        assert model.cyclicBC
        th = np.mod(np.arctan2(model.lcoords[:, 2], model.lcoords[:, 0]), 2 * np.pi)
        col = np.mod(np.rint(th / (2 * np.pi / 48)).astype(int), 48)
        row = np.rint(model.lcoords[:, 1] / 300.0).astype(int)
        seam_col = np.bincount(col[model.cyclicPts]).argmax()
        orig = model._advectStratiScheme

        def checked():
            marker = (col == seam_col).astype(float)
            model.stratZ[:, : model.stratStep] = marker[:, None]
            orig()
            z = model.stratZ[:, 0]
            nxt = (seam_col + 1) % 48
            for r in np.unique(row):
                in_row = row == r
                best = col[in_row][np.argmax(z[in_row])]
                assert best == nxt, f"row {r}: marker at column {best}, expected {nxt}"
                assert z[in_row & (col == nxt)].min() > 0.99
                assert z[in_row & (col == seam_col)].max() < 0.01
            seen.append(model.stratStep)

        model._advectStratiScheme = checked
        model.runProcesses()
    finally:
        if model is not None:
            model.destroy()
        os.chdir(cwd)
    assert seen, "no stratigraphy move happened"


def test_sphere_scheme_advection_moves_stratigraphy(tmp_path):
    """
    Protects: the semi-Lagrangian stratigraphy move on a global (spherical)
    mesh, where the departure point ``x - v*dt`` of a rotation lies slightly
    off the sphere; the 3-D nearest-node lookup must still pick the upstream
    nodes.

    Invariant: under a rigid rotation about the z axis (iioe2), a Gaussian
    bump of a layer field centred on the equator moves its centroid by
    omega*dt in longitude (within 15%), and a uniform field stays uniform.
    """
    from scipy.spatial import ConvexHull

    from gospl.model import Model

    n, R = 6000, 6378137.0
    k = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * k / n)
    th = np.pi * (1 + 5 ** 0.5) * k
    v = R * np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]
    cells = ConvexHull(v).simplices.astype(np.int32)
    dt = 1.0e4
    spacing = np.sqrt(4 * np.pi / n) * R               # ~ mean node spacing
    omega = spacing / R / dt                           # ~ one node spacing per step
    vel = np.c_[-omega * v[:, 1], omega * v[:, 0], np.zeros(n)]
    np.savez(tmp_path / "sph.npz", v=v, c=cells, z=np.full(n, 500.0), uv=vel)
    (tmp_path / "sph.yml").write_text(f"""name: sphere strata advection
domain:
    npdata: ['sph','v','c','z']
    flowdir: 5
    advect: 'iioe2'
time:
    start: 0.
    end: {3 * dt}
    tout: {3 * dt}
    dt: {dt}
    strat: {dt}
spl:
    K: 0.
    d: 0.
    m: 0.5
sea:
    position: -1000.
climate:
  - start: 0.
    uniform: 1.
tectonics:
    - start: 0.
      end: {3 * dt}
      hdisp: ['sph','uv']
output:
    dir: 'out'
""")
    cwd = os.getcwd()
    os.chdir(tmp_path)
    model = None
    shifts = []
    try:
        model = Model("sph.yml", verbose=False)
        lc = model.lcoords
        lon = np.arctan2(lc[:, 1], lc[:, 0])
        lat = np.arcsin(np.clip(lc[:, 2] / R, -1, 1))
        width = 6 * spacing / R
        orig = model._advectStratiScheme

        def checked():
            bump = np.exp(-((lon / width) ** 2 + (lat / width) ** 2))
            model.stratZ[:, : model.stratStep] = bump[:, None]
            model.phiS[:, : model.stratStep] = 0.3
            orig()
            z = model.stratZ[:, 0]
            near = np.abs(lon) < 0.5
            shifts.append(np.sum(z[near] * lon[near]) / np.sum(z[near])
                          - np.sum(bump[near] * lon[near]) / np.sum(bump[near]))
            assert np.allclose(model.phiS[:, : model.stratStep], 0.3, rtol=0, atol=1e-12)

        model._advectStratiScheme = checked
        model.runProcesses()
    finally:
        if model is not None:
            model.destroy()
        os.chdir(cwd)
    assert shifts, "no stratigraphy move happened"
    for s in shifts:
        assert abs(s - omega * dt) < 0.15 * omega * dt, (s, omega * dt)


def test_cyclic_plate_advection_wraps_across_seam(tmp_path):
    """
    Protects: plate-interpolation advection (`advect: interp`) on a cyclic
    (cylinder) model displaces the nodes AROUND the cylinder
    (`_cylinderDisplace`) when building the interpolation tables.

    Silent failure prevented: `_readAdvectionData` displaced the 3-D cylinder
    nodes along the flat-frame `hdisp` (a straight x shift), so the surface
    and the stratigraphy were smeared in place instead of moving around the
    seam (fixed 2026-11-11; docs/dev/FIXED_BUGS.md).

    Invariant: one plate stage moving 250 m (one of the 48 columns) carries a
    marker in the seam column, in BOTH the elevation and a layer field, to the
    next column across the seam in every row (>= 99%).
    """
    from gospl.model import Model

    yml = _cyclic_strata_model(tmp_path, "interp", 0.05)      # 250 m over 5 kyr
    text = (tmp_path / yml).read_text()
    old = "    - start: 0.\n      end: 1.e4\n      hdisp: ['cyl','uv']"
    assert old in text
    (tmp_path / yml).write_text(text.replace(
        old, "    - start: 0.\n      end: 5.e3\n      hdisp: ['cyl','uv']\n"
             "    - start: 5.e3\n      end: 1.e4"))
    cwd = os.getcwd()
    os.chdir(tmp_path)
    model = None
    seen = []
    try:
        model = Model(yml, verbose=False)
        th = np.mod(np.arctan2(model.lcoords[:, 2], model.lcoords[:, 0]), 2 * np.pi)
        col = np.mod(np.rint(th / (2 * np.pi / 48)).astype(int), 48)
        row = np.rint(model.lcoords[:, 1] / 300.0).astype(int)
        seam_col = np.bincount(col[model.cyclicPts]).argmax()
        nxt = (seam_col + 1) % 48
        orig = model._advectPlates

        def checked():
            marker = (col == seam_col).astype(float)
            model.hLocal.setArray(100.0 + 50.0 * marker)
            model.dm.localToGlobal(model.hLocal, model.hGlobal)
            model.stratZ[:, : model.stratStep] = marker[:, None]
            orig()
            h = (model.hLocal.getArray() - 100.0) / 50.0
            z = model.stratZ[:, 0]
            for fld, name in ((h, "elevation"), (z, "stratZ")):
                for r in np.unique(row):
                    in_row = row == r
                    best = col[in_row][np.argmax(fld[in_row])]
                    assert best == nxt, f"{name} row {r}: marker at column {best}, expected {nxt}"
                    assert fld[in_row & (col == nxt)].min() > 0.99, name
            seen.append(True)

        model._advectPlates = checked
        model.runProcesses()
    finally:
        if model is not None:
            model.destroy()
        os.chdir(cwd)
    assert seen, "no plate update happened"
