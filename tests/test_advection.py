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


def test_scheme_advection_leaves_stratigraphy_known_bug(tmp_path):
    """
    KNOWN BUG, pinned at its CURRENT behaviour (AGENTS.md > Known bugs):
    the advection schemes (`advect: upwind|iioe1|iioe2` + `hdisp`) move the
    surface fields in `tectonics._varAdvector` (elevation, cumED, flexure,
    soil) but never the stratigraphic record, so the layers stay where they
    were deposited while the surface above them moves. (Plate interpolation,
    `advect: interp`, does move them: `test_plate_advection_moves_stratigraphy`.)

    When fixed this test FAILS on purpose: replace it with a guard that the
    layers move with the surface.
    """
    import shutil

    from gospl.model import Model

    for f in ("flat_advect.yml", "flat_advect.npz"):
        shutil.copy(FIXTURES_DIR / f, tmp_path / f)
    yml = (tmp_path / "flat_advect.yml").read_text()
    yml = yml.replace("    dt: 2.5e3\n", "    dt: 2.5e3\n    strat: 2.5e3\n")
    assert "strat: 2.5e3" in yml
    (tmp_path / "adv_strat.yml").write_text(yml)

    cwd = os.getcwd()
    os.chdir(tmp_path)
    model = None
    seen = []
    try:
        model = Model("adv_strat.yml", verbose=False)
        orig = model._varAdvector

        def pinned():
            z0 = model.stratZ.copy()
            h0 = model.hLocal.getArray().copy()
            orig()
            if model.stratStep > 0:
                moved = np.abs(model.hLocal.getArray() - h0).max()
                seen.append((moved, np.array_equal(model.stratZ, z0)))

        model._varAdvector = pinned
        model.runProcesses()
    finally:
        if model is not None:
            model.destroy()
        os.chdir(cwd)
    assert seen, "no advection step with recorded layers"
    assert all(moved > 1.0 for moved, _ in seen), "the surface did not move"
    assert all(same for _, same in seen), "the layers moved: the known bug may be fixed"
