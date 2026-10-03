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
