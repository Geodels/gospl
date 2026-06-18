"""
Pytest configuration and shared fixtures for goSPL regression tests.

This file registers the `slow` marker and provides two model fixtures used
by the full-stack regression tests. Both fixtures pytest.skip() when their
backing YAML + mesh files are not present, so the suite can still run the
fast (parser-only) tests in environments without committed example data.

NEEDS_HUMAN_REVIEW
------------------
The fixtures below look for `tests/fixtures/minimal.yml` and
`tests/fixtures/incising.yml` plus the mesh npz files they reference. None
of those are committed today. Before this suite is useful in CI, the team
should commit a tiny example pair under `tests/fixtures/`:

  minimal.yml  -- small flat mesh (~100 nodes), short run (dt=10 yr,
                  end=20 yr), stratNb=0, iceOn=false, flexOn=false,
                  spl.G=0, domain.nodep=false (so sedChange runs).
                  Used by test_rcvIDi_is_copy_not_reference and
                  test_scratch_vec_trap.

  incising.yml -- tilted plane mesh with rain, spl.G=0, domain.nodep=true
                  (so flowAccumulation + erodepSPL run but sedChange does
                  not). The gradient must be steep enough that at least
                  one node sees |dz| > 1 µm over a single timestep.
                  Used by test_Eb_sign_convention.

Both should be small enough to run in < 10 s on a single rank.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def pytest_configure(config):
    """Register the `slow` and `benchmark` markers so the pytest
    `-m 'not slow'` and `-m 'not benchmark'` selectors work."""
    config.addinivalue_line(
        "markers",
        "slow: requires a full goSPL Model instantiation "
        "(mesh + PETSc DMPlex + MPI). Skipped when fixture inputs are absent.",
    )
    config.addinivalue_line(
        "markers",
        "benchmark: analytical benchmark against exact solutions — "
        "requires scipy, matplotlib, pandas; skipped automatically "
        "if dependencies unavailable.",
    )


def _instantiate(yml_name: str):
    """
    Try to load `tests/fixtures/<yml_name>` into a goSPL Model.

    Skips the calling test with a descriptive message when:
      - the YAML or its mesh are not committed,
      - goSPL itself cannot import (missing petsc4py / Fortran extension),
      - the model raises during init (so the test does not silently pass).

    Tests that need a model should request the `minimal_model` or
    `incising_model` fixture rather than calling this helper directly.
    """
    yml = FIXTURES_DIR / yml_name
    if not yml.exists():
        pytest.skip(
            f"{yml.relative_to(yml.parent.parent.parent)} not present. "
            f"NEEDS_HUMAN_REVIEW: commit a tiny mesh npz + YAML pair under "
            f"tests/fixtures/ so the regression suite can instantiate "
            f"Model. See tests/conftest.py header for the required shape."
        )
    try:
        from gospl.model import Model
    except Exception as exc:  # pragma: no cover - environment-dependent
        pytest.skip(f"Cannot import gospl.model: {exc!r}")

    # cd into the fixture dir so any relative `npdata:` path inside the
    # YAML resolves the same way it does for a user invocation.
    cwd = os.getcwd()
    os.chdir(yml.parent)
    try:
        model = Model(str(yml.name), verbose=False, showlog=False)
    finally:
        os.chdir(cwd)
    return model


@pytest.fixture
def minimal_model():
    """
    Smallest workable Model: short run, no stratigraphy, no ice/flex,
    deposition ON (so `sedChange` is called from `runProcesses`).

    See tests/conftest.py header for what the backing fixture file must
    look like.
    """
    return _instantiate("minimal.yml")


@pytest.fixture
def minimal_picard_model():
    """
    Same as minimal_model but with the opt-in lagged-diffusivity (Picard)
    marine/lake diffusion solver (`diffusion: marineSolver: picard`). Used to
    test the opt-in parsing and that the Picard deposit matches the default TS
    on the minimal fixture. See minimal_picard.yml.
    """
    return _instantiate("minimal_picard.yml")


@pytest.fixture
def minimal_strat_model():
    """
    Minimal model with stratigraphy recording ON (single-fraction). Baseline
    for the dual-lithology all-coarse parity test. See minimal_strat.yml.
    """
    return _instantiate("minimal_strat.yml")


@pytest.fixture
def minimal_dual_model():
    """
    Minimal model with dual-lithology (coarse/fine) stratigraphy enabled
    (realistic fine contrast). Exercises the full dual path end-to-end:
    erodeStrat split, _getSedFlux fine routing, deposeStrat fineFrac,
    per-fraction compaction, fine-pile advection. See minimal_dual.yml.
    """
    return _instantiate("minimal_dual.yml")


@pytest.fixture
def minimal_dual_coarse_model():
    """
    Dual lithology enabled but configured all-coarse (bedrock_coarse_frac=1.0,
    no contrast). Must reproduce minimal_strat_model. See minimal_dual_coarse.yml.
    """
    return _instantiate("minimal_dual_coarse.yml")


@pytest.fixture
def minimal_ice_model():
    """
    Minimal model with the diagnostic glacial model (opt-in via the `ice`
    section). Exercises the routing-based ice diagnostic. See minimal_ice.yml
    and docs/DESIGN_ICE_SHEET.md.
    """
    return _instantiate("minimal_ice.yml")


@pytest.fixture
def minimal_ice_till_model():
    """
    Minimal glacial model with abrasion + till on (`ice.abrasion.Kg`,
    `ice.till.on`). Exercises the till production/transport/deposition path.
    See minimal_ice_till.yml.
    """
    return _instantiate("minimal_ice_till.yml")


@pytest.fixture
def minimal_ice_flex_model():
    """
    Minimal glacial model with global flexural isostasy on. Confirms the ice
    thickness feeds the existing ice-loading path in applyFlexure. See
    minimal_ice_flex.yml.
    """
    return _instantiate("minimal_ice_flex.yml")


@pytest.fixture
def flat_fem_flex_model():
    """
    Flat (2D) model with the parallel FV biharmonic flexure solver
    (`flexure: method: fem`) — solved directly on the DMPlex (no gFlex/regular
    grid). 16 km domain, thin Te so the flexural response is well resolved and
    decays inside the domain. See flatbig_fem.yml.
    """
    return _instantiate("flatbig_fem.yml")


@pytest.fixture
def flat_wall_model():
    """
    Flat (2D) model with ALL edges set to wall/closed (`bc: 'wwww'`) and the sea
    far below the domain — a fully closed continental box. Used to verify that
    wall boundaries do not leak flow/sediment (mass must be conserved, as on a
    closed sphere). See flat_wall.yml.
    """
    return _instantiate("flat_wall.yml")


@pytest.fixture
def cyclic_cyl_model():
    """
    Flat model with a CYCLIC (periodic) boundary pair (`bc: '0c0c'` → E/W
    cyclic, N/S open). The mesh is a cylinder (the periodic axis mapped to a
    circle), which is intrinsically flat, so its FV geometry equals a periodic
    strip's and the seam cells link the two edges in the neighbour graph. See
    cyclic_cyl.yml / cyclic_cyl.npz.
    """
    return _instantiate("cyclic_cyl.yml")


@pytest.fixture
def cyclic_advect_model():
    """
    Cyclic (periodic E/W) cylinder model with horizontal advection: a localised
    elevation bump just inside the seam and a uniform around-seam displacement
    (flat vx > 0) that carries it across the periodic boundary. Exercises the
    flat→cylinder velocity transform and the cyclic-seam handling of the
    advection solver. See cyclic_cyl_advect.yml / cyclic_cyl_advect.npz.

    Keeps the fixtures dir as the working directory for the whole test: the
    tectonics `hdisp` displacement npz is loaded lazily during runProcesses(),
    so its relative path must still resolve after model construction.
    """
    yml = FIXTURES_DIR / "cyclic_cyl_advect.yml"
    if not yml.exists():
        pytest.skip(f"{yml.name} not present under tests/fixtures/")
    try:
        from gospl.model import Model
    except Exception as exc:  # pragma: no cover - environment-dependent
        pytest.skip(f"Cannot import gospl.model: {exc!r}")
    cwd = os.getcwd()
    os.chdir(yml.parent)
    try:
        model = Model(str(yml.name), verbose=False, showlog=False)
        yield model
    finally:
        os.chdir(cwd)


@pytest.fixture
def minimal_ice_dual_model():
    """
    Minimal glacial model with abrasion + till on AND dual-lithology
    stratigraphy. Exercises the strata coupling of glacial till
    (iceplex._glacialTillStrata): abraded rock removed from the pile and
    re-deposited as a moraine layer, split coarse/fine. See minimal_ice_dual.yml.
    """
    return _instantiate("minimal_ice_dual.yml")


@pytest.fixture
def minimal_ice_soil_model():
    """
    Minimal model combining the diagnostic ('mfd') glacial driver with the
    soil-aware non-linear SPL (`soil:` block, cptSoil) + glacial abrasion + till.
    Confirms glacial erosion/till coexists with soil production. See
    minimal_ice_soil.yml.
    """
    return _instantiate("minimal_ice_soil.yml")


@pytest.fixture
def minimal_prov_model():
    """
    Minimal model with in-model provenance tracers on (`provenance:`, Phase 0).
    Exercises stratP allocation/seeding and the opt-in gating. See
    minimal_prov.yml and docs/DESIGN_PROVENANCE.md §6.
    """
    return _instantiate("minimal_prov.yml")


@pytest.fixture
def minimal_prov_multi_model():
    """
    Minimal provenance model with TWO source classes from a per-vertex map.
    The source map (`prov_src.npz`) is generated here from the mesh coordinates
    (class 0 / 1 split by x) so no binary fixture is committed. See
    minimal_prov2.yml.
    """
    import numpy as np

    mesh = FIXTURES_DIR / "mesh.npz"
    if not mesh.exists():
        pytest.skip("mesh.npz fixture not present")
    v = np.load(mesh)["v"]
    rock = (v[:, 0] > v[:, 0].mean()).astype(np.int64)   # 2 classes, split by x
    np.savez(FIXTURES_DIR / "prov_src.npz", rock=rock)
    return _instantiate("minimal_prov2.yml")


@pytest.fixture
def incising_model():
    """
    Tilted-plane Model with rain and `nodep: true` so the river-incision
    branch runs without marine deposition rearranging the sign budget.

    See tests/conftest.py header for what the backing fixture file must
    look like.
    """
    return _instantiate("incising.yml")


@pytest.fixture
def minimal_model_with_evap():
    """
    Minimal model with uniform evaporation injected post-construction
    (= half of the YAML uniform rain rate). Skips if the backing YAML
    does not declare uniform rainfall (the test it serves cannot reason
    about the FA reduction otherwise).

    See DESIGN_EVAPORATION.md §4 T1.
    """
    import pandas as pd

    model = _instantiate("minimal.yml")
    if model.raindata is None:
        pytest.skip("minimal.yml has no rainfall; evap-reduces-FA test needs rain")
    rUni = model.raindata.at[0, "rUni"]
    if pd.isnull(rUni):
        pytest.skip("minimal.yml rain is not uniform; evap test needs scalar rain")
    model.evapdata = pd.DataFrame(
        [{"start": 0.0, "eUni": 0.5 * float(rUni), "eMap": None, "eKey": None}]
    )
    # Reset evap counters / per-node fields so the first applyForces call
    # in runProcesses populates them from the injected evapdata.
    model.evapNb = -1
    model.evapVal = None
    model.evapMesh = None
    model.evapLoss = 0.0
    return model
