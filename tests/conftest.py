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
    """Register the `slow` marker so pytest -m 'not slow' works."""
    config.addinivalue_line(
        "markers",
        "slow: requires a full goSPL Model instantiation "
        "(mesh + PETSc DMPlex + MPI). Skipped when fixture inputs are absent.",
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
