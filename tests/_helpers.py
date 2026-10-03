"""
Helpers shared by more than one tests/test_*.py module.

Imported as `from _helpers import ...`; tests/conftest.py puts this directory on
sys.path so the import also works under `--import-mode=importlib`.
"""

from __future__ import annotations

import atexit
import os
import shutil
import tempfile
from pathlib import Path

import pytest

_SRC_FIXTURES = Path(__file__).parent / "fixtures"


def _private_fixtures_dir() -> Path:
    """Per-process copy of the tracked fixture INPUTS (yml/npz, ~2 MB).

    Every fixture YAML writes its model output relative to its own directory
    (`output: dir:`), and `Model.__init__` rmtree's + recreates that directory
    at step 0. Running models straight out of tests/fixtures therefore (a)
    litters the source tree with output dirs and (b) races under
    `pytest -n N`: two workers instantiating `minimal.yml` delete each
    other's output (FileExistsError / FileNotFoundError). Each process -- the
    serial run, or each xdist worker -- gets its own copy, removed at exit.

    Set GOSPL_TEST_FIXTURES_INPLACE=1 to run in tests/fixtures instead (e.g.
    to inspect a test's model output afterwards).
    """
    if os.environ.get("GOSPL_TEST_FIXTURES_INPLACE"):
        return _SRC_FIXTURES
    worker = os.environ.get("PYTEST_XDIST_WORKER", "main")
    dst = Path(tempfile.mkdtemp(prefix=f"gospl-fixtures-{worker}-"))
    for f in _SRC_FIXTURES.iterdir():
        if f.is_file():                 # inputs only; subdirs are old outputs
            shutil.copy2(f, dst / f.name)
    atexit.register(shutil.rmtree, dst, ignore_errors=True)
    return dst


# The directory every test must use for fixture inputs (never the literal
# tests/fixtures path -- see _private_fixtures_dir).
FIXTURES_DIR = _private_fixtures_dir()





def _gw_model(fixture):
    """Load a groundwater fixture from tests/fixtures (skip if absent)."""
    import os
    from gospl.model import Model

    fx = str(FIXTURES_DIR)
    if not os.path.exists(os.path.join(fx, fixture)):
        pytest.skip(f"{fixture} fixture not present")
    cwd = os.getcwd()
    os.chdir(fx)
    try:
        return Model(fixture, verbose=False, showlog=False)
    finally:
        os.chdir(cwd)


def _strata_parser(stratNb):
    """
    Bare parser primed for `_extraStrata`: it reads `self.input`,
    `self.phi0s`/`self.z0s` (set by `_readCompaction` upstream), and
    `self.stratNb` (set by `_readTime` upstream). See `_bare_parser`.
    """
    # Imported lazily: conftest imports this module, and a module-level
    # importorskip here would abort the whole session when gospl is absent.
    inputparser = pytest.importorskip("gospl.tools.inputparser")
    parser = inputparser.ReadYaml.__new__(inputparser.ReadYaml)
    parser.input = {}
    parser.phi0s = 0.49
    parser.z0s = 3700.0
    parser.stratNb = stratNb
    return parser


# Wall-clock budget for one nested `mpirun` model run in tests/test_parallel.py.
# Each takes ~5-15 s locally, so a run that exceeds this has almost certainly
# DEADLOCKED (AGENTS.md > MPI contract, the #1 parallel deadlock class); failing
# fast here beats waiting out the old 600 s per test. Raise it on a very slow
# or oversubscribed machine with GOSPL_TEST_MPI_TIMEOUT=<seconds>.
MPI_TIMEOUT = float(os.environ.get("GOSPL_TEST_MPI_TIMEOUT", "240"))


def mpi_child_env():
    """Environment for a nested `mpirun` launched from inside pytest.

    pytest imported gospl (-> petsc4py.init -> MPI_Init), which under OpenMPI
    exports OMPI_*/PMIX_*/PRTE_*/OPAL_* into this process. If they leak into
    the child `mpirun` it believes it is already inside an MPI job and silently
    refuses to launch (rc=1, empty output). MPICH is unaffected, but the conda
    package and HPC container both use OpenMPI. OPAL_PREFIX is preserved so
    mpirun can still locate its own libraries.
    """
    env = {
        k: v for k, v in os.environ.items()
        if not k.startswith(("OMPI_", "PMIX_", "PRTE_", "OPAL_"))
    }
    if "OPAL_PREFIX" in os.environ:
        env["OPAL_PREFIX"] = os.environ["OPAL_PREFIX"]
    # GitHub's ubuntu runners: 4 vCPUs, 2 physical cores, one OpenMPI slot per
    # physical core; allow oversubscription so np>2 tests can launch.
    env["OMPI_MCA_rmaps_base_oversubscribe"] = "1"
    return env
