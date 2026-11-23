"""
Tests for the ``gospl`` command-line entry point (:mod:`gospl.cli`).
"""

import os

import pytest

from _helpers import FIXTURES_DIR

pytestmark = pytest.mark.cli

FIXTURES = FIXTURES_DIR


def test_cli_parse_aliases():
    """`--i`/`--v` (explicit aliases) and `-i`/`-v` map to input/verbose."""
    from gospl.cli import _parse_args

    a = _parse_args(["--i", "in.yml", "--v"])
    assert a.input == "in.yml" and a.verbose is True
    # Short forms + defaults (verbose/log/profile off unless requested).
    b = _parse_args(["-i", "in.yml"])
    assert b.input == "in.yml"
    assert b.verbose is False and b.log is False and b.profile is False
    # `--v` is unambiguous despite `--version` also existing.
    c = _parse_args(["--input", "in.yml", "--verbose"])
    assert c.verbose is True


def test_cli_version_flag(capsys):
    """`--version` prints the goSPL version and exits 0."""
    from gospl.cli import _parse_args
    import gospl

    with pytest.raises(SystemExit) as exc:
        _parse_args(["--version"])
    assert exc.value.code == 0
    assert gospl.__version__ in capsys.readouterr().out


def test_cli_requires_input():
    """The input file is mandatory — argparse exits when it is missing."""
    from gospl.cli import _parse_args

    with pytest.raises(SystemExit):
        _parse_args([])


@pytest.mark.slow
def test_cli_runs_minimal_model():
    """`gospl -i minimal.yml` reads, runs and tears down a model (rc 0)."""
    yml = FIXTURES / "minimal.yml"
    if not yml.exists():
        pytest.skip("minimal.yml fixture / mesh absent")
    from gospl.cli import main

    cwd = os.getcwd()
    os.chdir(FIXTURES)               # fixture YAML uses run-relative paths
    try:
        rc = main(["-i", "minimal.yml"])
    finally:
        os.chdir(cwd)
    assert rc == 0


def test_output_dir_absolute_path(tmp_path):
    """An absolute ``output: dir:`` is used as given. It used to be appended
    to the working directory (``cwd + "/" + dir``), so ``/data/run`` ended up in
    ``<cwd>/data/run`` and an absolute path was never honoured."""
    import re
    import shutil

    from gospl.model import Model

    run = tmp_path / "run"
    run.mkdir()
    shutil.copy(FIXTURES / "mesh.npz", run / "mesh.npz")
    out = tmp_path / "elsewhere" / "abs_out"
    text = (FIXTURES / "minimal.yml").read_text()
    text = re.sub(r"dir: '[^']*'", "dir: '%s'" % out, text)
    (run / "abs.yml").write_text(text)
    cwd = os.getcwd()
    os.chdir(run)
    try:
        m = Model("abs.yml", verbose=False, showlog=False)
        try:
            assert os.path.realpath(m.outputDir) == os.path.realpath(str(out))
        finally:
            m.destroy()
    finally:
        os.chdir(cwd)
    assert out.is_dir()
    assert not (run / str(out).lstrip("/")).exists()
