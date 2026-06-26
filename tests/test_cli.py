"""
Tests for the ``gospl`` command-line entry point (:mod:`gospl.cli`).
"""

import os
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"


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
