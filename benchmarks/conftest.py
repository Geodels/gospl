"""
Pytest configuration and fixtures for the goSPL analytical benchmark suite.

Each fixture copies the benchmark's input files into a per-test temp
directory and chdirs there. This isolates concurrent test runs, prevents
goSPL output (sim_output/, sims_outputs/, results/, *.pdf, *.md) from
landing inside the source tree, and lets pytest's tmp_path cleanup
handle teardown automatically.

Each fixture also runs a teardown step that copies the benchmark's
report artefacts (PDF figures, Markdown reports) into
`<repo_root>/results/<benchmark_name>/` so they survive pytest's
tmp_path cleanup AND are picked up by the GitHub Actions
upload-artifact@v4 step in tests-slow.yml (which uploads `results/`).

See AGENTS.md > Analytical benchmark suite.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

BENCHMARKS_DIR = Path(__file__).parent


def _copy_artefacts(src_dir: Path, dst_dir: Path, patterns: tuple[str, ...]) -> None:
    """
    Copy files matching the given glob patterns from src_dir into
    dst_dir. Used by every benchmark fixture's teardown to surface
    report files out of pytest's per-test tmp dir.

    Skips silently when a file doesn't exist (e.g. the test failed
    before generating it) so teardown never masks the underlying error.
    """
    dst_dir.mkdir(parents=True, exist_ok=True)
    for pattern in patterns:
        for path in src_dir.glob(pattern):
            if path.is_file():
                shutil.copy2(path, dst_dir / path.name)


@pytest.fixture(scope="function")
def spl_tmp_path(tmp_path, monkeypatch, request):
    """
    Copy benchmarks/spl/ into tmp_path and chdir there.

    Teardown copies `spl_benchmark.pdf` and `spl_benchmark.md` from
    the tmp dir to `<rootpath>/results/spl/`.
    """
    src = BENCHMARKS_DIR / "spl"
    dest = tmp_path / "spl"
    shutil.copytree(src, dest)
    monkeypatch.chdir(dest)
    yield dest

    artefacts = Path(request.config.rootpath) / "results" / "spl"
    _copy_artefacts(dest, artefacts, ("spl_benchmark*.pdf", "spl_benchmark*.md"))


@pytest.fixture(scope="function")
def hillslope_tmp_path(tmp_path, monkeypatch, request):
    """
    Copy benchmarks/hillslope/ into tmp_path, chdir there, and prepend
    the copied directory to sys.path so `from scripts import analysis`
    resolves to the COPIED scripts/ folder (not the source tree's).
    `monkeypatch.syspath_prepend` auto-restores sys.path on teardown
    so it does not leak into other tests.

    Teardown copies everything under the test's `results/` subfolder
    (Markdown report + per-run PNGs + convergence plots) to
    `<rootpath>/results/hillslope/`.
    """
    src = BENCHMARKS_DIR / "hillslope"
    dest = tmp_path / "hillslope"
    shutil.copytree(src, dest)
    monkeypatch.chdir(dest)
    monkeypatch.syspath_prepend(str(dest))
    yield dest

    artefacts = Path(request.config.rootpath) / "results" / "hillslope"
    _copy_artefacts(dest / "results", artefacts, ("*.pdf", "*.md", "*.png"))


@pytest.fixture(scope="function")
def knickpoint_tmp_path(tmp_path, monkeypatch, request):
    """
    Copy benchmarks/knickpoint/ into tmp_path and chdir there.

    Skips with a clear message if `boundary_condition/gospl_mesh.npz` is
    absent (the committed mesh is the source of truth; if a future
    contributor removes it, the test must skip rather than fail).

    `drop_baselevel.npz` is written by the test itself during phase 2
    of the knickpoint workflow; the committed copy is overwritten in
    the tmp directory and discarded on teardown.

    Teardown copies the `results/` subfolder (PDF + Markdown) to
    `<rootpath>/results/knickpoint/`.
    """
    src = BENCHMARKS_DIR / "knickpoint"
    mesh = src / "boundary_condition" / "gospl_mesh.npz"
    if not mesh.exists():
        pytest.skip(
            "benchmarks/knickpoint/boundary_condition/gospl_mesh.npz "
            "not found. See benchmarks/knickpoint/boundary_condition/"
            "README.md to (re)generate it."
        )
    dest = tmp_path / "knickpoint"
    shutil.copytree(src, dest)
    monkeypatch.chdir(dest)
    yield dest

    artefacts = Path(request.config.rootpath) / "results" / "knickpoint"
    _copy_artefacts(dest / "results", artefacts, ("*.pdf", "*.md"))
