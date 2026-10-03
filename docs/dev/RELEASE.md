# Versioning and release

The procedure is scripted in the `release` skill (`.claude/skills/release/SKILL.md`);
this file holds the background. Moved out of `AGENTS.md` (2026-10).

## `__version__` and the version literals

**`__version__` is set in `gospl/__init__.py`** via `importlib.metadata` reading the installed package metadata:
```python
try:
    from importlib.metadata import version, PackageNotFoundError
    __version__ = version("gospl")
except PackageNotFoundError:
    __version__ = "unknown"
```
The metadata version is driven by `meson.build` line 4 (`version: '2026.7.14'` as of this review). `gospl.__version__` derives from it via `importlib.metadata` — never hardcode the version in `__init__.py`. The `PackageNotFoundError` fallback covers the case where the package is cloned but not installed (e.g. bare `git clone` without `pip install -e .`).

**There is one other version literal that MUST be kept in sync** with `meson.build` at every bump: `conda/meta.yaml` line 2 (`{% set version = "..." %}`). conda-build does NOT introspect `meson.build`; it resolves the jinja `version` literally and embeds it in the `.conda` artefact filename. A drift between the two means the conda channel publishes a different version number than the installed `gospl.__version__`, and the workflow's `anaconda upload --skip-existing` silently no-ops if the conda value matches an already-published release. This bit us once between `v2026.6.12` (PyPI+Docker only) and `v2026.6.13`; the inline comment in `conda/meta.yaml` flags it for future contributors. **Bump checklist: change `meson.build:4` AND `conda/meta.yaml:2` together, always.**

**Version spelling convention (adopted 2026-06-12): no leading zeros on month or day** — e.g. `2026.6.13`, not `2026.06.13`. PyPI auto-normalizes per PEP 440 (strips leading zeros for display and in the wheel/sdist filename), so a `2026.06.13` `meson.build` would show up on PyPI as `2026.6.13` while conda artifacts retained the `2026.06.13` spelling — the two channels would visually diverge for the same release. Writing the no-zero form everywhere keeps PyPI display, conda display, `.conda` filename, `.tar.gz` sdist filename, git tag (`v2026.6.13`), and `gospl.__version__` all bitwise-identical. Past tags (`v2026.06.08`, `v2026.06.11`) stay as historical record; do NOT retroactively re-spell them.

**Since 2026-07 there are FOUR version locations to bump together:** `meson.build:4`,
`conda/meta.yaml:2`, `docs/conf.py` (`version = "..."`; `release` mirrors it), and
`docs/_static/version_switch.json` (stable + a new entry). `docs/conf.py` is the one
that gets forgotten (the 2026.6.24 bump missed it, so Read the Docs rendered a stale
version for `v2026.6.30`).
