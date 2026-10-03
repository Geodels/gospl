---
name: release
description: Cut a goSPL release - bump the version in all four places, merge dev into master, tag, watch the gated publish to PyPI/conda/Docker, rebuild Read the Docs. Use when asked to release, tag, publish or bump the version.
---

# Cut a goSPL release

Background and history: `docs/dev/RELEASE.md`, `conda/AGENTS.md`, `docker/AGENTS.md`.
A release is **merge `dev` → `master` (fast-forward) + push a `v*` tag**. The tag
fires `conda-build.yml`, `pypi-publish.yml` and `docker-build.yml`; each `needs:` the
`_release-gate.yml` test job, so nothing publishes unless the suite passes. **A PyPI
version string can never be reused once published**, so stop and ask before every
outward-facing step (push, tag).

## 1. Pick the version
CalVer `YYYY.M.D`, **no leading zeros** (`2026.10.3`, not `2026.10.03`). Release
candidates: `2026.10.3rc1` (PEP 440 canonical; a bare `rc` normalises differently on
PyPI and conda).

## 2. Pre-flight on `dev`
```bash
git checkout dev && git pull
python scripts/lint_mpi_collectives.py
pytest tests/ -n 4                      # full suite
pytest benchmarks/ -m benchmark         # analytical benchmarks (the gate runs these)
```
Do not tag a red tree: the gate will block the publish anyway, but a failed tag
wastes a ~45-minute CI cycle.

## 3. Bump the version in ALL FOUR places, in ONE commit
| File | What |
|---|---|
| `meson.build` line 4 | `version: '<v>'` (drives `gospl.__version__` via importlib.metadata) |
| `conda/meta.yaml` line 2 | `{% set version = "<v>" %}` (conda-build does NOT read meson.build) |
| `docs/conf.py` | `version = "<v>"` (`release` mirrors it). The one that gets forgotten |
| `docs/_static/version_switch.json` | add the new entry, point `stable` at it |

Verify: `grep -n "<v>" meson.build conda/meta.yaml docs/conf.py docs/_static/version_switch.json`
shows all four. Also add a `docs/dev/CHANGELOG_DEV.md` row for the release and update
the "Released packages" table in `conda/AGENTS.md`. Commit messages carry no Claude
co-author trailer (maintainer preference).

## 4. Merge and tag (ASK FIRST)
```bash
git checkout master && git merge --ff-only dev
git push origin master                  # ask before pushing
git tag v<v> && git push origin v<v>    # ask before pushing: this starts the publish
git checkout dev
```
`master` is normally a strict ancestor of `dev`; if `--ff-only` refuses, stop and
ask (someone committed to master directly).

## 5. Watch the gated publish
The three workflows each run the gate (`pytest tests/` on ubuntu + macOS py3.11,
benchmarks on ubuntu). `gh` lives at `/opt/anaconda3/bin/gh` and is not on PATH in a
plain shell; call it by full path, e.g.
`/opt/anaconda3/bin/gh run list --workflow pypi-publish.yml --limit 3`.
If the gate fails: nothing was published, so the version is still free. Fix on
`dev`, then move the tag: `git push origin :refs/tags/v<v>`, retag the fix commit,
push again. Do not re-tag for trivial post-tag cleanup while a publish is running.

## 6. After publish
- Read the Docs: a force-moved tag does NOT rebuild automatically. Builds → build
  `v<v>` and `stable`; if still stale, wipe the version and rebuild.
- GitHub release: create it as a **draft** for the maintainer to publish.
- `docker/slurm/*.pbs|.slurm` `CONTAINER=` paths reference the new `.sif` name.
