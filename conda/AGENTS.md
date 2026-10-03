# Conda package (conda/)

Read before editing `conda/meta.yaml` or `.github/workflows/conda-build.yml`, or when a conda build/install misbehaves.

This file is part of the goSPL agent guide. Read the root `AGENTS.md` first; it holds
the rules that apply everywhere (MPI contract, solver lifecycles, scratch Vecs,
conventions, commit checklist). The text below was moved verbatim from the root
file in the 2026-10 split, so dates and `## ...` cross-references still refer to
the original section names (most now live in the files listed in the root
`AGENTS.md` > Map).

## Conda Package Validation

### Released packages
| Version | Date | Channel | Install |
|---|---|---|---|
| `v2026.06.08` | 2026-06-08 | `geodels` | `mamba install -c geodels -c conda-forge gospl=2026.06.08` |
| `v2026.06.11` | 2026-06-11 | `geodels` | `mamba install -c geodels -c conda-forge gospl=2026.06.11` |
| `v2026.6.13` | 2026-06-12 | `geodels` | `mamba install -c geodels -c conda-forge gospl=2026.6.13` |
| `v2026.6.30` | 2026-06-30 | `geodels` | `mamba install -c geodels -c conda-forge gospl=2026.6.30` |
| `v2026.7.14` | 2026-07-14 | `geodels` | `mamba install -c geodels -c conda-forge gospl` |

`v2026.6.12` is intentionally absent from this table: that tag fired `pypi-publish` and `docker-build` correctly, but `conda/meta.yaml:2` had not yet been bumped from `2026.06.11`, so the conda-build run on the `v2026.6.12` tag produced a `gospl-2026.06.11-*.conda` artefact, which `anaconda upload --skip-existing` silently no-op'd against the already-published 2026-06-11 release. The recovery was `v2026.6.13` with `meson.build` + `conda/meta.yaml` synced. PyPI users on `gospl==2026.6.12` are functionally equivalent to `2026.6.13` (the only diff is the AGENTS.md / docs sweep below); conda users went straight from `2026.06.11` to `2026.6.13`.

### Local build and smoke-test procedure (osx-arm64)
Run this sequence from the repository root before pushing a release tag to
`geodels`. It validates the conda recipe, builds the package, and exercises
the full fast test suite against the *installed* package (not the source tree).
```bash
mamba install -n base -y conda-build        # one-time; skip if already present
conda build purge                           # clear stale build cache
conda build conda/ \
  -c conda-forge \
  --override-channels \
  --python 3.11 \
  --variants '{"python": ["3.11"]}' \
  2>&1 | tee build.log
mamba create -n gospl-smoke python=3.11 -c conda-forge -y
mamba install -n gospl-smoke \
   $PKG.conda \
  -c local -c conda-forge -y
mamba install -n gospl-smoke pytest -c conda-forge -y
cd /tmp
mamba run -n gospl-smoke python -c \
  "from gospl.model import Model; print('ok')"
mamba run -n gospl-smoke python -m pytest \
  /path/to/gospl/tests/ \
  -v --tb=short \
  --import-mode=importlib
cd -
mamba env remove -n gospl-smoke -y
```

### Known platform constraints (osx-arm64)
- **vtk → vtk-base**: the full `vtk` package on osx-arm64 pulls in `gtk3` /
  `gdk-pixbuf`, whose post-link script fails inside the conda-build sandbox.
  The recipe uses `vtk-base` (headless subset) instead. goSPL only uses VTK
  for unstructured mesh I/O, not rendering.
- **petsc4py ABI mismatch (historical)**: conda-forge previously published
  only a py310 (`np2py310`) osx-arm64 build of `petsc4py`, which worked for
  all computation but segfaulted during MPI finalization when spawned as a
  subprocess via `mpirun`. `test_parallel_correctness` detects this at runtime
  via `_petsc4py_abi_mismatch()` and skips rather than failing. As of the
  2026-06-11 recipe fix, conda-forge **does** publish py311 and py312
  openmpi-linked `petsc4py 3.21.2` builds (`py311h196a43b_0`,
  `py312ha15fc32_0`), which the recipe now pins — so the ABI-mismatch skip no
  longer fires on osx-arm64 and the test actually runs.
- **`test_parallel_correctness` nested-mpirun env leak (FIXED, regression-guarded)**:
  under OpenMPI, `import gospl` → `petsc4py.init()` → `MPI_Init` in the pytest
  parent exports `OMPI_*`/`PMIX_*`/`PRTE_*` env vars; the test's
  `subprocess.run(["mpirun", ...])` inherited them, so OpenMPI thought it was
  already inside an MPI job and silently refused to launch the nested run
  (`rc=1`, empty output). This was masked in CI until 2026-06 (ubuntu used
  MPICH, which is immune; osx-arm64 *skipped* the test on the old py310-only
  petsc4py). It surfaced on all cells once `environment.yml` pinned OpenMPI 4.x
  + petsc4py 3.21.x (real py311/py312 builds → no skip). **Fix:** `run_at_rank`
  now scrubs `OMPI_*`/`PMIX_*`/`PRTE_*`/`OPAL_*` from the child env (preserving
  `OPAL_PREFIX`) before spawning `mpirun`. Harmless under MPICH. The same leak
  affects the published conda package and HPC container (both OpenMPI), so this
  is a real fix, not just a CI patch.
- **Multi-version render**: pass `--variants '{"python": ["3.11"]}'` to
  conda-build to prevent it rendering the recipe for all Python versions in
  the conda-forge global pinnings file (which includes 3.13, incompatible
  with `numpy=1.26`).
- **OpenMPI 4.x pin (osx-arm64)**: the recipe pins `openmpi >=4.0,<5.0` in
  both `host` and `run`. `openmpi` is a pure transitive dependency (pulled in
  via `mpi4py`/`petsc4py`), so without an explicit pin conda-forge resolves to
  `openmpi 5.x`, which fails at `MPI_Init` on macOS with `PML add procs
  failed / Not found (-13)`. The 4.x line initialises cleanly. Because
  `openmpi` is otherwise invisible in the recipe, the pin must be listed
  explicitly — relying on a downstream package to constrain it does not work.
- **petsc4py 3.21.x ceiling (osx-arm64)**: the recipe pins both `petsc
  >=3.21,<3.22` and `petsc4py >=3.21,<3.22` in `host` and `run`. conda-forge
  `petsc4py >=3.22` is **mpich-only** on osx-arm64 — no openmpi-linked variant
  exists — so allowing >=3.22 silently drags in mpich and conflicts with the
  openmpi pin above. `3.21.2` is the last openmpi-linked `petsc4py` build
  published for osx-arm64.
- **h5py must be MPI-linked**: the recipe pins `h5py * mpi_openmpi*` and
  `hdf5 * mpi_openmpi*` in `run`. goSPL performs parallel (collective) HDF5
  writes via h5py; the conda-forge `nompi` variant of h5py is otherwise a
  valid solve and the solver will pick it, but it silently fails on collective
  writes under MPI. Pinning the `mpi_openmpi*` build string forces the
  MPI-linked variant and keeps it consistent with the openmpi 4.x line.
