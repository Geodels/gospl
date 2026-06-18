# Local conda recipe for goSPL

This directory contains a [conda-build](https://docs.conda.io/projects/conda-build/) recipe that produces a local goSPL conda package. It is intended for **local testing only** — submitting to conda-forge would use a different layout (under their `staged-recipes` repo) and would need pinning against the conda-forge global pinnings file.

## Supported platforms

Built and tested on:

- `linux-64`
- `osx-arm64`

Windows is skipped because the MPI toolchain that goSPL depends on (`mpi4py`, `petsc4py`) is not well supported by conda-forge on Windows.

## Prerequisites

You need conda and conda-build:

```bash
conda install -n base conda-build
```

The recipe pulls all build- and host-time dependencies from `conda-forge`, so make sure conda-forge is in your channel list:

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
```

## Building locally

From the repository root:

```bash
conda build conda/ -c conda-forge --override-channels --python 3.11
```

The `-c conda-forge --override-channels` is required because `petsc4py`, `mpi4py` and `pyshtools` are only available on conda-forge. The `--python 3.11` is required because some host dependencies (notably the `numpy 1.26` builds that match `petsc4py`'s ABI) are not yet built for newer Python versions on osx-arm64 at the time of writing.

This will:

1. Pull the build-stage compilers (`{{ compiler('c') }}`, `{{ compiler('fortran') }}`) from conda-forge.
2. Stage host dependencies (`python`, `numpy`, `cython`, `meson-python`, `mpi4py`, `petsc4py`, ...).
3. Run `python -m pip install . --no-build-isolation --no-deps -vv` to drive the meson-python backend.
4. Verify the install by importing `gospl` and `gospl._fortran`.

The finished `.conda` file is placed under `${CONDA_PREFIX}/conda-bld/<platform>/`. The exact path is printed by conda-build at the end.

## Validating the recipe without building

A full build takes 15–30 minutes and a few hundred MB of downloads. For a quick syntax-and-solver check, use `conda render`:

```bash
conda render conda/ -c conda-forge --override-channels --python 3.11
```

This parses the recipe, runs the Jinja templating, and asks the solver to resolve all build/host/run dependencies — but does not compile anything.

## Trying different Python versions

```bash
conda build conda/ -c conda-forge --override-channels --python 3.11
conda build conda/ -c conda-forge --override-channels --python 3.12
```

## Testing the built package

```bash
conda install -c local gospl
python -c "import gospl; print(gospl.__file__)"
```

For a quick MPI smoke test:

```bash
mpirun -n 2 python -c "from mpi4py import MPI; print(MPI.COMM_WORLD.Get_rank())"
```

## Notes

- The version string in `meta.yaml` (currently `2026.06.08`) is **hard-coded** to match `meson.build` at the time of writing. If you bump the version in `meson.build`, update `meta.yaml` AND `docs/conf.py` to match — there is no automatic syncing across these files.
- `pyshtools` is required by the global (spherical) flexure backend (spherical-harmonic implementation). The flat-model flexure is solved directly on the unstructured mesh and needs no external library. `meshio`, `pyproj` and the former `gflex` dependency are no longer required.
- The recipe builds against the `petsc4py` package on conda-forge; if you have a custom PETSc build you want to link against instead, build goSPL from source via `pip install --no-build-isolation -e .` rather than through this recipe.
