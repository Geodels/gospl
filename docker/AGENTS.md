# HPC container (docker/)

Read before editing anything under `docker/` or `.github/workflows/docker-build.yml`.

This file is part of the goSPL agent guide. Read the root `AGENTS.md` first; it holds
the rules that apply everywhere (MPI contract, solver lifecycles, scratch Vecs,
conventions, commit checklist). The text below was moved verbatim from the root
file in the 2026-10 split, so dates and `## ...` cross-references still refer to
the original section names (most now live in the files listed in the root
`AGENTS.md` > Map).

## Docker / HPC Container

This section covers the Singularity/Apptainer container that packages goSPL for
NCI Gadi and Pawsey Setonix. Read it when working on anything under `docker/`.

**Note the MPI-stack difference from the conda package above.** The conda
package (osx-arm64) links **OpenMPI 4.x** and pins `h5py * mpi_openmpi*`. The
HPC container links **MPICH** (hybrid-MPI bind-mount, below) and builds h5py
against that MPICH. The `mpi_openmpi*` build string is OpenMPI-specific and
MUST NOT be carried into the container — OpenMPI is not MPICH-ABI-compatible.

### Target systems
| System | Scheduler | Container runtime | MPI stack | Notes |
|---|---|---|---|---|
| NCI Gadi | PBS | Singularity (module: `singularity`) | Intel MPI (MPICH ABI) | Cannot build images on Gadi — root required |
| Pawsey Setonix | Slurm | Singularity 4.1.0 (module: `singularity/4.1.0-mpi`) | Cray MPI (MPICH ABI) | Ubuntu 24.04 base mandatory (CPE 25.03+) |

### File layout
```
docker/
├── Dockerfile            Multi-stage: MPICH → PETSc → goSPL venv
├── Singularity.def       Native Apptainer definition (alternative build path)
├── build.sh              docker build → .sif conversion; run on local Linux
└── slurm/
    ├── gadi.pbs          PBS job script (NCI Gadi)
    └── setonix.slurm     Slurm job script (Pawsey Setonix)
```
A repo-root `.dockerignore` keeps `.git`/`results/`/caches out of the build
context (the Dockerfile `COPY . /opt/gospl-src`s the whole tree).

### MPI contract inside the container
The container ships MPICH as a build-time placeholder. At runtime, the cluster's
native high-performance MPI (Cray on Setonix, Intel MPI on Gadi) is bind-mounted
over it by the `-mpi` Singularity module flavour. This is the **hybrid MPI**
pattern — do not deviate from it. Consequences:

- `mpi4py`, `petsc4py`, **and parallel HDF5 + h5py** inside the container **must**
  be compiled from source against the container's MPICH (not pre-built
  conda/pip binary wheels). The Dockerfile and Singularity.def both do this
  explicitly (`--no-binary=mpi4py`, `--no-binary=petsc4py`, source HDF5 +
  `HDF5_MPI=ON --no-binary=h5py`, and `pip install --no-deps` for goSPL so the
  binary wheels are never pulled). Do not swap to binary wheels.
- The MPICH version in the Dockerfile (`ARG MPICH_VERSION`, currently `4.2.3`)
  must remain ABI-compatible with the cluster MPI. As of 2026-06, `mpich 4.2.x`
  is compatible with both Cray MPI and Intel MPI. If you bump it, verify ABI
  compatibility with Setonix's Cray MPI via Pawsey docs before merging.
- On Setonix, load `singularity/4.1.0-mpi`; this module sets `SINGULARITY_BINDPATH`
  and `SINGULARITYENV_LD_LIBRARY_PATH` automatically — no manual `--bind` for MPI
  paths is needed.

### goSPL-specific container invariants
These are required for correctness inside the container and must not be changed:

1. **`petsc4py.init` is called exactly once**, in `gospl/__init__.py`. The container
   does not alter this. Do not add `petsc4py.init()` calls to entrypoint scripts,
   Slurm wrappers, or any other layer around goSPL.

2. **Threading must be disabled inside each MPI rank.** Both the Dockerfile/def
   (`ENV`) and the job scripts export:
   ```bash
   export OMP_NUM_THREADS=1
   export OPENBLAS_NUM_THREADS=1
   export MKL_NUM_THREADS=1
   ```
   Do not remove these.

3. **Collective calls must span all ranks.** Nothing in the job scripts or
   entrypoint may cause one rank to diverge from the MPI contract described in
   `## MPI contract` above. In particular, do not redirect `stdout` of a subset
   of ranks before a `Barrier` or `Allreduce`.

4. **`destroy_DMPlex` must run on clean exit.** goSPL model runs should always
   call `model.destroy()` (the job scripts' `-c` driver does:
   `Model(...).runProcesses(); m.destroy()`). Do not `SIGKILL` the job before
   PETSc finalization unless debugging — leaked Vecs cause HPC sysadmin
   complaints and false leak reports.

5. **The scratch Vecs (`self.tmp`, `self.tmpL`, `self.tmp1`, `self.h`, `self.hl`,
   `self.dh`, etc.) are not safe to read after the step that wrote them.** If a
   post-processing or checkpointing script runs inside the container alongside
   goSPL, do not read these Vecs for elevation or cumulative ED — use
   `self.hLocal/hGlobal` and `self.cumED/cumEDLocal` respectively (see
   `## Scratch vector contract`).

### Known Setonix issue: parallel I/O inside Singularity
Pawsey has an active known issue: MPI-parallel HDF5 collective I/O inside
Singularity on Setonix can fail. goSPL uses collective HDF5 writes via the
MPI-linked h5py built into the container (against the container's MPICH — see
the MPI contract above). If you hit this on Setonix:
- First check Pawsey's Known Issues page for resolution status.
- Workaround: switch goSPL output to independent I/O mode if supported, or
  post-process outside the container after the simulation.
- Do NOT work around this by switching to the `nompi` h5py variant — goSPL
  requires MPI-linked h5py for collective writes. (Inside the container that
  means the MPICH-linked build, NOT the conda `mpi_openmpi*` build — see the
  note at the top of this section.)

### Build workflow (local Linux machine)
```bash
cd docker/
./build.sh v2026.7.14           # build + smoke-test + convert to .sif (no push)
./build.sh v2026.7.14 --push    # build + push Docker Hub + convert to .sif
```
The `build.sh` script:
1. Runs `docker build --platform linux/amd64` from the repo root with
   `--build-arg GOSPL_VERSION=<tag without leading v>`.
2. Smoke-tests the image (`import gospl; import petsc4py; import mpi4py` →
   `goSPL container build OK`) and asserts `gospl.__version__` matches the tag.
3. Optionally pushes to `docker.io/geodels/gospl-hpc` (override the repo with
   `REGISTRY=... ./build.sh ...`).
4. Converts to `.sif` via `apptainer build` / `singularity build` from the
   local Docker daemon (`docker-daemon://`), or from the registry (`docker://`)
   when `--push` was given.

**CI**: `.github/workflows/docker-build.yml` builds the same image on every `v*`
tag push and pushes `docker.io/geodels/gospl-hpc:<tag>` + `:latest` to Docker
Hub (needs `DOCKERHUB_USERNAME` / `DOCKERHUB_TOKEN` secrets). `workflow_dispatch`
builds + smoke-tests only (no push). The from-source MPICH/PETSc/HDF5 build is
slow (~30-60 min cold); the workflow uses the gha build cache.

**You cannot build on Gadi or Setonix.** Build the `.sif` locally and `scp` it,
or — once the image is on Docker Hub — `singularity pull docker://geodels/gospl-hpc:v2026.7.14`
on a login node (a conversion, not a root build).

numpy (`1.26.4`) and cython (`<3.1`) are pinned for the whole venv via a
constraints file wired to both `PIP_CONSTRAINT` and `PIP_BUILD_CONSTRAINT`
(both env vars because pip ≥26.2 stops honouring `PIP_CONSTRAINT` for build
deps): numpy stays on the conda/CI baseline (no transitive drift to 2.x), and
petsc4py 3.21.x cannot be cythonized by Cython ≥3.1 (`cyautodoc`
"ExpressionWriter" crash; h5py is fine with <3.1 too).

**setuptools is deliberately NOT pinned globally** — petsc4py 3.21.x and h5py
have *opposite* requirements: petsc4py's `confpetsc.py` needs the classic
`distutils.util.execute(dry_run=...)` that setuptools ≥74 dropped (py3.12 has no
stdlib distutils), while h5py's build requires setuptools ≥77.0.1. They're
resolved per-package: **petsc4py gets a dedicated build-constraint file** pinning
`setuptools<74` (passed as `PIP_BUILD_CONSTRAINT` for just that `pip install`),
while **h5py keeps build isolation** under the global (setuptools-free) constraint
and pulls its own ≥77. (Note: `PIP_BUILD_CONSTRAINT`/`--build-constraint` cannot
be combined with `--no-build-isolation` — pip rejects it — so petsc4py stays in
build isolation.) These build-toolchain contortions are a consequence of building
the older petsc4py 3.21.x from source; a future PETSc bump may let them be relaxed.

The base layers are split into cached stages — `mpich-build` → `petsc-build` →
`hdf5-build` (HDF5 chained `FROM petsc-build`). `docker-build.yml` warms them
with `--target hdf5-build` into a dedicated gha cache scope so a failure in the
python layers never triggers a cold ~45-min base rebuild.

### Version bumping
When bumping the goSPL version in the container, change **only** the
`GOSPL_VERSION` build arg (`build.sh` derives it from the `<tag>` you pass). Do
not hardcode version strings in the Dockerfile/def — goSPL's version derives
from `gospl.__version__`, set via `importlib.metadata` reading `meson.build`
line 4 (see `## MPI contract > __version__`). The build asserts that the
installed `gospl.__version__` equals `GOSPL_VERSION`, so the build arg and
`meson.build` must match on any published `.sif`.

### Checklist before publishing a new .sif
1. `./build.sh <tag>` completes without errors and the smoke test prints
   `goSPL container build OK`.
2. The MPICH version inside the container is ABI-compatible with the cluster MPI
   on both Gadi and Setonix.
3. `OMP_NUM_THREADS=1` / `OPENBLAS_NUM_THREADS=1` are exported in both job
   scripts.
4. The `.sif` filename includes the version tag (e.g. `gospl-hpc-v2026.7.14.sif`).
5. Update the `## Milestones` table in this file with the new tag and `.sif`
   publication date.
6. The `docker/slurm/gadi.pbs` and `docker/slurm/setonix.slurm` `CONTAINER=`
   paths reference the new `.sif` filename.
