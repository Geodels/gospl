# syntax=docker/dockerfile:1
###############################################################################
# goSPL HPC container — hybrid-MPI (MPICH ABI) build for NCI Gadi & Pawsey
# Setonix.  See AGENTS.md > "Docker / HPC Container" for the full contract.
#
# Multi-stage:
#   1. mpich-build  — MPICH compiled from source (placeholder MPI)
#   2. petsc-build  — PETSc compiled against that MPICH
#   3. runtime      — goSPL venv; mpi4py / petsc4py / parallel-HDF5 + h5py all
#                     compiled FROM SOURCE against the container's MPICH
#
# At runtime the cluster's native MPI (Cray MPI on Setonix, Intel MPI on Gadi —
# both MPICH-ABI-compatible) is bind-mounted over /opt/mpich by the `-mpi`
# Singularity module flavour.  This is the documented hybrid-MPI pattern.
#
# DO NOT swap MPICH for OpenMPI: OpenMPI is NOT MPICH-ABI-compatible and the
# runtime bind-mount would break.  This is also why h5py here is built against
# the container's MPICH and NOT the conda `mpi_openmpi*` build (that build
# string is OpenMPI-specific and applies only to the osx-arm64 conda package).
###############################################################################

# Ubuntu 24.04 base is MANDATORY for Pawsey Setonix (CPE 25.03+).
ARG BASE_IMAGE=ubuntu:24.04

# MPICH 4.2.x is ABI-compatible with Cray MPI (Setonix) and Intel MPI (Gadi)
# as of 2026-06.  Verify against Pawsey docs before bumping (AGENTS.md).
ARG MPICH_VERSION=4.2.3
ARG PETSC_VERSION=3.21.6
ARG HDF5_VERSION=1.14.6
# numpy pinned to the conda/CI baseline (1.26.x). goSPL + pyshtools/vtk are
# validated against numpy 1.26; do not let pip drift to 2.x (see AGENTS.md).
ARG NUMPY_VERSION=1.26.4
# Full HDF5 source URL is parameterised because the HDF Group churns its
# release-asset naming; override with --build-arg if the default 404s.
ARG HDF5_URL=https://github.com/HDFGroup/hdf5/releases/download/hdf5_${HDF5_VERSION}/hdf5-${HDF5_VERSION}.tar.gz

###############################################################################
# Stage 1 — MPICH from source
###############################################################################
FROM ${BASE_IMAGE} AS mpich-build
ARG MPICH_VERSION
ENV DEBIAN_FRONTEND=noninteractive
# python3 is required by MPICH configure to generate the F90 bindings.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gfortran wget ca-certificates python3 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp/build
RUN wget -q "https://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz" \
    && tar xf "mpich-${MPICH_VERSION}.tar.gz" \
    && cd "mpich-${MPICH_VERSION}" \
    # ch4:ofi is the device the Cray/Intel MPI ABI shim expects.
    && ./configure \
        --prefix=/opt/mpich \
        --with-device=ch4:ofi \
        FFLAGS=-fallow-argument-mismatch \
        FCFLAGS=-fallow-argument-mismatch \
    && make -j"$(nproc)" \
    && make install \
    && cd / && rm -rf /tmp/build

###############################################################################
# Stage 2 — PETSc against the container MPICH
###############################################################################
FROM mpich-build AS petsc-build
ARG PETSC_VERSION
ENV PATH=/opt/mpich/bin:${PATH}
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gfortran wget ca-certificates python3 \
        libopenblas-dev liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp/build
RUN wget -q "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${PETSC_VERSION}.tar.gz" \
    && tar xf "petsc-${PETSC_VERSION}.tar.gz" \
    && cd "petsc-${PETSC_VERSION}" \
    && ./configure \
        --prefix=/opt/petsc \
        --with-mpi-dir=/opt/mpich \
        --with-debugging=0 \
        --with-fortran-bindings=0 \
        --with-shared-libraries=1 \
        COPTFLAGS="-O3" CXXOPTFLAGS="-O3" FOPTFLAGS="-O3" \
    && make PETSC_DIR="$(pwd)" PETSC_ARCH=arch-linux-c-opt all \
    && make PETSC_DIR="$(pwd)" PETSC_ARCH=arch-linux-c-opt install \
    && cd / && rm -rf /tmp/build

###############################################################################
# Stage 3 — parallel HDF5 against the container MPICH
# (chained FROM petsc-build so one `--target hdf5-build` warms the whole base:
#  MPICH + PETSc + HDF5. NOT the nompi variant — goSPL needs collective I/O.)
###############################################################################
FROM petsc-build AS hdf5-build
ARG HDF5_VERSION
ARG HDF5_URL
RUN apt-get update && apt-get install -y --no-install-recommends \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /tmp/build
RUN wget -q "${HDF5_URL}" -O hdf5.tar.gz \
    && tar xf hdf5.tar.gz \
    && cd "hdf5-${HDF5_VERSION}" \
    && CC=/opt/mpich/bin/mpicc ./configure \
        --prefix=/opt/hdf5 --enable-parallel --enable-shared \
    && make -j"$(nproc)" && make install \
    && cd / && rm -rf /tmp/build

###############################################################################
# Stage 4 — goSPL runtime venv
###############################################################################
FROM ${BASE_IMAGE} AS runtime
ARG PETSC_VERSION
ARG HDF5_VERSION
ARG HDF5_URL
ARG NUMPY_VERSION
# GOSPL_VERSION must match meson.build line 4; build.sh derives it from the tag.
ARG GOSPL_VERSION
ENV DEBIAN_FRONTEND=noninteractive

# Runtime + build toolchain (we compile mpi4py/petsc4py/h5py in this stage).
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-dev python3-venv python3-pip \
        build-essential gfortran pkg-config wget ca-certificates \
        libopenblas-dev liblapack-dev libfftw3-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Bring in the MPICH + PETSc + HDF5 trees built above.
COPY --from=mpich-build /opt/mpich /opt/mpich
COPY --from=petsc-build /opt/petsc /opt/petsc
COPY --from=hdf5-build  /opt/hdf5  /opt/hdf5

ENV PETSC_DIR=/opt/petsc \
    PATH=/opt/mpich/bin:/opt/venv/bin:${PATH} \
    LD_LIBRARY_PATH=/opt/mpich/lib:/opt/petsc/lib:/opt/hdf5/lib

# One thread per MPI rank — AGENTS.md container invariant 2. Job scripts
# re-export these; keeping them here makes a bare `singularity exec` safe too.
ENV OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

# Python venv (PATH already prepends /opt/venv/bin above).
ENV VIRTUAL_ENV=/opt/venv
# Pin numpy and cython for the whole venv, including build-isolation envs, via a
# constraints file. numpy stays at 1.26.x (scipy/pandas/vtk/pyshtools can't pull
# 2.x); cython is held <3.1 because petsc4py 3.21.x's cyautodoc cannot be
# cythonized by Cython >=3.1 (h5py is fine with <3.1 too). setuptools is NOT
# pinned globally — petsc4py 3.21.x needs an OLD setuptools (<74, for the classic
# distutils.util.execute(dry_run=...) signature) but h5py needs a NEW one
# (>=77.0.1), so the two can't share a global pin. We resolve this per-package:
# petsc4py gets a dedicated build-constraint file pinning setuptools<74 (below),
# while h5py keeps build isolation and pulls its own setuptools>=77.
# PIP_CONSTRAINT covers current pip; PIP_BUILD_CONSTRAINT covers pip >=26.2.
ENV PIP_CONSTRAINT=/opt/pip-constraints.txt \
    PIP_BUILD_CONSTRAINT=/opt/pip-constraints.txt
RUN python3 -m venv "${VIRTUAL_ENV}" \
    && printf 'numpy==%s\ncython<3.1\n' "${NUMPY_VERSION}" > /opt/pip-constraints.txt \
    && pip install --no-cache-dir --upgrade pip setuptools wheel \
    # numpy/cython for petsc4py & mpi4py builds; meson-python+ninja are the
    # goSPL build backend, needed because we install it --no-build-isolation.
    && pip install --no-cache-dir "numpy==${NUMPY_VERSION}" "cython<3.1" "meson-python>=0.15.0" ninja

# --- mpi4py from SOURCE against the container MPICH (never a binary wheel) ---
RUN MPICC=/opt/mpich/bin/mpicc \
    pip install --no-cache-dir --no-binary=mpi4py "mpi4py>=4.0"

# --- petsc4py from SOURCE against /opt/petsc (version matched to PETSc) ---
# petsc4py 3.21.x's confpetsc.py is incompatible with setuptools >=74 (the
# classic distutils.util.execute(dry_run=...) was dropped). Pin setuptools<74
# via a petsc4py-ONLY build constraint, so h5py below still gets its required
# setuptools>=77 from the global (setuptools-free) build constraint.
RUN printf 'numpy==%s\ncython<3.1\nsetuptools<74\n' "${NUMPY_VERSION}" > /opt/petsc4py-bc.txt \
    && PETSC_DIR=/opt/petsc PIP_BUILD_CONSTRAINT=/opt/petsc4py-bc.txt \
       pip install --no-cache-dir --no-binary=petsc4py "petsc4py==${PETSC_VERSION}"

# --- MPI-linked h5py against the HDF5 built in the hdf5-build stage ---
# (build isolation pulls h5py's own setuptools>=77; NOT the nompi variant.)
RUN CC=/opt/mpich/bin/mpicc HDF5_MPI=ON HDF5_DIR=/opt/hdf5 \
    pip install --no-cache-dir --no-binary=h5py h5py

# --- remaining pure/wheel deps (no direct MPI linkage) ---
RUN pip install --no-cache-dir \
        scipy pandas numpy-indexed "ruamel.yaml" \
        "gflex>=1.1.0,<2.0" pyshtools vtk

# --- install goSPL from the build context (repo root) ---
WORKDIR /opt/gospl-src
COPY . /opt/gospl-src
RUN pip install --no-cache-dir --no-build-isolation --no-deps .
# Version is driven solely by meson.build line 4; assert it matches the tag.
RUN test "$(python -c 'import gospl; print(gospl.__version__)')" = "${GOSPL_VERSION}" \
    || { echo "ERROR: container goSPL != ${GOSPL_VERSION}" >&2; exit 1; }

# Build-time smoke test (build.sh repeats this at runtime).
RUN python -c "import gospl, petsc4py, mpi4py; print('goSPL container build OK')"

LABEL org.opencontainers.image.title="gospl-hpc" \
      org.opencontainers.image.description="goSPL hybrid-MPI (MPICH ABI) HPC container for NCI Gadi / Pawsey Setonix" \
      org.opencontainers.image.source="https://github.com/Geodels/gospl" \
      org.opencontainers.image.version="${GOSPL_VERSION}"

WORKDIR /work
CMD ["python", "-c", "import gospl; print('goSPL', gospl.__version__)"]
