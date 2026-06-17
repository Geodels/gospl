.. _installSrc:

=========================
Installation via Source
=========================

.. note::

  The instructions below describe how to build goSPL from the git source tree.
  This approach is mainly for **experienced users** on a Linux environment.
  For most use-cases, the :ref:`install.geodels-channel` one-line conda install
  or the Docker image are simpler and faster.


Update system
--------------

::

      apt-get update -qq
      apt-get install -yq --no-install-recommends \
          bash-completion build-essential gfortran cmake wget \
          python3-minimal python3-dev python3-pip python3-venv \
          python3-tk python3-dbg python3-setuptools \
          pkg-config proj-bin


MPICH
------

goSPL requires an MPI implementation. MPICH 4.2.x is recommended; it is
ABI-compatible with Intel MPI (Gadi) and Cray MPI (Setonix) for HPC deployments.

::

      MPICH_VERSION=4.2.3
      wget https://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz
      tar xvzf mpich-${MPICH_VERSION}.tar.gz
      cd mpich-${MPICH_VERSION}
      ./configure \
          --enable-fast=all,O3 \
          --prefix=/opt/mpich \
          --enable-shared \
          --enable-romio \
          FFLAGS="-fallow-argument-mismatch"
      make -j$(nproc)
      make install
      ldconfig
      cd /tmp && rm -rf mpich-${MPICH_VERSION}*

      export MPI_DIR=/opt/mpich
      export PATH=/opt/mpich/bin:$PATH
      export LD_LIBRARY_PATH=/opt/mpich/lib:$LD_LIBRARY_PATH


PETSc
------

goSPL requires PETSc ≥ 3.21 (3.21.x recommended for compatibility with
conda-forge ``petsc4py`` builds).

::

      PETSC_VERSION=3.21.6
      wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${PETSC_VERSION}.tar.gz
      tar zxf petsc-${PETSC_VERSION}.tar.gz
      cd petsc-${PETSC_VERSION}
      python3 ./configure \
          --prefix=/opt/petsc \
          --with-mpi-dir=/opt/mpich \
          --with-debugging=0 \
          --with-shared-libraries=1 \
          --download-fblaslapack \
          COPTFLAGS="-O3" \
          CXXOPTFLAGS="-O3" \
          FOPTFLAGS="-O3"
      make PETSC_DIR=$(pwd) PETSC_ARCH=arch-linux-c-opt all
      make PETSC_DIR=$(pwd) PETSC_ARCH=arch-linux-c-opt install
      cd /tmp && rm -rf petsc-${PETSC_VERSION}*

      export PETSC_DIR=/opt/petsc
      export PATH=/opt/petsc/bin:$PATH
      export LD_LIBRARY_PATH=/opt/petsc/lib:$LD_LIBRARY_PATH


Dependencies
-------------

goSPL has the following required dependencies. Install them into a virtual
environment (recommended) or your system Python::

      python3 -m venv /opt/gospl-env
      source /opt/gospl-env/bin/activate
      pip install --upgrade pip wheel setuptools

Install the build toolchain (the version pins matter — see the note below) and
build ``mpi4py`` and ``petsc4py`` from source against the MPICH/PETSc installed
above::

      pip install "cython<3.1" "setuptools<74" "meson-python>=0.15.0" ninja numpy
      MPICC=/opt/mpich/bin/mpicc pip install --no-binary mpi4py "mpi4py>=4.0"
      PETSC_DIR=/opt/petsc pip install --no-build-isolation "petsc4py==3.21.6"

.. note::

   ``petsc4py`` 3.21.x cannot be built with **Cython ≥ 3.1** (its ``cyautodoc``
   hook crashes with ``'ExpressionWriter' object has no attribute
   'emit_string'``) nor with **setuptools ≥ 74** (its ``confpetsc.py`` calls the
   classic ``distutils.util.execute(dry_run=...)`` that newer setuptools dropped,
   and Python 3.12 has no stdlib ``distutils`` fallback). Pin ``cython<3.1`` and
   ``setuptools<74`` in the environment and build petsc4py with
   ``--no-build-isolation`` so it uses them. The ``petsc4py`` version **must**
   match your PETSc install (``3.21.6`` here).

Then install the remaining runtime dependencies::

      pip install \
          numpy scipy \
          h5py netCDF4 \
          pandas ruamel.yaml \
          numpy-indexed \
          pyshtools \
          vtk pyevtk \
          stripy meshplex \
          xarray matplotlib \
          pytest

.. note::

   ``meshio`` and ``pyproj`` are **no longer required** by goSPL. The flexure
   module uses a ``pyshtools``-based spherical-harmonic implementation; if your
   environment still lists ``isoFlex``, it can be removed.


Minimum supported versions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

========================= ================== =============================================================
Dependency                Minimum Version    Notes
========================= ================== =============================================================
NumPy                     1.26               Numerical computing tools.
SciPy                     1.11               Optimisation, linear algebra, integration, interpolation.
Cython                    3.0                Superset of Python; used by goSPL's extensions. Pin ``<3.1`` to build petsc4py 3.21.x from source (see note above).
mpi4py                    4.0                Bindings for MPI; **must** be compiled against your MPICH/OpenMPI.
petsc4py                  3.21               PETSc Python interface; **must** match your PETSc install.
h5py                      3.10               HDF5 interface; **must** be the MPI-linked build for parallel I/O.
pandas                    2.0                Data analysis; forcing DataFrames use named column access.
ruamel.yaml               0.18               YAML parser for goSPL input files.
vtk                       9.3                VTK for unstructured mesh I/O (headless; ``vtk-base`` on conda).
numpy-indexed             0.3.7              Indexed operations on NumPy arrays.
pyshtools                 4.13               Spherical-harmonic backend for the flexure module.
meson-python              0.16               Build backend (replaces ``setup.py``).
========================= ================== =============================================================


goSPL installation
-------------------

Once all dependencies are installed, clone the goSPL repository::

      git clone https://github.com/Geodels/gospl.git
      cd gospl

Install in editable mode::

      pip install --no-build-isolation -e .

The ``--no-build-isolation`` flag is required so that pip uses the MPI, PETSc,
and compiler toolchain already present in your environment rather than fetching a
fresh isolated build environment (which would not find your local PETSc install).

Verify the installation::

      python -c "from gospl.model import Model; import gospl; print(gospl.__version__)"

.. note::

   The version string is read from the installed package metadata, which is
   driven by ``meson.build`` line 4. If goSPL was cloned but not installed
   (bare ``git clone`` without ``pip install -e .``), ``gospl.__version__``
   returns ``"unknown"`` — this is expected.
