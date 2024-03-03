.. _installSrc:

=========================
Installation via Source
=========================

.. note::

  Below are some of the main instructions to build :mod:`gospl` from the git source tree. This approach is mainly for experienced users working on a Linux environment. It is highly recommended to use ``docker`` or ``conda``, for quick installation and for package and dependency updates.


Update System
--------------

  apt-get update -qq
  apt-get install -yq --no-install-recommends bash-completion build-essential
  apt-get install -yq --no-install-recommends python3-minimal python3-dev python3-pip
  apt-get install -yq --no-install-recommends python3-tk python3-dbg cmake
  apt-get install -yq --no-install-recommends python3-setuptools wget gfortran


MPICH
-------

::

      mkdir /tmp/mpich-build
      wget http://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-3.3.tar.gz
      tar xvzf mpich-3.3.tar.gz
      cd mpich-3.3
      ./configure --enable-fast=all,O3 --prefix=/opt/mpich
      make -j4
      make install
      ldconfig
      cd /tmp
      rm -fr *

      export MPI_DIR=/opt/mpich
      export PATH=/opt/mpich/bin:$PATH


PETSc
-------

::

      mkdir /tmp/petsc-build
      wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.13.5.tar.gz
      tar zxf petsc-lite-3.13.5.tar.gz
      cd petsc-3.13.5
      ./configure --with-debugging=0 --prefix=/opt/petsc
                  --COPTFLAGS="-g -O3" --CXXOPTFLAGS="-g -O3" --FOPTFLAGS="-g -O3"
                  --with-zlib=1
                  --download-fblaslapack=1
                  --download-ctetgen=1
                  --download-triangle=1
                  --download-hdf5=1
                  --download-mumps=1
                  --download-parmetis=1
                  --download-eigen=1
                  --download-metis=1
                  --download-hypre=1
                  --download-scalapack=1
                  --download-pragmatic=1
                  --useThreads=1
                  --with-shared-libraries
                  --with-cxx-dialect=C++11
      make PETSC_DIR=/tmp/petsc-build/petsc-3.13.5 PETSC_ARCH=arch-linux-c-opt all
      make PETSC_DIR=/tmp/petsc-build/petsc-3.13.5 PETSC_ARCH=arch-linux-c-opt install
      make PETSC_DIR=/opt/petsc PETSC_ARCH="" check
      cd /tmp
      rm -fr *
      export PETSC_DIR=/opt/petsc
      export PATH=/opt/petsc/bin:$PATH


Dependencies
----------------------

:mod:`gospl` has many required dependencies. If a
dependency is not installed, :mod:`gospl` will raise an ``ImportError`` when
the method/class requiring that dependency is called.

A dependency ``XXXX`` is installed via the following command in a terminal::

      pip install XXXX


========================= ================== =============================================================
Dependency                Minimum Version    Notes
========================= ================== =============================================================
NumPy                     1.22.2             Numerical computing tools.
SciPy                     1.5.2              Optimization, linear algebra, integration, interpolation
Cython                    0.29.21            Superset of the Python programming language
mpi4py                    3.1.1              Bindings for the Message Passing Interface standard
petsc4py                  3.13.0             Interface to PETSc libraries
h5py                      2.10.0             Interface to the HDF5 binary data format
pandas                    1.1.2              Data analysis and manipulation tool
ruamel.yaml               0.16.12            Parsing YAML to Python objects
meshio                    4.2.2              I/O for mesh files.
meshplex                  0.13.4             Fast tools for simplex meshes
pre-commit                2.7.1              Managing and maintaining multi-language pre-commit hooks
vtk                       9.0.3              Toolkit for 3D computer graphics and image processing
numpy-indexed             0.3.5              Functionality for indexed operations on numpy ndarrays
========================= ================== =============================================================


Setup install
----------------------

Once all the listed dependencies above have been installed, :mod:`gospl`
source files are available through `GitHub <https://github.com/Geodels/gospl>`_::

      git clone https://github.com/Geodels/gospl

It can then be installed locally on your system using::

      python setup.py install --user

If you wish to uninstall **gospl** you can do::

      python3 setup.py install --record gospl-files.txt

To record a list of installed files in ``gospl-files.txt``. Once you want to uninstall you can
use ``xargs`` to proceed with the uninstall::

      xargs rm -rf < gospl-files.txt
