.. _installPip:

==========================
Installation via PyPI
==========================

.. _install.pypi:

goSPL is published on `PyPI <https://pypi.org/project/gospl/>`_ as a **source
distribution** (sdist). The pip path compiles goSPL's Fortran extension on
install against the user's existing scientific Python stack — there is **no
binary wheel**, by design.

.. important::

    Pip will **not** resolve goSPL's compiled MPI dependencies for you.
    ``petsc4py``, ``mpi4py``, and a parallel-HDF5-linked ``h5py`` must already
    be installed and linked against the same MPI implementation. The clean
    way to satisfy that is a conda environment (see :ref:`installConda`); the
    pip install then drops goSPL on top of it. If you do not already have a
    working MPI / PETSc stack, **start with the conda install instead** — it
    is faster, more reliable, and avoids hand-pinning compiler toolchain
    versions.

When to use the pip install
----------------------------

This approach is the right fit when:

- You already have a working ``mpi4py`` / ``petsc4py`` / parallel ``h5py``
  stack (typically from conda) and want to add goSPL to it.
- You are deploying on an HPC login node that exposes a system MPI + PETSc
  via modules.
- You want the smallest possible artefact (no extra runtime deps beyond
  what is already in your env).

It is **not** the right fit when:

- You are on a fresh machine with no MPI installed — use the conda path.
- You are on Windows without WSL — pip cannot build the Fortran extension
  natively on Windows.
- You need a fully reproducible binary build — use the conda package.


Quick install
-------------

In a conda environment that already supplies the goSPL runtime dependencies
(see :ref:`installConda` for the ``environment.yml`` recipe), install goSPL
from PyPI::

    pip install --no-deps --no-build-isolation gospl

Then verify::

    python -c "import gospl; print(gospl.__version__)"

Why each flag matters:

- ``--no-deps`` — the conda environment already supplies every runtime
  dependency with the right MPI / HDF5 build strings (``mpi_openmpi*`` on
  conda-forge). Letting pip resolve dependencies from PyPI would overwrite
  them with non-MPI variants and silently break collective HDF5 writes.
- ``--no-build-isolation`` — uses the conda environment's existing
  ``numpy``, ``meson-python``, and Fortran compiler instead of building a
  fresh isolated venv. Faster and avoids ABI drift between the build-time
  and run-time numpy.


Installing a specific version
------------------------------

To pin to a specific release::

    pip install --no-deps --no-build-isolation gospl==2026.6.13

Pre-release / **release candidate** builds (suffixed ``rcN``) are not
installed by default; request them explicitly with ``--pre``::

    pip install --pre --no-deps --no-build-isolation gospl

.. note::

    **Latest released package on PyPI:**

    ========================= ============= ====================================
    Version                   Date          Install command
    ========================= ============= ====================================
    ``2026.6.30``             2026-06-30    ``pip install --no-deps --no-build-isolation gospl``
    ``2026.6.13``             2026-06-12    ``pip install --no-deps --no-build-isolation gospl==2026.6.13``
    ``2026.6.12``             2026-06-12    ``pip install --no-deps --no-build-isolation gospl==2026.6.12``
    ``2024.9.1``              2024-08-24    *historical release*
    ========================= ============= ====================================

    Earlier releases (``< 2026.6.x``) predate the current Fortran/Cython
    toolchain pins and are not recommended for new work.


Bring-your-own dependency setup
-------------------------------

If you are not using conda, you must ensure that the following are present
**before** running ``pip install gospl``:

- A working MPI implementation (OpenMPI ≥ 4 or MPICH ≥ 4.2), with
  ``mpicc`` / ``mpifort`` on ``PATH``.
- PETSc ≥ 3.21 built against that MPI, with ``PETSC_DIR`` set.
- HDF5 built with parallel I/O (``--enable-parallel``).
- ``mpi4py`` and ``petsc4py`` compiled against your MPI / PETSc (NOT the
  pre-built PyPI wheels — they will not match your MPI ABI).
- ``h5py`` built against the parallel HDF5 (``HDF5_MPI=ON``).
- ``gfortran`` (or compatible Fortran compiler) on ``PATH`` so meson-python
  can build goSPL's ``_fortran`` extension at install time.

This bring-your-own setup is the same constraint that the HPC container
satisfies (see :ref:`installHPC`); if you find yourself building all of
the above by hand, consider using the published Singularity image instead.

See :ref:`installSrc` for a complete from-scratch source build that
includes the MPICH + PETSc + HDF5 + numpy / petsc4py / h5py build steps.


Upgrading
---------

To upgrade to the latest PyPI release inside an existing environment::

    pip install --no-deps --no-build-isolation --upgrade gospl

To install from a TestPyPI pre-release (during release validation)::

    pip install --no-deps --no-build-isolation \
        --index-url https://test.pypi.org/simple/ \
        --extra-index-url https://pypi.org/simple/ \
        gospl

.. note::

    The PyPI sdist is built by the ``pypi-publish`` GitHub Actions workflow
    on every release tag (``v*``) and uses **Trusted Publishing** (OIDC,
    no API tokens). The version on PyPI is always normalized per
    PEP 440 — leading zeros on month / day are stripped — so
    ``gospl 2026.6.13`` and ``gospl 2026.06.13`` resolve to the same
    distribution. goSPL writes the no-leading-zero form throughout the
    project so the PyPI display, conda artefact name, git tag, and
    ``gospl.__version__`` all match exactly.


Editable install for development
---------------------------------

For local development against the git source tree (rather than a released
sdist), use the source-checkout path::

    git clone https://github.com/Geodels/gospl.git
    cd gospl
    pip install --no-deps --no-build-isolation -e .

This is the same command CI runs on every PR. The ``--no-deps`` and
``--no-build-isolation`` flags carry the same meaning as the PyPI install
above. Editable mode means changes to the Python sources take effect
immediately; the Fortran extension still has to be rebuilt
(``pip install --no-build-isolation --force-reinstall --no-deps -e .``)
when you change ``fortran/functions.F90`` or ``fortran/functions.pyf``.

See :ref:`installSrc` for the full from-scratch build (system Python,
manual MPICH + PETSc + HDF5) used when no conda environment is available.
