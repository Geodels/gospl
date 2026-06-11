.. _installConda:

=========================
Installation via Conda
=========================

.. _install.anaconda:


.. important::

    This is the preferred approach to install goSPL on a local machine.


.. warning::

    For **Windows users**, first install **Linux on Windows** with
    `WSL <https://learn.microsoft.com/en-us/windows/wsl/install>`_. This works
    with Windows 10 and 11.


.. _install.geodels-channel:


Quick install from the ``geodels`` conda channel
-------------------------------------------------

The simplest way to install goSPL is from the pre-built conda package published
on the ``geodels`` channel on `anaconda.org <https://anaconda.org/geodels/gospl>`_.
This avoids building goSPL's Fortran/Cython extensions locally and pulls every
runtime dependency from ``conda-forge`` in one step.

Create a fresh environment and install ``gospl`` in one command::

    mamba create -n gospl -c geodels -c conda-forge gospl

or with classic ``conda``::

    conda create -n gospl -c geodels -c conda-forge gospl

Then activate the environment::

    mamba activate gospl     # or: conda activate gospl

and verify the install::

    python -c "from gospl.model import Model; print('goSPL ready')"

.. note::

    ``-c conda-forge`` is required alongside ``-c geodels`` so that
    ``petsc4py``, ``mpi4py``, ``pyshtools``, ``vtk`` and the rest of the
    scientific stack resolve from ``conda-forge``. Without it, conda falls back
    to the ``defaults`` channel, which does not carry compatible builds of these
    packages.

To install a specific version (e.g. for reproducibility)::

    mamba create -n gospl -c geodels -c conda-forge gospl=2026.06.11

.. note::

    **Latest released packages on the** ``geodels`` **channel:**

    ========================= ============= ========================================
    Version                   Date          Install command
    ========================= ============= ========================================
    ``v2026.06.11``           2026-06-11    ``mamba install -c geodels -c conda-forge gospl``
    ``v2026.06.08``           2026-06-08    ``mamba install -c geodels -c conda-forge gospl=2026.06.08``
    ========================= ============= ========================================

The two sections below cover the alternative ``environment.yml`` approach, which
is useful if you want to add extra packages to the environment in the same step,
work from an editable source checkout, or track the ``master`` branch between
releases.


Installing Anaconda / Miniforge
--------------------------------

One of the simplest ways to install goSPL and its dependencies is with
`Anaconda <https://docs.continuum.io/anaconda/>`__ or
`Miniforge <https://github.com/conda-forge/miniforge>`__, cross-platform
(Linux, macOS, Windows) Python distributions for data analytics and scientific
computing.

.. note::

    **Recommended: use Miniforge.** Miniforge ships ``mamba`` (a fast drop-in
    replacement for ``conda``) and defaults to the ``conda-forge`` channel,
    which is where all of goSPL's dependencies live. Download from
    `github.com/conda-forge/miniforge <https://github.com/conda-forge/miniforge/releases>`__.

.. note::

    For **Windows users**, install Anaconda or Miniforge via the WSL terminal.
    Several guides are available online (e.g.
    `this one <https://emilykauffman.com/blog/install-anaconda-on-wsl>`_).

    For other approaches, installation instructions for
    `Anaconda <https://docs.continuum.io/anaconda/>`__ can be found
    `here <https://docs.continuum.io/anaconda/install.html>`__.

Another advantage of installing Anaconda/Miniforge is that you don't need admin
rights — they can be installed in your home directory.


.. _install.miniconda:

Installing Miniconda
---------------------

.. warning::

    If you have installed Anaconda or Miniforge already, skip this section.

If you want a minimal installation,
`Miniconda <https://conda.pydata.org/miniconda.html>`__ installs only the conda
package manager and a base Python environment. Download and run the installer
from `conda.pydata.org/miniconda.html <https://conda.pydata.org/miniconda.html>`__.
**Windows users** must run the installer from within the WSL terminal.


Building the goSPL environment
-------------------------------

The next step is to download the conda environment file for goSPL.

For the latest (``master`` branch) of goSPL, download ``environment.yml`` from
the repository::

  curl https://raw.githubusercontent.com/Geodels/gospl/master/environment.yml --output environment.yml

or::

  wget https://raw.githubusercontent.com/Geodels/gospl/master/environment.yml

Alternatively, download it directly from your browser:
`environment.yml <https://raw.githubusercontent.com/Geodels/gospl/master/environment.yml>`_.

The current ``environment.yml`` for the goSPL runtime environment is:

.. code-block:: yaml

    name: gospl
    channels:
        - conda-forge
    dependencies:
        - python>=3.11,<3.13
        - meson-python>=0.15.0
        - setuptools>=61.0
        - pkg-config
        - numpy
        - petsc>=3.21,<3.22
        - petsc4py>=3.21,<3.22
        - openmpi>=4.0,<5.0
        - mpi4py>=4.0
        - scipy
        - matplotlib
        - numpy-indexed
        - pandas
        - h5py * mpi_openmpi*
        - hdf5 * mpi_openmpi*
        - vtk-base
        - ruamel.yaml
        - gflex
        - pyshtools
        - cython
        - c-compiler
        - cxx-compiler
        - fortran-compiler
        - pytest

.. note::

    **Key dependency constraints** (do not remove these pins without reading
    the rationale):

    - ``petsc>=3.21,<3.22`` and ``petsc4py>=3.21,<3.22``: on macOS (osx-arm64)
      ``conda-forge`` stopped publishing ``openmpi``-linked ``petsc4py`` builds
      from 3.22 onward; pinning to 3.21.x keeps the OpenMPI-linked variant that
      is required for correct MPI finalisation.
    - ``openmpi>=4.0,<5.0``: OpenMPI 5.x fails at ``MPI_Init`` on macOS with
      a *PML add procs failed* error. Pin to the 4.x line on all platforms for
      consistency.
    - ``h5py * mpi_openmpi*`` and ``hdf5 * mpi_openmpi*``: goSPL performs
      collective (parallel) HDF5 writes; the ``nompi`` variant of h5py silently
      fails on collective writes under MPI.
    - ``vtk-base`` (not ``vtk``): the full ``vtk`` package on macOS/arm64 pulls
      in ``gtk3``/``gdk-pixbuf``, whose post-link script fails in the conda-build
      sandbox. goSPL uses VTK for unstructured mesh I/O only, not rendering.
    - Python is pinned ``>=3.11,<3.13`` (not a strict ``=3.11``). A strict pin
      conflicts with the CI matrix override on Python 3.12 cells.

.. note::

    ``pyshtools`` is required: the flexural-isostasy module uses a
    ``pyshtools``-based spherical-harmonic implementation. ``meshio`` and
    ``pyproj`` are **no longer** required dependencies and can be removed from
    any existing environment file.

.. note::

    goSPL is also available as a pre-built conda package on the ``geodels``
    channel — see :ref:`install.geodels-channel` above for the one-line install.
    The ``environment.yml`` approach is useful when you want to follow the
    ``master`` branch directly, work from an editable source checkout, or extend
    the environment with additional packages in one step.

Navigate to the directory containing ``environment.yml`` and run::

    mamba env create -f environment.yml

or with classic ``conda``::

    conda env create -f environment.yml

This creates an environment with all dependencies required to run goSPL.

Activate the environment::

    conda activate gospl

Install additional packages if needed::

    conda install jupyter

If a package is only available via ``pip``::

    pip install <package-name>

Remove the environment::

    conda remove --name gospl --all

Verify removal::

    conda info --envs


Alternative goSPL installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To try another branch or version of goSPL, activate the environment and run::

    pip install git+https://github.com/Geodels/gospl.git@NAME

where ``NAME`` is the branch name or version tag (e.g. ``v2026.06.11``,
``master``).

Alternatively, clone the repository and install in editable mode::

    git clone https://github.com/Geodels/gospl.git
    cd gospl
    pip install --no-build-isolation -e .

The ``--no-build-isolation`` flag is required so that pip uses the MPI, PETSc,
and compiler toolchain already present in the conda environment rather than
fetching a fresh isolated build environment.


Building the goSPL full-stack environment
------------------------------------------

The full-stack environment adds pre- and post-processing libraries on top of
the runtime dependencies. It is used to run all of the
`goSPL-examples <https://github.com/Geodels/goSPL-examples>`_ notebooks.

.. code-block:: yaml

    name: gospl-smoke
    channels:
        - conda-forge
    dependencies:
        - python>=3.11,<3.13
        - meson-python>=0.15.0
        - setuptools>=61.0
        - pkg-config
        - numpy
        - petsc>=3.21,<3.22
        - petsc4py>=3.21,<3.22
        - openmpi>=4.0,<5.0
        - mpi4py>=4.0
        - pip
        - scipy
        - matplotlib
        - numpy-indexed
        - pandas
        - h5py * mpi_openmpi*
        - hdf5 * mpi_openmpi*
        - meshio
        - vtk-base
        - pre-commit
        - ruamel.yaml
        - cython
        - c-compiler
        - cxx-compiler
        - fortran-compiler
        - meshplex
        - gflex
        - netcdf4
        - xarray
        - pyshtools
        - uxarray
        - pyinterp
        - jigsawpy
        - stripy
        - xesmf
        - mpas_tools
        - pygmt
        - rasterio
        - pysheds
        - seaborn
        - pyevtk
        - numba
        - shapely
        - pyvista
        - pyproj
        - triangle
        - pytest

Install with::

    mamba env create -f environment.yml

followed by::

    pip install --no-build-isolation -e /path/to/gospl

This environment allows you to run all examples provided in the
`goSPL-examples repository <https://github.com/Geodels/goSPL-examples>`_.
