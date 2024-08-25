.. _installConda:

Installing goSPL and its dependencies stack from source can be **tedious and difficult**.


=========================
Installation via Conda
=========================

.. _install.anaconda:


.. important::

    This is the preferred approach to install goSPL on a local machine.


.. warning::

    For **Windows users**, first install **Linux on Windows** with `WSL <https://learn.microsoft.com/en-us/windows/wsl/install>`_. This should work with most recent versions of Windows (from Windows 10 to 11). 


Installing Anaconda
--------------------------

One of the simplest way to locally install not only goSPL, but required dependencies  is with `Anaconda <https://docs.continuum.io/anaconda/>`__, a cross-platform (Linux, Mac OS X, Windows) Python distribution for data analytics and scientific computing.


.. note::

    For **Windows users**, you will need to install Anaconda via the Windows Ubuntu Terminal from WSL. There are several articles on the web to do so (such as this `one <https://emilykauffman.com/blog/install-anaconda-on-wsl>`_)

    For other approaches, some installation instructions for `Anaconda <https://docs.continuum.io/anaconda/>`__ can be found `here <https://docs.continuum.io/anaconda/install.html>`__.

A full list of the packages available as part of the `Anaconda <https://docs.continuum.io/anaconda/>`__ distribution can be found `here <https://docs.continuum.io/anaconda/packages/pkg-docs/>`__.

Another advantage of installing Anaconda is that you don't need admin rights and it can be installed in the user's home directory, which makes it trivial to delete Anaconda if you decide (just delete that folder).

After getting Anaconda installed, the user will have already access to some essential Python packages and will be able to install a functioning goSPL environment by following the directives below.


.. _install.miniconda:

Installing Miniconda
----------------------------


.. warning::

    If you have installed Anaconda already no need to install Miniconda and you can skip this section.
    

The previous section outlined how to download the `Anaconda <https://docs.continuum.io/anaconda/>`__ distribution. However this approach means you will install well over one hundred packages and involves downloading the installer which is a few hundred megabytes in size.

If you want to have more control on which packages, or have a limited internet
bandwidth, then installing goSPL with `Miniconda <https://conda.pydata.org/miniconda.html>`__ may be a better solution.

`Conda <https://conda.pydata.org/docs/>`__ is the package manager that the `Anaconda <https://docs.continuum.io/anaconda/>`__ distribution is built upon. It is a package manager that is both cross-platform and language agnostic (it can play a similar role to a ``pip`` and ``virtualenv`` combination).

`Miniconda <https://conda.pydata.org/miniconda.html>`__ allows you to create a minimal self contained Python installation, and then use the `conda` command to install additional packages.


First you will need `Conda <https://conda.pydata.org/docs/>`__ to be installed and downloading and running the `Miniconda <https://conda.pydata.org/miniconda.html>`__
will do this for you. The installer can be found `here <https://conda.pydata.org/miniconda.html>`__. It is worth mentioning that for **Windows users**, the miniconda installation will need to be done through the WSL terminal.

Building goSPL environment
-------------------------------

The next step consists in downloading the conda environment for goSPL. A conda environment is like a virtual environment that allows you to install a specific flavor of Python and set of libraries. For the latest version (`master` branch) of goSPL, this is done by downloading the ``environment.yml`` `file <https://raw.githubusercontent.com/Geodels/gospl/master/environment.yml>`_. To do this you can use the ``curl``::

  curl https://raw.githubusercontent.com/Geodels/gospl/master/environment.yml --output environment.yml

or ``wget`` command::

  wget https://raw.githubusercontent.com/Geodels/gospl/master/environment.yml

This will save the file locally under the same name as it was on github: ``environment.yml``.

.. code-block:: yaml

    name: gospl
    channels:
        - conda-forge
        - defaults
    dependencies:
        - python=3.11
        - meson-python>=0.15.0
        - setuptools>=61.0
        - pkg-config
        - build
        - numpy
        - petsc4py
        - scipy
        - numpy-indexed
        - pandas
        - h5py
        - meshio
        - vtk
        - ruamel.yaml
        - gflex
        - mpi4py
        - xarray
        - rioxarray
        - pyproj
        - cython
        - compilers
        - pip:
            - git+https://github.com/Geodels/isoFlex.git
            - git+https://github.com/Geodels/gospl.git


Alternatively you can get it from your preferred web browser by clicking on the following link: `environment.yml <https://raw.githubusercontent.com/Geodels/gospl/master/environment.yml>`_ and saving it under the following name ``environment.yml``.

.. note::

  goSPL is not directly packaged as a `Conda <https://conda.pydata.org/docs/>`__ library because some of its dependencies are not available via this installation. The use of the environment file however provides an easy installation approach.

Once the `environment.yml <https://raw.githubusercontent.com/Geodels/gospl/master/environment.yml>`_ file has been downloaded on your system. The following directives provide a step-by-step guide to create a local conda environment for goSPL.

Navigate to the directory containing the `environment.yml <https://raw.githubusercontent.com/Geodels/gospl/master/environment.yml>`_ file and run the following commands from a terminal window::

    conda env create -f environment.yml

This will create an environment with the dependencies and packages required to run goSPL.

To put yourself inside this environment run::

    source activate gospl


To install other packages, jupyter for example::

    conda install jupyter

If you need packages that are available via ``pip`` but not ``conda``, then
the ``pip`` library is already installed, and can be used to install those packages::

    pip install django

To remove the environment, in your terminal window or an Anaconda Prompt, run::

    conda remove --name gospl --all


To verify that the environment was removed, in your terminal window or an Anaconda Prompt, run::

    conda info --envs


The ``gospl`` conda environment should not be in your list of environment anymore.


Alternative goSPL installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You might want to try another branch/version of goSPL. To do so you could once your conda environment has been activated run the following::

    pip install git+https://github.com/Geodels/gospl.git@NAME


where ``NAME`` needs to be replaced by the branch/version you want to try.

Alternatively, you could clone or download the goSPL `repository <https://github.com/Geodels/gospl/archive/refs/heads/master.zip>`_ and run the following command in the repository directory::

    pip install --no-build-isolation -e .



Building goSPL full stack environment
----------------------------------------

The full stack environment contains not only the libraries used for running goSPL, but also the ones for making the **pre- and post-processing**. It can be installed in a similar fashion as the one described above using the following ``environment.yml`` file.

.. code-block:: yaml

    name: gospl
    channels:
        - conda-forge
        - defaults
    dependencies:
        - python=3.11
        - meson-python>=0.15.0
        - setuptools>=61.0
        - pkg-config
        - build
        - numpy
        - petsc4py
        - pip
        - scipy
        - numpy-indexed
        - pandas
        - h5py
        - meshio
        - vtk
        - pre-commit
        - ruamel.yaml
        - mpi4py
        - cython
        - compilers
        - meshplex
        - gflex
        - netcdf4
        - xarray
        - rioxarray
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
        - basemap
        - numba
        - shapely
        - pyvista
        - pyproj 
        - triangle
        - pip:
            - git+https://github.com/Geodels/isoFlex.git
            - git+https://github.com/Geodels/gospl.git

This conda environment will allow you to run all the examples provided in this `repository <https://github.com/Geodels/goSPL-examples>`_. 
