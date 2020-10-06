.. _installConda:

=========================
Installation via Conda
=========================

.. _install.anaconda:

Installing Anaconda
--------------------------

Installing :mod:`gospl` and its dependencies stack can be a tedious and
difficult.

One of the simplest way to install not only :mod:`gospl`, but required Python
packages  is with `Anaconda <https://docs.continuum.io/anaconda/>`__, a cross-platform (Linux, Mac OS X, Windows) Python distribution for data analytics and
scientific computing.

After running the installer, the user will have already access to some essential Python packages and will be able to install a functioning :mod:`gospl` environment by following the directives below.

Installation instructions for `Anaconda <https://docs.continuum.io/anaconda/>`__
`can be found here <https://docs.continuum.io/anaconda/install.html>`__.

A full list of the packages available as part of the
`Anaconda <https://docs.continuum.io/anaconda/>`__ distribution
`can be found here <https://docs.continuum.io/anaconda/packages/pkg-docs/>`__.

Another advantage to installing Anaconda is that you don't need
admin rights to install it. Anaconda can install in the user's home directory,
which makes it trivial to delete Anaconda if you decide (just delete
that folder).

.. _install.miniconda:

Installing Miniconda
----------------------------

The previous section outlined how to get download the
`Anaconda <https://docs.continuum.io/anaconda/>`__ distribution.
However this approach means you will install well over one hundred packages
and involves downloading the installer which is a few hundred megabytes in size.

If you want to have more control on which packages, or have a limited internet
bandwidth, then installing :mod:`gospl` with
`Miniconda <https://conda.pydata.org/miniconda.html>`__ may be a better solution.

`Conda <https://conda.pydata.org/docs/>`__ is the package manager that the
`Anaconda <https://docs.continuum.io/anaconda/>`__ distribution is built upon.
It is a package manager that is both cross-platform and language agnostic
(it can play a similar role to a pip and virtualenv combination).

`Miniconda <https://conda.pydata.org/miniconda.html>`__ allows you to create a
minimal self contained Python installation, and then use the
`Conda <https://conda.pydata.org/docs/>`__ command to install additional packages.


First you will need `Conda <https://conda.pydata.org/docs/>`__ to be installed and
downloading and running the `Miniconda
<https://conda.pydata.org/miniconda.html>`__
will do this for you. The installer
`can be found here <https://conda.pydata.org/miniconda.html>`__

Building ``gospl`` environment
-------------------------------

The next step consists in downloading the conda environment for :mod:`gospl`.
A conda environment is like a virtualenv that allows you to specify a specific version of Python and set of libraries.
This is done by downloading the `conda-env.yml file <https://raw.githubusercontent.com/Geodels/gospl/master/conda-env.yml>`_. To do this you can use the ``curl``::

  curl https://raw.githubusercontent.com/Geodels/gospl/master/conda-env.yml --output conda-env.yml

or ``wget`` command::

  wget https://raw.githubusercontent.com/Geodels/gospl/master/conda-env.yml

This will save the file locally under the same name as it was on github: ``conda-env.yml``.

Alternatively you can get it from your preferred web browser by clicking on the following link: `conda-env.yml <https://raw.githubusercontent.com/Geodels/gospl/master/conda-env.yml>`_ and saving it under the following name ``conda-env.yml``.

.. note::

  :mod:`gospl` is not directly packaged as a `Conda <https://conda.pydata.org/docs/>`__ library because some of its dependencies are not available via this installation. The use of the environment file however provides an easy installation approach.

Once the `conda-env.yml <https://raw.githubusercontent.com/Geodels/gospl/master/conda-env.yml>`_ file has been downloaded on your system. The following directives provide a step-by-step guide to create a local conda environment for :mod:`gospl`.

Navigate to the directory containing the `conda-env.yml <https://raw.githubusercontent.com/Geodels/gospl/master/conda-env.yml>`_ file and run the following commands from a terminal window::

    conda env create -f conda-env.yml

This will create an environment with the dependencies and packages required to run :mod:`gospl`.

To put your self inside this environment run::

    source activate gospl-package

On Windows the command is::

    activate gospl-package

To install other packages, IPython for example::

    conda install ipython

To install the full `Anaconda <https://docs.continuum.io/anaconda/>`__
distribution::

    conda install anaconda

If you need packages that are available to ``pip`` but not ``conda``, then
the ``pip`` library is already installed, and can be used to install those packages::

    pip install django

To remove the environment, in your terminal window or an Anaconda Prompt, run::

    conda remove --name gospl-package --all


To verify that the environment was removed, in your terminal window or an Anaconda Prompt, run:

    conda info --envs


The ``gospl-package`` package should not be in the list of environment anymore.
