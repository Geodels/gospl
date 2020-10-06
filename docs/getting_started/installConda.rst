.. _installConda:

=========================
Installation via Conda
=========================

Python version support
----------------------

Officially Python 3.7.1 and above, 3.8, and 3.9.

Installing pandas
-----------------

.. _install.anaconda:

Installing with Anaconda
~~~~~~~~~~~~~~~~~~~~~~~~

Installing pandas and the rest of the `NumPy <https://www.numpy.org/>`__ and
`SciPy <https://www.scipy.org/>`__ stack can be a little
difficult for inexperienced users.

The simplest way to install not only pandas, but Python and the most popular
packages that make up the `SciPy <https://www.scipy.org/>`__ stack
(`IPython <https://ipython.org/>`__, `NumPy <https://www.numpy.org/>`__,
`Matplotlib <https://matplotlib.org/>`__, ...) is with
`Anaconda <https://docs.continuum.io/anaconda/>`__, a cross-platform
(Linux, Mac OS X, Windows) Python distribution for data analytics and
scientific computing.

After running the installer, the user will have access to pandas and the
rest of the `SciPy <https://www.scipy.org/>`__ stack without needing to install
anything else, and without needing to wait for any software to be compiled.

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

Installing with Miniconda
~~~~~~~~~~~~~~~~~~~~~~~~~

The previous section outlined how to get pandas installed as part of the
`Anaconda <https://docs.continuum.io/anaconda/>`__ distribution.
However this approach means you will install well over one hundred packages
and involves downloading the installer which is a few hundred megabytes in size.

If you want to have more control on which packages, or have a limited internet
bandwidth, then installing pandas with
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

The next step is to create a new conda environment. A conda environment is like a
virtualenv that allows you to specify a specific version of Python and set of libraries.
Run the following commands from a terminal window::

    conda create -n name_of_my_env python

This will create a minimal environment with only Python installed in it.
To put your self inside this environment run::

    source activate name_of_my_env

On Windows the command is::

    activate name_of_my_env

The final step required is to install pandas. This can be done with the
following command::

    conda install pandas

To install a specific pandas version::

    conda install pandas=0.20.3

To install other packages, IPython for example::

    conda install ipython

To install the full `Anaconda <https://docs.continuum.io/anaconda/>`__
distribution::

    conda install anaconda

If you need packages that are available to pip but not conda, then
install pip, and then use pip to install those packages::

    conda install pip
    pip install django
