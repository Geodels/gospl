.. _installPypi:


The following installation guidelines are mainly for experienced users working on a Linux environment. It is highly recommended to use ``docker`` or ``conda``, for quick installation and for package and dependency updates.
You can find simple installation instructions for :mod:`gospl` in this document.

=========================
Installation via PyPI
=========================

:mod:`gospl` can be installed via pip from
`PyPI <https://pypi.org/project/gospl>`__.

::

    pip install gospl

The installation however requires additional dependencies and is likely to
failed if some initial packages are not already available. When using ``pip``,
it is recommended to use ``virtualenv``. For a guide on how to use it, readers
are invited to look at visit the `Python documentation <https://docs.python.org/3/tutorial/venv.html>`_.


Installing dependencies
-----------------------------


On linux, you will need to install ``MPI``, ``PETSc`` and ``Hdf5`` via the your Linux distribution's package manager.


For example, the commands in this table will install ``MPICH`` from your distribution.

.. csv-table::
    :header: "Distribution", "Status", "Download / Repository Link", "Install method"
    :widths: 10, 10, 20, 50


    Debian, stable, `official Debian repository <https://packages.debian.org/search?keywords=mpich&searchon=names&suite=all&section=all>`__ , ``sudo apt-get install mpich``
    Debian & Ubuntu, unstable (latest packages), `NeuroDebian <http://neuro.debian.net/index.html#how-to-use-this-repository>`__ , ``sudo apt-get install mpich``
    Ubuntu, stable, `official Ubuntu repository <https://packages.ubuntu.com/search?keywords=mpich&searchon=names&suite=all&section=all>`__ , ``sudo apt-get install mpich``
    OpenSuse, stable, `OpenSuse Repository  <https://software.opensuse.org/download/package?package=mpich&project=openSUSE%3A12.1>`__ , ``zypper in mpich``
    Fedora, stable, `official Fedora repository  <https://fedora.pkgs.org/30/fedora-x86_64/mpich-3.2.1-9.fc30.i686.rpm.html>`__ , ``dnf install mpich``
    Centos/RHEL, stable, `EPEL repository <https://centos.pkgs.org/6/centos-x86_64/mpich-3.1-5.el6.x86_64.rpm.html>`__ , ``yum install mpich``



Using the ``apt-get`` command. This will likely be done by::

    sudo apt-get update
    sudo apt-get install mpich libmpich-dev libhdf5-mpich-dev
    sudo apt-get install libblas-dev liblapack-dev
    sudo apt-get install libscotchparmetis-dev libmetis-dev


If you have ``PETSc`` already installed, you can point to your ``PETSC_DIR`` and  ``PETSC_ARCH`` or you can unset these paths before installing the ``petsc4py`` module via ``pip``::

    unset PETSC_DIR
    unset PETSC_ARCH


Installing ``PETSc`` and ``gospl``
-----------------------------------

Then the remaining library are installed locally by running::

    pip install --upgrade pip setuptools wheel
    pip install petsc --user
    pip install petsc4py --user
    pip install gospl --user
