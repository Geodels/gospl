.. _getting_started:

===============
Getting started
===============


Python version support
----------------------

goSPL supports Python :mod:`3.11` and :mod:`3.12`.

Installation instructions
-------------------------

The easiest way to start using goSPL is to install it from the pre-built
``geodels`` conda package or from its Docker image. Both approaches avoid
building goSPL's Fortran/Cython extensions locally. Instructions for
installing from source and deploying on HPC systems (NCI Gadi, Pawsey
Setonix) are also provided.


.. grid:: 1 1 2 2
    :gutter: 2 3 4 4

    .. grid-item-card::
        :img-top: ../_static/docker.svg
        :text-align: center

        **Docker image**
        ^^^

        The `geodels/gospl-examples <https://hub.docker.com/r/geodels/gospl-examples>`_
        image packages the full ``gospl-smoke`` conda environment (goSPL, PETSc,
        GMT, VTK, JupyterLab) and is the recommended way to run the examples
        without a local conda install.

        +++

        .. button-ref:: installDocker
            :color: secondary
            :click-parent:

            Install via Docker

    .. grid-item-card::
        :img-top: ../_static/anaconda.svg
        :text-align: center

        **Conda environment**
        ^^^

        Install goSPL from the pre-built ``geodels`` conda channel in one
        command, or build a full environment from ``environment.yml`` for
        editable/source installs.

        +++

        .. button-ref:: installConda
            :color: secondary
            :click-parent:

            Install via Conda


    .. grid-item-card::
        :img-top: ../_static/python.svg
        :text-align: center

        **Install from source**
        ^^^

        Build goSPL from the git source tree. Requires MPICH, PETSc ≥ 3.21,
        and the full scientific Python stack. Recommended for experienced users
        on Linux.

        +++

        .. button-ref:: installSrc
            :color: secondary
            :click-parent:

            Learn more

    .. grid-item-card::
        :img-top: ../_static/hpc.svg
        :text-align: center

        **Install on HPC**
        ^^^

        Deploy goSPL on NCI Gadi (PBS) or Pawsey Setonix (Slurm) — either
        as a virtual-environment native install or inside a Singularity/
        Apptainer container.

        +++

        .. button-ref:: installHPC
            :color: secondary
            :click-parent:

            Learn more

.. toctree::
    :maxdepth: 2
    :hidden:

    installDocker
    installConda
    installSrc
    installHPC
