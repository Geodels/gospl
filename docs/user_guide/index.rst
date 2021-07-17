.. _user_guide:

================
User Guide
================

This user guide covers essential features of :mod:`gospl`, mostly in the form of interactive Jupyter notebooks and Python scripts. Reading this guide, you will learn:

- data structure used in gospl input file,
- how to generate initial conditions like topography, precipitation and tectonic maps to force a simulation,
- how to extract some of the output from the model results to visualise them in Jupyter notebooks,
- how to run sequence of backward/forward gospl models using Python functions,
- how to set a script for running gospl on HPC.

Notebooks cover just a small selection of functions as an illustration of principles. For a full overview of ``gospl`` capabilities, head to the `API reference <https://gospl.readthedocs.io/en/latest/api_ref/index.html>`_. For additional examples, you might be interested in the following set of examples available from the `Stellar-SFM project <https://geodels.github.io/stellar-sfm/welcome.html>`_.

Step 1 - The input file
------------------------------


.. toctree::
    :maxdepth: 3
    :hidden:

    inputfile


.. raw:: html

    <div class="container">
        <div class="row">
            <div class="col-12 d-flex install-block">
                <div class="card install-card shadow w-100">
                <div class="card-header">
                    Understanding the input file structure
                </div>
                <div class="card-body">
                  <p class="card-text"> Imposing initial conditions, specifying physical processes parameters and understanding how the input file is structured... </p>

.. container:: custom-button

    :ref:`Learn more <inputfile>`

.. raw:: html

                </div>
                </div>
            </div>

        </div>
    </div>

Step 2 - Tutorials via Jupyter notebooks
-------------------------------------------

Installing additional libraries & the examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here we assume that you have followed one of the methods described in the `Getting Started guide <https://gospl.readthedocs.io/en/latest/getting_started/index.html>`_ and have successfully installed :mod:`gospl` either via ``pip`` or ``conda``.

.. note::

  If you are using the docker environment then the additional libraries required to run the pre & post processing files are already installed as well as the notebooks examples and you can skip this step.

Pre/post processing libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are using ``conda``, you will first put your self inside this environment run::

    source activate gospl-package

On Windows the command is::

    activate gospl-package

We will now install some additional libraries. For ``conda``::

    conda install pyvista pyevtk panel netCDF4 gdown

or via ``pip``::

    pip install pyvista pyevtk panel netCDF4 gdown


Meshing libraries
~~~~~~~~~~~~~~~~~~~~~~~~

stripy
***************

`stripy <https://github.com/underworldcode/stripy>`_ is a Python interface to TRIPACK and STRIPACK Fortran code for (constrained) triangulation in Cartesian coordinates and on a sphere.


The library can be installed as a ``pip`` package::

    pip install stripy


JIGSAW
********

`JIGSAW Python <https://github.com/dengwirda/jigsaw-python>`_ is an unstructured mesh generator and tessellation library; designed to generate high-quality triangulations and polyhedral decompositions of general planar, surface and volumetric domains.

The library can be installed as a ``conda`` package::

    conda install jigsaw


Installing it from source is also relatively straightforward once you have a ``C++`` compiler and the ``cmake`` utility installed::

    # Clone/download + unpack this repository.
    git clone https://github.com/dengwirda/jigsaw-python.git
    # Install the library...
    cd jigsaw-python
    python3 setup.py build_external
    python3 setup.py install



Notebooks examples
~~~~~~~~~~~~~~~~~~~~~

The notebooks are available from the `Github repository <https://github.com/Geodels/gospl/tree/master/notebooks>`_ but you can also directly download them as a **tar** file from `here <https://drive.google.com/file/d/1SvRj27NBF4aA2E8svyniysQtDHuiVNIf/view?usp=sharing>`_ or using the following command in your terminal::

    gdown https://drive.google.com/uc?id=1SvRj27NBF4aA2E8svyniysQtDHuiVNIf
    tar xvf notebooks.tar


Running the paleo-constrained example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 3

    bfModel/bfModel


The above example is a simpler version (smaller temporal extent and coarse resolution) of the simulation presented :ref:`here <examples>`.

.. toctree::
    :maxdepth: 3
    :hidden:

    examples


Running the stratigraphic example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 3

    dualLith/dualLithology

Step 3 - Advanced workflows
-----------------------------------

This section does not provide any dataset but some of Jupyter notebooks, post-processing functions and scripts that one can use to run and analyse complex paleo-forcing global scale models.


.. raw:: html

    <div class="container">
        <div class="row">
            <div class="col-12 d-flex install-block">
                <div class="card install-card shadow w-100">
                <div class="card-header">
                    Workflows for backward / forward model on moving mesh
                </div>
                <div class="card-body">
                  <p class="card-text"> A set of scripts and proposed workflows to run
                  a model with plate motion and surface remeshing conditions. </p>

.. container:: custom-button

    :ref:`Look at the workflows <advance>`

.. raw:: html

                </div>
                </div>
            </div>

        </div>
    </div>



.. toctree::
    :maxdepth: 3

    advance

Step 4 - Setting up gospl for HPC
-----------------------------------

.. toctree::
    :maxdepth: 3

    HPC.rst
