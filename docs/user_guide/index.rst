.. _user_guide:

================
User Guide
================

This user guide covers essential features of goSPL, mostly in the form of interactive Jupyter notebooks and Python scripts. Reading this guide, you will learn:

- data structure used in the gospl input file,
- how to generate initial conditions like topography, precipitation and tectonic maps to force a simulation,
- how to extract some of the output from the model results to visualise them in Jupyter notebooks,
- how to run a sequence of backward/forward gospl models using Python functions,
- how to set a script for running gospl on HPC.

Notebooks cover just a small selection of functions as an illustration of principles. For a full overview of goSPL capabilities, head to the `API reference <https://gospl.readthedocs.io/en/latest/api_ref/index.html>`_. For additional examples, you might be interested in the following set of examples available from the `Stellar-SFM project <https://geodels.github.io/stellar-sfm/welcome.html>`_.


Step 1 - The input file
------------------------------

.. toctree::
    :maxdepth: 3
    :hidden:

    inputfile


.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        :text-align: center
        
        **Understanding the input file structure**
        ^^^

        Imposing initial conditions, specifying physical processes parameters and understanding how the input file is structured...

        +++

        .. button-ref:: inputfile
            :color: secondary
            :click-parent:

            Learn more.

