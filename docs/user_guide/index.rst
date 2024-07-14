.. _user_guide:

================
User Guide
================

This user guide explains the different parameters available in goSPL input file.

.. note::
    For examples using goSPL, the user is invited to download the **goSPL-examples** `repository <https://github.com/Geodels/goSPL-examples>`_ which covers some of the basic functionalities of the code: 

    - data structure used in the gospl input file,
    - how to generate initial conditions like topography, precipitation and tectonic maps to force a simulation,
    - how to extract some of the output from the your results.

    Those examples highlight just a small selection of functions as an illustration of principles. 

    For a full overview of goSPL capabilities, head to the `API reference <https://gospl.readthedocs.io/en/latest/api_ref/index.html>`_. 

    For additional examples, you might be interested in the following set of examples available from the `Stellar-SFM project <https://geodels.github.io/stellar-sfm/welcome.html>`_.

    .. warning::

        `Stellar-SFM project <https://geodels.github.io/stellar-sfm/welcome.html>`_ is based on a previous version of goSPL and some of the new features will not be available. There might also be some of the oldest features that might not be fully functional with the newest version. If you are interested in reproducing the examples from this Stellar-SFM project, you could use the branch `v2023 <https://github.com/Geodels/gospl/tree/v2023>`_ from goSPL. 


Input file
---------------------

The code is primarily a **parallel global scale landscape evolution model**, built to simulate **topography and basins** dynamics. The following processes are considered:

- **river incision** and **deposition** using stream power law,
- continental **deposition** in depressions,
- **marine deposition** at river mouth,
- **hillslope processes** in both marine and inland areas,
- **sediment compaction** as stratigraphic layers geometry and properties change, 
- spatially and temporally varying **tectonics** (horizontal and vertical displacements).
- spatially and temporally varying **precipitation** grids as well as **orographic** rain and sea-level fluctuations, 
- possibility to account for **flexural** isostasy driven by changes in surface loading.


Required parameters
^^^^^^^^^^^^^^^^^^^^

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        :text-align: center
        
        **Mesh and temporal definition**
        ^^^

        Imposing initial nesh conditions and simulation duration.

        +++

        .. button-ref:: inputfile
            :color: secondary
            :click-parent:

            Learn more.


Parameters related to surface processes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        :text-align: center
        
        **Surface processes parameters**
        ^^^

        Definition of the stream power law and hillslope parameters.

        +++

        .. button-ref:: surfproc
            :color: secondary
            :click-parent:

            Learn more.


Optional parameters related to climate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        :text-align: center
        
        **Climate-related conditions**
        ^^^

        Imposing precipitation and sea-level.

        +++

        .. button-ref:: optfile1
            :color: secondary
            :click-parent:

            Learn more.


Optional parameters related to tectonics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        :text-align: center
        
        **Adding tectonic forcing**
        ^^^

        Imposing tectonic conditions on your mesh.

        +++

        .. button-ref:: optfile2
            :color: secondary
            :click-parent:

            Learn more.

.. toctree::
    :maxdepth: 3
    :hidden:

    inputfile
    surfproc
    optfile1
    optfile2
