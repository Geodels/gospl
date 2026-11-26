.. _user_guide:

================
User Guide
================

This user guide explains the different parameters available in goSPL input file.

.. note::

    **Start here: which modules do I need?** Choose a regional or a global model and the
    processes you want to simulate: the chooser shows what each one needs, the closest
    example, the main parameters with their usual values, an estimate of the time and memory
    of your run, and writes a skeleton input file you can download. The units, defaults and
    usual ranges of every input key, and the meaning and units of every output field, are on
    the two reference pages.

    .. grid:: 1 1 3 3
        :gutter: 2

        .. grid-item::

            .. button-ref:: module_chooser
                :ref-type: ref
                :color: primary
                :expand:

                Which modules do I need?

        .. grid-item::

            .. button-ref:: input_reference
                :ref-type: ref
                :color: primary
                :outline:
                :expand:

                Input parameters: units and ranges

        .. grid-item::

            .. button-ref:: output_reference
                :ref-type: ref
                :color: primary
                :outline:
                :expand:

                Output fields: units and ranges

    Ready-to-run models, regional and global, are in the :ref:`examples gallery <examples>`
    (from the goSPL-examples :examples:`repository <.>`): how to build the inputs (meshes,
    rain, tectonic maps), run a model and look at the results. Every function of the code is
    described in the :ref:`API reference <api_ref>`.


Input file
---------------------

.. grid:: 1
    :padding: 2

    .. grid-item-card::  
        
        The code is primarily a **parallel global scale landscape evolution model**, built to simulate **topography and basins** dynamics. The following processes are considered:

        - **river incision** and **deposition** using stream power law (detachment- & transport-limited),
        - continental **deposition** in depressions,
        - **marine deposition** at river mouth,
        - **hillslope processes** (linear and non-linear) in both marine and inland areas,
        - **soil / regolith** production, weathering and transport,
        - **groundwater** (implicit water table) driving a generic **duricrust** that armors erodibility, with optional **conservative solute geochemistry** (multi-tracer dissolution → transport → precipitation → river/ocean export),
        - **ice sheets and glacial erosion** (glacial abrasion, till/moraine deposition, meltwater and ice loading),
        - **dual-lithology** (coarse/fine) sediment tracking and **sediment provenance** (source-to-sink attribution),
        - **sediment compaction** as stratigraphic layers geometry and properties change,
        - spatially and temporally varying **tectonics** (horizontal and vertical displacements),
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


Groundwater, duricrust and geochemistry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. grid:: 1
    :padding: 3

    .. grid-item-card::
        :text-align: center

        **Water table, duricrust and solute geochemistry**
        ^^^

        An implicit water table, a crust that hardens the surface, and the
        dissolution, transport and export of solutes down the rivers to the
        ocean.

        +++

        .. button-ref:: gw_user
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


Running goSPL
^^^^^^^^^^^^^

.. grid:: 1
    :padding: 3

    .. grid-item-card::
        :text-align: center

        **Launching a simulation**
        ^^^

        The ``gospl`` command (``mpirun -np 8 gospl -i input.yml``) and the
        equivalent Python ``Model`` API.

        +++

        .. button-ref:: running
            :color: secondary
            :click-parent:

            Learn more.


Model outputs
^^^^^^^^^^^^^

.. grid:: 1
    :padding: 3

    .. grid-item-card::
        :text-align: center

        **Output variables & visualisation**
        ^^^

        The fields goSPL writes, what they mean, and how to view them in ParaView.

        +++

        .. button-ref:: outputs
            :color: secondary
            :click-parent:

            Learn more.


Performance on HPC
^^^^^^^^^^^^^^^^^^^

.. grid:: 1
    :padding: 3

    .. grid-item-card::
        :text-align: center

        **Strong-scaling benchmarks**
        ^^^

        Speedup and parallel efficiency on global meshes (10, 8, 5 km).

        +++

        .. button-ref:: hpc_scaling
            :color: secondary
            :click-parent:

            Learn more.

.. toctree::
    :maxdepth: 3
    :hidden:

    reference/modules
    inputfile
    surfproc
    groundwater
    optfile1
    optfile2
    running
    outputs
    reference/inputs
    reference/outputs_ref
    hpc_scaling
