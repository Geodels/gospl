.. _tech_guide:

================
Technical Guide
================

This guide covers the implicit, iterative approaches used to solve the multiple flow direction water routing and the erosion deposition processes main algorithms implemented in goSPL.

goSPL is mostly written in ``Python`` with some functions in ``Fortran`` and takes advantage of ``PETSc`` solvers over parallel computing architectures using ``MPI``.

Further information on any specific methods can be obtained in the :ref:`api_ref`.

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Short Description**
        ^^^

        The code is primarily a **parallel global scale landscape evolution model**, built to simulate **topography and basins** dynamics. The following processes are considered:

        - **river incision** and **deposition** using stream power law,
        - continental **deposition** in depressions,
        - **marine deposition** at river mouth,
        - **hillslope processes** in both marine and inland areas,
        - **sediment compaction** as stratigraphic layers geometry and properties change, 
        - spatially and temporally varying **tectonics** (horizontal and vertical displacements).
        - spatially and temporally varying **precipitation** grids as well as **orographic** rain and sea-level fluctuations, 
        - possibility to account for **flexural** isostasy driven by changes in surface loading.

.. grid:: 1 1 2 2
    :gutter: 2 3 4 4

    .. grid-item-card::
        :img-top: ../_static/rivers.png
        :text-align: center

        **River flow implementation**
        ^^^

        Based on a parallel implicit drainage area (IDA) method. Want to gain insights on the implemented approach?

        +++

        .. button-ref:: flow
            :color: secondary
            :click-parent:

            Learn more about goSPL flow.

    .. grid-item-card::
        :img-top: ../_static/erosion.png
        :text-align: center

        **Erosion rate and sediment flux**
        ^^^

        Based on the stream power law (SPL), river erosion depends on local slope,discharge and erodibility coefficient.

        +++

        .. button-ref:: ero
            :color: secondary
            :click-parent:

            Learn more about the SPL.

    .. grid-item-card::
        :img-top: ../_static/depression.png
        :text-align: center

        **Inland depressions & deposition**
        ^^^

        Computes the evolution in internally drained basins using a priority-flood algorithm.

        +++

        .. button-ref:: dep
            :color: secondary
            :click-parent:

            Learn more about inland depressions.

    .. grid-item-card::
        :img-top: ../_static/deltaDiagram.jpg
        :text-align: center

        **Hillslope and marine deposition**
        ^^^

        Change in elevation induced by creep law and transport of river sediment in the marine realm based on diffusion equations.

        +++

        .. button-ref:: hill
            :color: secondary
            :click-parent:

            Learn more about hillslope processes.

    .. grid-item-card::
        :img-top: ../_static/strati.png
        :text-align: center

        **Stratigraphy and compaction**
        ^^^

        Record stratigraphic layers through time, track two sediment types and compute porosity changes induced by deposition.

        +++

        .. button-ref:: strat
            :color: secondary
            :click-parent:

            Stratigraphic and compaction implementation.

    .. grid-item-card::
        :img-top: ../_static/tectonic.png
        :text-align: center

        **Tectonic forcing**
        ^^^

        Displacements either induced by lithospheric or mantle forcing are used to move the surface both horizontally and vertically.

        +++

        .. button-ref:: tecto
            :color: secondary
            :click-parent:

            Learn more about the tectonics implementation.


.. toctree::
    :maxdepth: 3
    :hidden:

    flow
    ero
    dep
    hill
    strat
    tecto
