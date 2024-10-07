.. _api_ref:

==============
API reference
==============


.. note::
    
    This section provides an overview of goSPL main objects, functions and methods.

.. grid:: 1 1 2 2
    :padding: 3
    :gutter: 2 3 4 4

    .. grid-item-card::
        :text-align: center

        **Class Model**
        ^^^

        Instantiates model object and performs surface processes evolution. 

        +++

        .. button-ref:: model_ref
            :color: secondary
            :click-parent:

            See functions and source code.


    .. grid-item-card::
        :text-align: center

        **Class Mesh**
        ^^^

        Defines spherical mesh characteristics and builds PETSc DMPlex.

        +++

        .. button-ref:: mesh_ref
            :color: secondary
            :click-parent:

            See functions and source code.


    .. grid-item-card::
        :text-align: center

        **Class Pit**
        ^^^

        Performs parallel depression filling of the surface.

        +++

        .. button-ref:: pit_ref
            :color: secondary
            :click-parent:

            See functions and source code.

    .. grid-item-card::
        :text-align: center

        **Class Flow**
        ^^^

        Flow accumulation computation for unstructured mesh.

        +++

        .. button-ref:: flow_ref
            :color: secondary
            :click-parent:

            See functions and source code.

    .. grid-item-card::
        :text-align: center

        **Class SPL**
        ^^^

        Stream power law linear slope-dependency.

        +++

        .. button-ref:: spl_ref
            :color: secondary
            :click-parent:

            See functions and source code.

    .. grid-item-card::
        :text-align: center

        **Class SPL (non-linear)**
        ^^^

        Stream power law non-linear slope-dependency.

        +++

        .. button-ref:: nlspl_ref
            :color: secondary
            :click-parent:

            See functions and source code.

    .. grid-item-card::
        :text-align: center

        **Class SPL with soil**
        ^^^

        Stream power law accounting for soil production.

        +++

        .. button-ref:: soilspl_ref
            :color: secondary
            :click-parent:

            See functions and source code.

    .. grid-item-card::
        :text-align: center

        **Class Sediment Continent**
        ^^^

        Functions related to sediment transport and deposition for continental regions. 

        +++

        .. button-ref:: sed_ref
            :color: secondary
            :click-parent:

            See functions and source code.


    .. grid-item-card::
        :text-align: center

        **Class Sediment Marine**
        ^^^

        Functions related to sediment transport and deposition for marine continental regions. 

        +++

        .. button-ref:: sea_ref
            :color: secondary
            :click-parent:

            See functions and source code.
    
    .. grid-item-card::
        :text-align: center

        **Class Hillslope**
        ^^^

        Functions to evaluate linear and non-linear hillslope processes.

        +++

        .. button-ref:: hill_ref
            :color: secondary
            :click-parent:

            See functions and source code.
    
    .. grid-item-card::
        :text-align: center

        **Class Tectonics**
        ^^^

        Functions to evaluate vertical and horizontal interpolation when considering horizontal displacements.

        +++

        .. button-ref:: plate_ref
            :color: secondary
            :click-parent:

            See functions and source code.

    .. grid-item-card::
        :text-align: center

        **Class Stratigraphy**
        ^^^

        Functions related to stratigraphic architecture and compaction.

        +++

        .. button-ref:: stra_ref
            :color: secondary
            :click-parent:

            See functions and source code.

    .. grid-item-card::
        :text-align: center

        **Class IceMesh**
        ^^^

        Functions related to glacial erosion.

        +++

        .. button-ref:: ice_ref
            :color: secondary
            :click-parent:

            See functions and source code.

    .. grid-item-card::
        :text-align: center

        **Class Input**
        ^^^

        Input methods declaration.

        +++

        .. button-ref:: in_ref
            :color: secondary
            :click-parent:

            See functions and source code.


    .. grid-item-card::
        :text-align: center

        **Class Grid Processes**
        ^^^

        Functions related to additional processes performed on a regular grid.

        +++

        .. button-ref:: grid_ref
            :color: secondary
            :click-parent:

            See functions and source code.

    .. grid-item-card::
        :text-align: center

        **Class Output**
        ^^^

        Output methods declaration.

        +++

        .. button-ref:: out_ref
            :color: secondary
            :click-parent:

            See functions and source code.

.. toctree::
    :maxdepth: 3
    :hidden:

    model_ref
    mesh_ref
    pit_ref
    flow_ref
    spl_ref
    nlspl_ref
    soilspl_ref
    sed_ref
    sea_ref
    hill_ref
    plate_ref
    ice_ref
    stra_ref
    in_ref
    grid_ref
    out_ref
