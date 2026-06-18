.. _inputfile:

==============================
Required input file parameters
==============================

.. note::

  Input files for  goSPL are based on `YAML`_ syntax.
  The YAML structure is shown through indentation (one or more spaces) and sequence items are denoted by a dash.


Initial mesh definition and simulation declaration
---------------------------------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: yaml

            name: Description of your simulation run

            domain:
                npdata: ['input/mesh','v','c','z']
                flowdir: 5
                bc:'0101'
                fast: False
                seadepo: True
                nperodep: 'input/erodep20Ma'
                npstrata: 'input/sed20Ma'
                advect: 'iioe2'

        The following parameters are **required**:

        a. the initial spherical or 2D surface mesh ``npdata`` (**.npz** file). This file contains the following keys: ``v`` the mesh vertices coordinates, ``z`` the vertice elevations, ``c`` the mesh cells.
        b. the flow direction method to be used ``flowdir`` that takes an integer value between 1 (for SFD) and 6 (for MFD)
        
        In addition the following optional parameters could be set:

        c. boundary conditions (``bc``) when not running a global model: a 4-character string, one per edge in the order **north, east, south, west** (N at ``y = ymax``, E at ``x = xmax``, S at ``y = ymin``, W at ``x = xmin``). Each character is:

           - ``o`` — **open**: a deep base-level outlet (the edge is held below sea level); water and the sediment it carries drain out across it.
           - ``f`` — **fixed**: a base-level outlet at the natural edge elevation; flow and sediment still leave the domain there (this is *not* a no-flux wall).
           - ``w`` — **wall**: a true no-flux wall; flow is contained and sediment deposits against the edge instead of draining out.
           - ``c`` — **cyclic** (periodic): flow and sediment wrap from one edge to the opposite one.

           Legacy digits are accepted: ``0`` is read as ``o`` (open) and ``1`` as ``f`` (fixed). The default is ``'ffff'``. For example ``'ofof'`` opens north/south with fixed east/west; ``'wfwf'`` walls north/south with fixed east/west.

           Cyclic edges must be set as an **opposite pair** — both ``south`` and ``north``, *or* both ``east`` and ``west`` — and **at most one pair** may be cyclic (so up to two periodic edges, never all four). For example ``'ococ'`` makes east/west periodic with open north/south.

           .. important::

              A cyclic run **requires a periodic input mesh** in that direction: the mesh must be a **cylinder** (the periodic axis wrapped onto a circle whose circumference is the domain width), so that its cells genuinely connect the two seam edges. A cylinder is intrinsically flat, so its finite-volume geometry is identical to a periodic flat strip's; goSPL detects the cylinder's two open ends as the (non-periodic) boundary and routes flow/sediment across the seam through the wrapping cells. goSPL does **not** synthesise the wrap from an ordinary flat mesh (a planar wrap would have incorrect seam geometry). This works in parallel.

           .. note::

              The ``w`` (wall) boundary contains flow, and deposits the sediment that reaches the closed edge. A domain that is *fully* enclosed by walls (no ``o``/``f``/marine exit anywhere) conserves sediment over a single step, but **not** indefinitely once its basins fill: goSPL's spill-based sediment cascade cannot aggrade a fully-closed basin (it assumes excess always spills toward an outlet). Combine walls with at least one ``o``/``f`` edge (or a marine area) for a fully mass-conserving long run.
        d. the ``fast`` key allows you to run a model without applying any surface processes on top. This is used to check your input files prior to run your simulation with all options. By default it is set to *False*.
        e. ``seadepo`` performing marine deposition or not. By default it is set to *True*.
        f. to start a simulation using a previous erosion/deposition map use the ``nperodep`` key and specify a file (**.npz** format with the erosion deposition defined with the key ``ed``) containing for each vertex of the mesh the cumulative erosion deposition values in metres. 
        g. to start a simulation using an initial stratigraphic layer use the ``npstrata`` key (**.npz** file) and specify a file containing for each vertex of the mesh the stratigraphic layer thickness ``strataH``, the elevation at time of deposition ``strataZ``, and the porosities of the sediment ``phiS``. An optional ``stratK`` key (same ``(mesh_points, init_layers)`` shape) may be provided to impose a **per-layer erodibility multiplier**. When the layer is exposed at the surface, the effective SPL erodibility becomes ``K * stratK`` where ``K`` is the value declared in the ``spl`` block. Values ``< 1`` model softer-than-default bedrock; values ``> 1`` model more resistant lithologies; the default (key absent, or value ``1.0``) keeps the SPL behaviour unchanged. Sediment eroded from such a layer and re-deposited downstream loses the multiplier and reverts to ``1.0`` (i.e. the YAML-default ``K``).

        h. when dual lithology is enabled (``strata: dual: True``), two further optional keys (same ``(mesh_points, init_layers)`` shape) let each initial layer carry its **own coarse/fine composition**: ``strataHf`` is the fine-fraction bulk thickness of each layer (so the coarse thickness is ``strataH - strataHf``, and the per-layer fine fraction is ``strataHf / strataH``), and ``phiF`` is the per-layer fine porosity. ``strataHf`` is clamped to ``[0, strataH]`` on load; if it is absent the initial column is all-coarse, and if ``phiF`` is absent it defaults to the fine surface porosity ``phi0f``. These keys are ignored (with a warning) when dual lithology is off.
        h. ``advect`` define the advection scheme used when applying horizontal displacements. Choices are ``upwind``, ``iioe1``, ``iioe2`` and ``interp``  (go to the technical `information <https://gospl.readthedocs.io/en/latest/tech_guide/tecto.html#horizontal-advection>`_ in the documentation for more information). 

.. warning::

  It is worth noting that all the input files require to run a goSPL simulation must be defined as numpy zip array (**.npz**). This allows to directly and efficiently load the dataset during initialisation. This is specially efficient when running large models. 
  
.. note::

  It is also possible to have only one **.npz** file containing the required keys. For example you could set an input file containing the following keys: ``v`` the mesh vertices coordinates, ``z`` the vertice elevations, ``c`` the mesh cells, and ``ed`` the erosion deposition node values and call it for both the ``npdata`` and ``nperodep`` parameters.


Setting model temporal evolution
--------------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        

        **Declaration example**:

        .. code:: yaml

            time:
                start: -20000000.
                end: 0.
                dt: 250000.
                tout: 1000000.
                rstep: 25
                strat: 500000.

        ``time`` is also a required component of every input file. The following parameters are needed:

        a. ``start`` is the model start time in years,
        b. ``end`` is the model end time in years,
        c. ``dt`` is the model internal time step (the approach in goSPL uses an implicit time step).
        d. ``tout`` is the output interval used to create model outputs.

        The following parameters are **optional**:

        e. to restart a simulation use the ``rstep`` key and specify the time step number.
        f. ``strat`` is the stratigraphic timestep interval used to update the stratigraphic record.


.. important::

  In cases where the specify ``dt`` and ``strat`` parameters are greater than ``tout``, they will automatically be rescaled to match with the output interval. When turned-on the stratal records will be output at the same time as the output ones, but the file will potentially contain multiple stratigraphic layers per output if ``tout`` is a multiple of ``strat``.

Output folder definition
-------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: yaml

            output:
                dir: 'forward'
                makedir: False

        Finally, you will need to specify the output folder, with 2 possible parameters:

        a. ``dir`` gives the output directory name.

        The following parameter is **optional**:

        b. the option ``makedir`` gives the ability to delete any existing output folder with the same name (if set to *False* - default value) or to create a new folder with the given `dir` name plus a number at the end (*e.g.* outputDir_XX if set to *True* with XX the run number). It allows you to avoid overwriting on top of previous runs.

.. _`YAML`: https://circleci.com/blog/what-is-yaml-a-beginner-s-guide/
