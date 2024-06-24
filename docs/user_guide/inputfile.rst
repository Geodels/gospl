.. _inputfile:

==============================
Required input file parameters
==============================

.. note::

  Input files for  goSPL are based on `YAML`_ syntax.
  The YAML structure is shown through indentation (one or more spaces) and sequence items are denoted by a dash. At the moment the following component are available:


Initial mesh definition and simulation declaration
---------------------------------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: python

            name: Global model from 20 Ma to present

            domain:
                npdata: 'input8/elev20Ma'
                flowdir: 5
                flowexp: 0.5
                bc:'0101'
                fast: False
                seadepo: True
                nperodep: 'strat8/erodep20Ma'
                npstrata: 'strat8/sed20Ma'

        The following parameters are **required**:

        a. the initial spherical surface mesh ``npdata`` as well as
        b. the flow direction method to be used ``flowdir`` that takes an integer value between 1 (for SFD) and 6 (for MFD)
        c. the exponent (``flowexp``) used in the flow direction approach. Default value is set to 0.42.
        d. boundary conditions (``bc``) when not running a global model. Each integer corresponds to an edge defined in the following order: south, east, north, and west. The integer is set to 0 for open and 1 closed boundaries.

        In addition the following optional parameters can be set:

        e. the ``fast`` key allows you to run a model without applying any surface processes on top. This is used to run backward model in a quick way, but can also potential be set to *True* if you want to check your input files prior to running a forward model with all options.
        f. ``seadepo`` performing marine deposition or not
        g. to start a simulation using a previous erosion/deposition map use the ``nperodep`` key and specify a file containing for each vertex of the mesh the cumulative erosion deposition values in metres.
        h. to start a simulation using an initial stratigraphic layer use the ``npstrata`` key and specify a file containing for each vertex of the mesh the stratigraphic layer thickness and the porosities of the sediments.

.. warning::

  It is worth noting that all the input files require to run a goSPL simulation must be defined as numpy zip array (**.npz**). This allows to directly and efficiently load the dataset during initialisation. This is specially efficient when running large models.


Setting model temporal evolution
--------------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        

        **Declaration example**:

        .. code:: python

            time:
                start: -20000000.
                end: 0.
                dt: 250000.
                tout: 1000000.
                rstep: 25
                tec: 1000000.
                strat: 500000.

        ``time`` is also a required component of every input file. The following parameters are needed:

        a. ``start`` is the model start time in years,
        b. ``end``` is the model end time in years,
        c. ``dt`` is the model internal time step (the approach in goSPL uses an implicit time step).
        d. ``tout`` is the output interval used to create model outputs,
        e. to restart a simulation use the ``rstep`` key and specify the time step number.
        f. ``tec`` is the tectonic timestep interval used to update the tectonic meshes and perform the required horizontal displacements (vertical displacements are done every ``dt``).
        g. ``strat`` is the stratigraphic timestep interval used to update the stratigraphic record.


.. important::

  In cases where the specify ``dt``, ``strat`` and ``tec`` parameters are greater than ``tout``, they will automatically be rescaled to match with the output interval. The ``tec`` parameter should be set to similar to the temporal time step used in your reconstruction (usually around 1Ma). This time step is used to perform the horizontal displacements. The vertical displacements are updated for each time step. When turn-on the stratal records will be output at the same time as the output ones, but the file will potentially contain multiple stratigraphic layers per output if ``strat`` is lower than ``tout``.

Output folder definition
-------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: python

            output:
                dir: 'forward'
                makedir: False

        Finally, you will need to specify the output folder, with 2 possible parameters:

        a. ``dir`` gives the output directory name and
        b. the option ``makedir`` gives the ability to delete any existing output folder with the same name (if set to False) or to create a new folder with the given `dir` name plus a number at the end (*e.g.* outputDir_XX if set to True with XX the run number). It allows you to avoid overwriting on top of previous runs.

.. _`Paraview`: https://www.paraview.org/download/
.. _`YAML`: https://circleci.com/blog/what-is-yaml-a-beginner-s-guide/
