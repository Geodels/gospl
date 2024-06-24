.. _inputfile:


.. important::

    The code is primarily a **parallel global scale landscape evolution model**, built to simulate **topography and basins** dynamics. The following processes are considered:

    - **river incision** and **deposition** using stream power law,
    - continental **deposition** in depressions,
    - **marine deposition** at river mouth,
    - **hillslope processes** in both marine and inland areas,
    - **sediment compaction** as stratigraphic layers geometry and properties change, 
    - spatially and temporally varying **tectonics** (horizontal and vertical displacements).
    - spatially and temporally varying **precipitation** grids as well as **orographic** rain and sea-level fluctuations, 
    - possibility to account for **flexural** isostasy driven by changes in surface loading.

==============================
Input file parameters
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


Stream Power Law parameters
---------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
                
        **Declaration example**:

        .. code:: python

            spl:
                K: 3.e-8
                d: 0.42
                m: 0.4
                fDa: 10.
                fDm: 40.
                mthd = 1

        This part of the input file define the parameters for the fluvial surface processes based on the *Stream Power Law* (SPL) and is composed of:

        a. ``K`` representing the erodibility coefficient which is scale-dependent and its value depend on lithology and mean precipitation rate, channel width, flood frequency, channel hydraulics. It is used in the SPL law: :math:`E = K (\bar{P}A)^m S^n`

        .. warning::
            It is worth noting that the coefficient *n* is fixed and take the value *1*.

        b. Studies have shown that the physical strength of bedrock which varies with the degree of chemical weathering, increases systematically with local rainfall rate. Following `Murphy et al. (2016) <https://doi.org/10.1038/nature17449>`_, the stream power equation is adapted to explicitly incorporate the effect of local mean annual precipitation rate, P, on erodibility: :math:`E = (K_i P^d) (\bar{P}A)^m S^n`. ``d`` (:math:`d` in the equation) is a positive exponent that has been estimated from field-based relationships to 0.42. Its default value is set to 0.
        c. ``m`` is the coefficient from the SPL law: :math:`E = K (\bar{P}A)^m S^n` and takes the default value of 0.5.
        d. ``fDa`` vg
        e. ``fDm`` vg
        f. ``mthd`` chosen approach 

Hillslope and marine deposition parameters
-------------------------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
                
        **Declaration example**:

        .. code:: python

            diffusion:
                hillslopeKa: 0.02
                hillslopeKm: 0.2
                smthDep: 20.0
                clinSlp: 5.e-5

        Hillslope processes in goSPL is defined using a classical *diffusion law* in which sediment deposition and erosion depend on slopes (*simple creep*). The following parameters can be tuned based on your model resolution:

        a. ``hillslopeKa`` is the diffusion coefficient for the aerial domain,
        b. ``hillslopeKm`` is the diffusion coefficient for the marine domain,
        c. ``smthDep`` is the transport coefficient of freshly deposited sediments entering the ocean from rivers,
        d. ``clinSlp`` is the maximum slope of clinoforms (needs to be positive), this slope is then used to estimate the top of the marine deposition based on distance to shore.        

Sea-level (eustatic) forcing
-----------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
                
        **Declaration example**:

        .. code:: python

            sea:
                position: 0.
                curve: 'data/sealevel.csv'


        The sea-level declaration is defined with 2 optional parameters:

        a. the relative sea-level ``position`` in meters (optional),
        b. a sea-level ``curve`` *e.g.* a file containing 2 columns (time and sea-level position).


Tectonic forcing parameters
----------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: python

            tectonic:
                - start: -20000000.
                end: -19000000.
                mapH: 'input8/disp20Ma'
                - start: -19000000.
                end: -18000000.
                mapH: 'input8/disp19Ma'
                - start: -18000000.
                end: -17000000.
                mapH: 'input8/disp18Ma'
                - start: -17000000.
                end: -16000000.
                mapH: 'input8/disp17Ma'
                mapV: 'input8/dispv17Ma'
                - start: -16000000.
                end: -15000000.
                mapV: 'input8/dispv16Ma'

        Follows the tectonic forcing conditions with a sequence of events defined by a starting time (``start``) and either a vertical only forcing (*e.g.* uplift and/or subsidence defined with ``mapV``) or a fully 3D displacement mesh ``mapH``. These displacements are set in metres per year.


.. important::

  As mentioned above and for the next key parameter as well, these forcing files are defined as numpy zip array (**.npz**).


Compaction & porosity variables defintion
------------------------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: python

            compaction:
                phis: 0.49
                z0s: 3700.0

        We assume  a depth-porosity relationship for the sediment compaction based on the following parameters:

        a. porosity at the surface ``phis``,
        b. e-folding depth ``z0s`` (in metres)


Climatic (rainfall) forcing conditions
----------------------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: python

            climate:
                - start: -20000000.
                map: ['input8/rain20Ma','r']
                - start: -15000000.
                uniform: 1.


        The climatic forcing is defined in a similar fashion as the tectonic one with again a sequence of events by a starting time (``start``) and either an uniform rainfall over the entire mesh (``uniform``) or with a precipitation mesh ``map``. The rainfall values have to be in metres per year.

Orographic rain definition
---------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: python

            orography:
                latitude: 40.0  
                wind_speed: 10.0 
                wind_dir: 0 
                nm: 0.005 
                env_lapse_rate: -4
                moist_lapse_rate: -7 
                ref_density: 7.4e-3 
                hw:  5000 
                conv_time: 1000. 
                fall_time: 1000. 
                oro_precip_base: 7.0 
                oro_precip_min: 0.01
                rainfall_frequency: 1 
            
        This part of the input file define the parameters for the orographic rain:

        a. ``latitude``: average latitude used to compute the Coriolis factors [degrees btw -90 and 90]; default 0
        b. ``wind_speed``: wind speed in m/s; default 10
        c. ``wind_dir``: wind direction [0: north, 270: west]; default 0
        d. ``nm``: moist stability frequency [1/s]; default 0.01
        e. ``env_lapse_rate``: environmental lapse rate [degrees Celsius/km]; default -4.0
        f. ``moist_lapse_rate``: moist adiabatic lapse rate [degrees Celsius/km]; default -7.0
        g. ``ref_density``: reference saturation water vapor density [kg/m^3]; default 7.4e-3
        h. ``hw``:  water vapor scale height [m]; default 3400
        i. ``conv_time``: cloud water to hydrometeor conversion time [s]; default 1000
        j. ``fall_time``: hydrometeor fallout time [s]; default 1000
        k. ``oro_precip_base``: non-orographic, uniform precipitation rate [mm/h]; default 7.
        l. ``oro_precip_min``: minimum precipitation [mm/h] when precipitation rate <= 0; default 0.01
        m. ``rainfall_frequency``: number of storm of 1 hour duration per day; default 1


Forcing paleo-topography definition
-----------------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: python

            forcepaleo:
                dir: 'output-backward'
                steps: [5,10,5]

        For simulations that require to be forced with paleo-topography maps obtained from backward models, you will also have to set this key composed of 2 parameters:

        a. ``dir`` the directory containing the outputs of the backward model,
        b. ``steps`` the steps from the model outputs that will be used to force the forward model topography.

.. important::

  The ``steps`` often correspond to the time where you have a paleotopography dataset that you want to match for example from a Scotese paleotopography map.

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
