.. _optfile1:


=================================
Climate related parameters
=================================
      
Sea-level (eustatic) forcing
-----------------------------

.. note::
    
    By default, the sea-level position in goSPL is set to 0 m. If you wish to set it to another position you can use the ``position`` key that changes the sea-level to a new value relative to sea-level. Another option consists in defining your own sea-level curve (``curve``) or using a published one (e.g. Haq curve for example). 

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
                
        **Declaration example**:

        .. code:: python

            sea:
                position: 0.
                curve: 'sealevel.csv'


        The sea-level declaration is defined with 2 optional parameters:

        a. the relative sea-level ``position`` in meters (optional), its default value is set to 0.0
        b. a sea-level ``curve`` *e.g.* a file containing 2 columns (time and sea-level position). Not required in case no sea-level fluctuations needs to be specified. 

.. important::

    The sea-level curve is defined as a 2 columns **csv** file containing in the first column the time in years (it doesn't need to be regularly temporally spaced) and in the second the sea-level position for the given time. When goSPL interprets this file, it will interpolate linearly between the defined times to find the position of the sea-level for every time step.

Climatic (rainfall) forcing conditions
----------------------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: python

            climate:
                - start: 1000.
                  map: ['rain20Ma','r']
                - start: 20000.
                  uniform: 1.


        The climatic forcing is defined in a similar fashion as the tectonic one. It is based on a sequence of events with each event starting at a given time (``start`` in years) and corresponding to a given precipitation condition. This could either be an uniform rainfall over the entire mesh (``uniform``) or a precipitation mesh ``map``. The rainfall values have to be in metres per year and the precipitation is updated at every time step (defined by ``dt``).

.. important::

    When defining a precipitation grid, one needs to use the **npz** format and needs to specify the key corresponding to the precipitation variable in the file. In the above example this key is ``'r'``. The precipitation grid needs to define values for all vertices in the mesh.


Orographic precipitation definition
------------------------------------

goSPL implements the Linear Theory of Orographic Precipitation following `Smith & Barstad (2004) <https://journals.ametsoc.org/view/journals/atsc/61/12/1520-0469_2004_061_1377_altoop_2.0.co_2.xml>`_.

.. note::
    
    The model includes airflow dynamics, condensed water advection, and downslope evaporation. It consists of two vertically-integrated steady-state advection equations describing: 

    1. the cloud water density and 
    2. the hydrometeor density. 

    Solving these equations using Fourier transform techniques, derives a single formula relating terrain and precipitation.

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: python

            orography:
                regdx: 200.
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

        a. ``regdx``: the resolution of the regular grid used to perform the orographic rain calculation.
        b. ``latitude``: average latitude used to compute the Coriolis factors [degrees btw -90 and 90]; default 0
        c. ``wind_speed``: wind speed in m/s; default 10
        d. ``wind_dir``: wind direction [0: north, 270: west]; default 0
        e. ``nm``: moist stability frequency [1/s]; default 0.01
        f. ``env_lapse_rate``: environmental lapse rate [degrees Celsius/km]; default -4.0
        g. ``moist_lapse_rate``: moist adiabatic lapse rate [degrees Celsius/km]; default -7.0
        h. ``ref_density``: reference saturation water vapor density [kg/m^3]; default 7.4e-3
        i. ``hw``:  water vapor scale height [m]; default 3400
        j. ``conv_time``: cloud water to hydrometeor conversion time [s]; default 1000
        k. ``fall_time``: hydrometeor fallout time [s]; default 1000
        l. ``oro_precip_base``: non-orographic, uniform precipitation rate [mm/h]; default 7.
        m. ``oro_precip_min``: minimum precipitation [mm/h] when precipitation rate <= 0; default 0.01
        n. ``rainfall_frequency``: number of storm of 1 hour duration per day; default 1

.. warning::

    In case where **flexure** and **orographic rain** capabilities are defined in the same simulation, you will need to have the same grid resolution (``regdx``) for each definition.