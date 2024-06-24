.. _optfile1:


=================================
Climate related parameters
=================================
      
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