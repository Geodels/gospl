.. _surfproc:

=================================
Surface processes  parameters
=================================

Stream Power Law parameters
---------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
                
        **Declaration example**:

        .. code:: yaml

            spl:
                K: 3.e-8
                d: 0.42
                m: 0.4
                n: 1.0
                G: 1.

        This part of the input file define the parameters for the fluvial surface processes based on the *Stream Power Law* (SPL) and is composed of:

        a. ``K`` representing the erodibility coefficient which is scale-dependent and its value depend on lithology and mean precipitation rate, channel width, flood frequency, channel hydraulics. It is used in the SPL law: :math:`E = K (\bar{P}A)^m S^n`

        The following parameters are **optional**:

        b. Studies have shown that the physical strength of bedrock which varies with the degree of chemical weathering, increases systematically with local rainfall rate. Following `Murphy et al. (2016) <https://doi.org/10.1038/nature17449>`_, the stream power equation could be adapted to explicitly incorporate the effect of local mean annual precipitation rate, P, on erodibility: :math:`E = (K_i P^d) (\bar{P}A)^m S^n`. ``d`` (:math:`d` in the equation) is a positive exponent that has been estimated from field-based relationships to 0.42. Its default value is set to 0.0
        c. ``m`` is the flow accumulation coefficient from the SPL law: :math:`E = K (\bar{P}A)^m S^n` and takes the default value of 0.5.
        d. ``n`` is the slope coefficient from the SPL law: :math:`E = K (\bar{P}A)^m S^n` and takes the default value of 1.0.
        e. ``G`` dimensionless deposition coefficient for continental domain when accounting for sedimentation rate in the SPL following the model of `Yuan et al, 2019 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JF004867>`_. The default value is 0.0 (purely detachment-limited model).
        

Hillslope and marine deposition parameters
-------------------------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
                
        **Declaration example**:

        .. code:: yaml

            diffusion:
                hillslopeKa: 0.02
                hillslopeKm: 0.2
                nonlinKm: 500.0
                clinSlp: 5.e-5

        Hillslope processes in goSPL is defined using a classical *diffusion law* in which sediment deposition and erosion depend on slopes (*simple creep*). The marine deposition of freshly deposited sediments by rivers is obtained using a non-linear diffusion and the following parameters can be tuned based on your model resolution:

        a. ``hillslopeKa`` is the diffusion coefficient for the aerial domain,
        b. ``hillslopeKm`` is the diffusion coefficient for the marine domain,
        c. ``nonlinKm`` is the transport coefficient of freshly deposited sediments entering the ocean from rivers (non-linear diffusion),
        d. ``clinSlp`` is the maximum slope of clinoforms (needs to be positive), this slope is then used to estimate the top of the marine deposition based on distance to shore.        
                
        *Optional additions for non-linear diffusion model*

        A more complex version of the creep law involves a non-linear relationship between soil flux and topographic gradient. 
        
        .. note::
    
            Several non-linear creep-transport laws have been suggested in the literature and 2 non-linear formulations are available in goSPL. 

        **Either** by adding the following to the above parameters:

        .. code:: yaml

                hillslopenl: 2.5
                
        e. ``hillslopenl`` is the slope exponent in the non-critical hillslope model defined in the work of `Wang et al. (2024) <https://www.sciencedirect.com/science/article/pii/S0169555X24001053>`_. Here the coefficient of diffusion is set to the values of ``hillslopeKa`` and ``hillslopeKm``.

        **Or** by defining the following two parameters:

        .. code:: yaml

                hillslopeSc: 0.8
                hillslopeNb: 4

        f. ``hillslopeSc`` is the critical slope,
        g. ``hillslopeNb`` is the number of terms used in the truncated Taylor series formulation.

        In this last model, the non-linear creep formulation is described in `Barnhart et al. (2019) <https://gmd.copernicus.org/articles/12/1267/2019/gmd-12-1267-2019.pdf>`_ (section 3.4.3 **EQ. 14**).


Soil production, erosion, transport and deposition
-----------------------------------------------------


.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: yaml

            soil:
                soilK: 4.e-6
                maxProd: 50.e-6
                depthProd: 0.5
                roughnessL: 1.0
                decayDepth: 0.7
                bedrockConv: 0.0001
                uniform: 0.5
                map: ['test_mesh8/hsoil', 'soil']

        a. ``soilK`` is the erodibility coefficient for soil,
        b. ``maxProd`` is the soil production maximum rate (m/yr),
        c. ``depthProd`` is the soil production decay depth (m),
        d. ``roughnessL`` is the roughness length scale,
        e. ``decayDepth`` is the soil transport decay depth for non-linear diffusion,
        f. ``bedrockConv`` is the soil to bedrock conversion fraction, bedrock begins where soil production is a very small fraction of the maximum soil production (optional). 

        Then the user can specify the initial soil thickness if any by setting **either**:

        g. ``uniform`` a uniform soil thickness on the entire surface (m),

        **or**:

        h. ``map`` a soil thickness map. 

        .. important::

            When defining a soil thickness grid, one needs to use the **npz** format and needs to specify the key corresponding to the soil thickness value in the file. In the above example this key is ``'soil'``. The soil grid needs to define values for all vertices in the mesh in metres.
            

Sediment surface erodibility factor
-------------------------------------


.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: yaml

            sedfactor:
                - start: 200000.
                  uniform: 3
                - start: 400000.
                  map: ['facEro','fsed']

        One could choose to impose variable erodibility factors through space and time to reflect different surficial rock composition. For example, those maps could be set to represent different rock erodibility index as proposed in `Mossdorf et al. (2018) <https://www.sciencedirect.com/science/article/abs/pii/S0143622817306859>`_. The factor are then used in front of the erodibility coefficient (``K`` in the SPL).

        .. important::

            When defining your variable erodibility factors grid, you needs to use the **npz** format and your factors would be specified by a key corresponding to the factor values for each vertice of the mesh. In the above example this key is ``'fsed'``. 


Compaction & porosity variables definition
------------------------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: yaml

            compaction:
                phis: 0.49
                z0s: 3700.0

        We assume  a depth-porosity relationship for the sediment compaction based on the following parameters:

        a. porosity at the surface ``phis``, default value is set to 0.49,       
        b. e-folding depth ``z0s`` (in metres), default value is set to 3700.       

        .. note::
    
            See the technical `documentation <https://gospl.readthedocs.io/en/latest/tech_guide/strat.html>`_ for more information.
