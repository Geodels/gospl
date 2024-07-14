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

        .. code:: python

            spl:
                K: 3.e-8
                d: 0.42
                m: 0.4
                G: 1.

        This part of the input file define the parameters for the fluvial surface processes based on the *Stream Power Law* (SPL) and is composed of:

        a. ``K`` representing the erodibility coefficient which is scale-dependent and its value depend on lithology and mean precipitation rate, channel width, flood frequency, channel hydraulics. It is used in the SPL law: :math:`E = K (\bar{P}A)^m S^n`

        The following parameters are **optional**:

        b. Studies have shown that the physical strength of bedrock which varies with the degree of chemical weathering, increases systematically with local rainfall rate. Following `Murphy et al. (2016) <https://doi.org/10.1038/nature17449>`_, the stream power equation could be adapted to explicitly incorporate the effect of local mean annual precipitation rate, P, on erodibility: :math:`E = (K_i P^d) (\bar{P}A)^m S^n`. ``d`` (:math:`d` in the equation) is a positive exponent that has been estimated from field-based relationships to 0.42. Its default value is set to 0.0
        c. ``m`` is the coefficient from the SPL law: :math:`E = K (\bar{P}A)^m S^n` and takes the default value of 0.5.

        .. note::
            It is worth noting that the coefficient *n* in the SPL is fixed and take the value *1*.

        d. ``G`` dimensionless deposition coefficient for continental domain when accounting for sedimentation rate in the SPL following the model of `Yuan et al, 2019 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JF004867>`_. The default value is 0.0 (purely detachment-limited model).
        

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
                nonlinKm: 100.0
                clinSlp: 5.e-5

        Hillslope processes in goSPL is defined using a classical *diffusion law* in which sediment deposition and erosion depend on slopes (*simple creep*). The following parameters can be tuned based on your model resolution:

        a. ``hillslopeKa`` is the diffusion coefficient for the aerial domain,
        b. ``hillslopeKm`` is the diffusion coefficient for the marine domain,
        c. ``nonlinKm`` is the transport coefficient of freshly deposited sediments entering the ocean from rivers (non-linear diffusion),
        d. ``clinSlp`` is the maximum slope of clinoforms (needs to be positive), this slope is then used to estimate the top of the marine deposition based on distance to shore.       

Sediment surface erodibility factor
-------------------------------------


.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: python

            sedfactor:
                - start: 200000.
                  uniform: 3
                - start: 400000.
                  map: ['facErosion4k','fsed']

        One could choose to impose variable erodibility factors through space and time to reflect different surficial rock composition. For example, those maps could be set to represent different rock erodibility index as proposed in `Mossdorf et al. (2018) <https://www.sciencedirect.com/science/article/abs/pii/S0143622817306859>`_. The factor are then used in front of the erodibility coefficient (``K`` in the SPL).

        .. important::

            When defining your variable erodibility factors grid, you needs to use the **npz** format and your factors would be specified by a key corresponding to the factor values for each vertice of the mesh. In the above example this key is ``'fsed'``. 


Compaction & porosity variables definition
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

        a. porosity at the surface ``phis``, default value is set to 0.49,       
        b. e-folding depth ``z0s`` (in metres), default value is set to 3700.       

