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
                fDa: 10.
                fDm: 40.
                mthd = 1

        This part of the input file define the parameters for the fluvial surface processes based on the *Stream Power Law* (SPL) and is composed of:

        a. ``K`` representing the erodibility coefficient which is scale-dependent and its value depend on lithology and mean precipitation rate, channel width, flood frequency, channel hydraulics. It is used in the SPL law: :math:`E = K (\bar{P}A)^m S^n`

        .. warning::
            It is worth noting that the coefficient *n* is fixed and take the value *1*.

        b. Studies have shown that the physical strength of bedrock which varies with the degree of chemical weathering, increases systematically with local rainfall rate. Following `Murphy et al. (2016) <https://doi.org/10.1038/nature17449>`_, the stream power equation is adapted to explicitly incorporate the effect of local mean annual precipitation rate, P, on erodibility: :math:`E = (K_i P^d) (\bar{P}A)^m S^n`. ``d`` (:math:`d` in the equation) is a positive exponent that has been estimated from field-based relationships to 0.42. Its default value is set to 0.
        c. ``m`` is the coefficient from the SPL law: :math:`E = K (\bar{P}A)^m S^n` and takes the default value of 0.5.
        d. ``fDa`` dimensionless deposition coefficient for continental domain
        e. ``fDm`` dimensionless deposition coefficient for marine domain
        f. ``mthd`` chosen approach to account for sediment deposition (should be either 1 or 2)

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

        a. porosity at the surface ``phis``,
        b. e-folding depth ``z0s`` (in metres)

