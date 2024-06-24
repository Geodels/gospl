.. _optfile2:


==============================
Tectonic related parameters
==============================

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

        It defines the tectonic forcing conditions from a sequence of events defined by a starting and ending time (``start`` and ``end``) and either a vertical only forcing (*e.g.* uplift and/or subsidence defined with ``mapV``) or a fully 3D displacement mesh ``mapH``. **These displacement rates are set in metres per year**.


.. important::

  Here again, these forcing files are defined as numpy zip array (**.npz**). These files use specific keys to identify the tectonic forcing. For vertical only condition, the key ``z`` specified the displacements. For the horizontal condition, the ``xyz`` key is used where `xyz` is a 3D array containing the displacement rates along the x, y and z axis (in m/yr). 

.. note::

  There is no requirement to impose continuous tectonics forcing and you might chose to have periods without displacement by making discontinuous events using the ``start`` and ``end`` keys. 

Flexural isostasy definition
-----------------------------------

This function computes the flexural isostasy equilibrium based on topographic change. It is a simple routine that accounts for flexural isostatic rebound associated with erosional loading/unloading.

.. warning::
    This function assumes a value of 1011 Pa for Young's modulus, 0.25 for Poisson's ratio and 9.81 m/s2 for g, the gravitational acceleration.

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: python

            flexure: 
                regdx: 200.
                thick: 20.e3
                rhoc: 2800.0
                rhoa: 3150.0

        Used to consider flexural isostasy in 2D simulation (*i.e.* not global scale).  paleo-topography maps obtained from backward models, you will also have to set this key composed of 2 parameters:

        a. ``regdx``: the resolution of the regular grid used to perform the flexural isostasy calculation,
        b. ``thick`` effective elastic plate thickness in m,
        c. ``rhoc`` crust density in kg/m3,
        d. ``rhoa`` asthenospheric density in kg/m3.
