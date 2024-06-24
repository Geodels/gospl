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

        Follows the tectonic forcing conditions with a sequence of events defined by a starting time (``start``) and either a vertical only forcing (*e.g.* uplift and/or subsidence defined with ``mapV``) or a fully 3D displacement mesh ``mapH``. These displacements are set in metres per year.


.. important::

  As mentioned above and for the next key parameter as well, these forcing files are defined as numpy zip array (**.npz**).

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
