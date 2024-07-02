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
                  mapH: 'disp20Ma'
                - start: -19000000.
                  end: -18000000.
                  mapH: 'disp19Ma'
                - start: -18000000.
                  end: -17000000.
                  mapH: 'disp18Ma'
                - start: -17000000.
                  end: -16000000.
                  mapH: 'disp17Ma'
                  mapV: 'dispv17Ma'
                - start: -16000000.
                  end: -15000000.
                  mapV: 'dispv16Ma'

        It defines the tectonic forcing conditions from a sequence of events defined by a starting and ending time (``start`` and ``end``) and either a vertical only forcing (*e.g.* uplift and/or subsidence defined with ``mapV``) or a fully 3D displacement mesh ``mapH``. **These displacement rates are set in metres per year**.


.. important::

  Here again, these forcing files are defined as numpy zip array (**.npz**). These files use specific keys to identify the tectonic forcing. For vertical only condition, the key ``z`` specified the displacements. For the horizontal condition, the ``xyz`` key is used where ``xyz`` is a 3D array containing the displacement rates along the x, y and z axis (in m/yr). 

.. note::

  There is no requirement to impose continuous tectonics forcing and you might chose to have periods without displacement by making discontinuous events using the ``start`` and ``end`` keys. 


Plate forcing parameters
----------------------------

Alternatively to the horizontal advection velocity rates proposed in the previous section, one might use the following approach where the plate advection parameters required for interpolation are already set. 

.. note::
    This approach does not need to build a `SciPy cKDTree <https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html>`_ to perform the inverse weighting distance interpolation as the neighborhing fields are already pre-calculated.


.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: python

            plates:
                - start: -20000000.
                  plate: 'plate20Ma'
                  upsub: 'vdisp20Ma'
                - start: -15000000.
                  plate: 'plate15Ma'
                  upsub: 'vdisp15Ma'
                - start: -10000000.
                  plate: 'plate10Ma'
                  upsub: 'vdisp10Ma'

        Plate related horizontal displacements ``plate`` are performed at specified ``start`` time whereas vertical displacements (``upsub``) are done at ``dt`` intervals. Like above, the ``upsub`` are set in metres per year. Both files are numpy zip arrays (**.npz**) and require specific keys. 

        .. important::

            1. For the plate advection information file ``plate``:
                - ``clust``: the cluster of nodes used for interpolation,
                - ``cngbh`` the indices of the nodes in the considered cluster neighborhood,
                - ``dngbh`` the distances between the advected nodes and the mesh,
                - ``ingbh`` the nodes that remain at the same position after advection.
            2. For the vertical displacements mesh ``upsub``:
                - ``t`` the uplift/subsidence tectonic rates,
                - ``z`` the paleo-elevation values 

            There is no requirement to impose both of these files and in the ``upsub`` mesh you can specify either ``z`` or ``t`` or both. If you do define ``z`` then your simulation is forced to fit with the paleo-elevation values.


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

        Used to consider flexural isostasy in 2D simulation (*i.e.* not global scale) where:

        a. ``regdx``: the resolution of the regular grid used to perform the flexural isostasy calculation,
        b. ``thick`` effective elastic plate thickness in m,
        c. ``rhoc`` crust density in kg/m3,
        d. ``rhoa`` asthenospheric density in kg/m3.


.. Global flexural isostasy declaration
.. -------------------------------------

.. This section computes the flexural isostasy equilibrium based on topographic change at global scale. Like previous section, it uses a simple routine that accounts for flexural isostatic rebound associated with erosional loading/unloading.

.. .. note::

..     This function is a simple hack to compute global flexural response based on 2D local loading changes and is likely not the best way for solving this problem. Preferred methods would consist in using spherical harmonics instead...

.. .. warning::

..     This function uses the following default variable values:
..       - Acceleration due to gravity: 9.8
..       - Young's Modulus: 65e9  
..       - Poisson's Ratio: 0.25
..       - Mantle density: 3300.0
..       - Infill material density: 2300.0

.. .. grid:: 1
..     :padding: 3

..     .. grid-item-card::  
        
..         **Declaration example**:

..         .. code:: python

..             gflex: 
..                 interpS: 'data/sflex_info'
..                 interR: 'data/rflex_info'
..                 procs: 8
..                 step: 5.0e5
..                 young: 65e9
..                 poisson: 0.25
..                 rhom: 3300.


..         Used to consider flexural isostasy in 2D simulation (*i.e.* not global scale).  paleo-topography maps obtained from backward models, you will also have to set this key composed of 2 parameters:

..         a. ``interpS``: spherical mesh information used to interpolate goSPL variable to a lon/lat grid,
..         b. ``interR`` regular grid information used to perform interpolation from lon/lat grid to the spherical mesh,
..         c. ``procs`` number of CPUs to use when performing the flexural isostasy calculation,
..         d. ``step`` time interval in years used to compute flexural response.

..         .. important::

..             The interpolation information files are **.npz** files containing the following keys:  ``plate``:
..               - ``wghts``: the weights used for interpolation,
..               - ``ids`` the indices of the nodes in the neighborhood of a given vertex,
..               - ``sumwght`` the sum of the weights used for interpolation,
..               - ``oids`` the nodes that remain at the same position after advection.

..             As for now the interpolation assumes a resolution of **0.25 degrees** for the lont/lat grid.