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

        .. code:: yaml

            tectonics:
                - start: -20000000.
                  end: -19000000.
                  upsub: ['data/uplift20Ma','t']
                - start: -18000000.
                  end: -17000000.
                  upsub: ['data/uplift18Ma','t']
                  hdisp: ['data/hdisp18Ma', 'hxyz']

        It defines the tectonic forcing conditions from a sequence of events defined by a starting and ending time (``start`` and ``end``) and either a vertical rate only forcing (*e.g.* uplift and/or subsidence defined with ``upsub``) or a fully 3D displacement rate ``hdisp``. **These displacement rates are set in metres per year**. 

        .. important::

          For horizontal advection the user needs to define the ``advect`` key in the ``domain`` `section <https://gospl.readthedocs.io/en/latest/user_guide/inputfile.html#initial-mesh-definition-and-simulation-declaration>`_ of the input file. The advection scheme could either be ``upwind``, ``iioe1``, ``iioe2`` or ``interp``  (go to the technical `information <https://gospl.readthedocs.io/en/latest/tech_guide/tecto.html#horizontal-advection>`_ in the documentation for more information). 


.. important::

  Here again, these forcing files are defined as numpy zip array (**.npz**). These files use specific keys to identify the tectonic forcing that are specified by the user in the input file. For vertical only condition, the displacements (in m/yr) is a **1D vector** with values on each node of the grid. For the horizontal condition, the key is a **3D array** containing the displacement rates along the x, y and z axis (in m/yr). When the horizontal advection is for 2D grids, the provided displacements also need to be in 3D but with the third dimension set to 0.0.

.. note::

  There is no requirement to impose continuous tectonics forcing and you might chose to have periods without displacement by making discontinuous events using the ``start`` and ``end`` keys. 

.. note::

  When applying horizontal displacement using the ``interp`` scheme (mainly used in global simulation to represent plate movements), the horizontal movements are performed at the end of the period (i.e., at the specified ``end`` time). In the other cases, the horizontal displacement is done at every timestep (specified by ``dt``).


.. Plate forcing parameters
.. ----------------------------

.. Alternatively to the horizontal advection velocity rates proposed in the previous section, one might use the following approach where the plate advection parameters required for interpolation are already set. 

.. .. note::
..     This approach does not need to build a `SciPy cKDTree <https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html>`_ to perform the inverse weighting distance interpolation as the neighborhing fields are already pre-calculated.


.. .. grid:: 1
..     :padding: 3

..     .. grid-item-card::  
        
..         **Declaration example**:

..         .. code:: python

..             plates:
..                 - start: -20000000.
..                   plate: 'plate20Ma'
..                   upsub: 'vdisp20Ma'
..                 - start: -15000000.
..                   plate: 'plate15Ma'
..                   upsub: 'vdisp15Ma'
..                 - start: -10000000.
..                   plate: 'plate10Ma'
..                   upsub: 'vdisp10Ma'

..         Plate related horizontal displacements ``plate`` are performed at specified ``start`` time whereas vertical displacements (``upsub``) are done at ``dt`` intervals. Like above, the ``upsub`` are set in metres per year. Both files are numpy zip arrays (**.npz**) and require specific keys. 

..         .. important::

..             1. For the plate advection information file ``plate``:
..                 - ``clust``: the cluster of nodes used for interpolation,
..                 - ``cngbh`` the indices of the nodes in the considered cluster neighborhood,
..                 - ``dngbh`` the distances between the advected nodes and the mesh,
..                 - ``ingbh`` the nodes that remain at the same position after advection.
..             2. For the vertical displacements mesh ``upsub``:
..                 - ``t`` the uplift/subsidence tectonic rates,
..                 - ``z`` the paleo-elevation values 

..             There is no requirement to impose both of these files and in the ``upsub`` mesh you can specify either ``z`` or ``t`` or both. If you do define ``z`` then your simulation is forced to fit with the paleo-elevation values.


Flexural isostasy definition
-----------------------------------

This function computes the flexural isostasy equilibrium based on topographic change. It accounts for flexural isostatic rebound associated with erosional loading/unloading.

.. important::

  It is computed on a regular grid (in X,Y in 2D or lon/lat for global simulation) and is limited in terms of parallelisation.

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: yaml

            flexure: 
              method: 'FD'
              regdx: 1000.
              ninterp: 4
              thick: 30.e3
              rhoc: 2300.0
              rhoa: 3300.0
              young: 65e9
              nu: 0.25
              bcN: "Mirror"
              bcE: "0Slope0Shear"
              bcS: "Mirror"
              bcW: "0Slope0Shear"

        Used to consider flexural isostasy (*i.e.* not global scale) where:

        a. ``method``: the approach used is either 'FD' in 2D or 'global' for global model.  
        b. ``regdx``: the resolution of the regular grid used to perform the flexural isostasy calculation,
        c. ``ninterp``: the number of points used to perform the interpolation between goSPL unstructured mesh and the regular grid (not used for the 'global' ``method``) 
        d. ``thick`` effective elastic plate thickness in m,
        e. ``rhoc`` crust density in kg/m3,
        f. ``rhoa`` asthenospheric density in kg/m3.
        g. ``young``  Young's Modulus in Pa.
        h. ``nu`` Poisson ratio.
        i. ``bcN``, ``bcE``, ``bcS``, ``bcW`` North, East, South and West boundary conditions.

        .. note::
  
          For non-global simulation, the user needs to specify the boundary conditions for the flexural isostasy calculation. Similar conditions to `gFlex <https://gmd.copernicus.org/articles/9/997/2016/gmd-9-997-2016.pdf>`_ are possible:
            
            - `0Displacement0Slope` 0-displacement-0-slope boundary condition

            - `0Moment0Shear`: Broken plate boundary condition second and third derivatives of vertical displacement are 0. This is like the end of a diving board.
            - `0Slope0Shear`: First and third derivatives of vertical displacement are zero. While this does not lend itsellf so easily to physical meaning, it is helpful to aid in efforts to make boundary condition effects disappear (i.e. to emulate the NoOutsideLoads cases)
            - `Mirror`: Load and elastic thickness structures reflected at boundary.
            - `Periodic``: Wrap-around boundary condition: must be applied to both North and South and/or both East and West. This causes, for example, the edge of the eastern and western limits of the domain to act like they are next to each other in an infinite loop.

.. warning::

    In case where **flexure** and **orographic rain** capabilities are defined in the same simulation, you will need to have the same grid resolution (``regdx``) for each definition.


In addition, it is possible to define variable lithospheric elastic thicknesses by using the ``temap`` key below where specific maps could be added through time during the run.

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: yaml

            temap:
              - start: 250.e6
                map: ['input/temap1', 'te']
              - start: 251.e6
                map: ['input/temap2', 'te']

Here again, the elastic maps files are provided as numpy zip array (**.npz**).

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