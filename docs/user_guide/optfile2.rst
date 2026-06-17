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

This function computes the flexural isostasy equilibrium based on topographic change. It accounts for flexural isostatic rebound associated with erosional loading/unloading, solving the thin-elastic-plate biharmonic equation
:math:`D\,\nabla^4 w + \Delta\rho\,g\,w = q` with :math:`D = E\,T_e^3/[12(1-\nu^2)]`.

goSPL provides two solvers, selected by ``method``:

- ``'fem'`` (**flat / 2D models, default**) — a parallel finite-volume biharmonic solve **directly on the unstructured mesh**: no regular grid, no gather-to-root, no external dependency, and spatially-varying elastic thickness handled in a single linear solve. This is the only flat-model solver (it replaced the former gFlex ``'FD'`` and ``'FFT'`` approaches).
- ``'global'`` — a spherical-harmonic solve for global (sphere) models, run serially on the root rank.

.. grid:: 1
    :padding: 3

    .. grid-item-card::

        **Declaration example** (flat / 2D model):

        .. code:: yaml

            flexure:
              method: 'fem'
              thick: 30.e3
              rhoc: 2300.0
              rhoa: 3300.0
              young: 65e9
              nu: 0.25
              bcN: "Mirror"
              bcE: "0Slope0Shear"
              bcS: "Mirror"
              bcW: "0Slope0Shear"

        where:

        a. ``method``: ``'fem'`` for flat/2D models (default) or ``'global'`` for global (spherical) models.
        b. ``thick`` effective elastic plate thickness :math:`T_e` in m (uniform; spatially/temporally variable :math:`T_e` is set with the ``temap`` key below).
        c. ``rhoc`` sediment/crust (load) density in kg/m³,
        d. ``rhoa`` asthenospheric (mantle) density in kg/m³,
        e. ``young``  Young's Modulus in Pa,
        f. ``nu`` Poisson ratio,
        g. ``bcN``, ``bcE``, ``bcS``, ``bcW`` North, East, South and West boundary conditions (``'fem'`` only — see the boundary-condition note below).

        The following keys apply to **both** methods:

        h. ``interval``: apply flexure only every ``interval`` goSPL steps, accumulating the load in between (default ``1`` = every step). The elastic response is instantaneous, so this only coarsens the temporal resolution of the loading history.

        For the ``'global'`` ``method`` (spherical-harmonic solve) the following additional keys are available:

        i. ``ninterp``: number of source points for the mesh ↔ Driscoll–Healy grid interpolation (default ``4``; unused by ``'fem'``).
        j. ``res_deg``: resolution (degrees) of the Driscoll–Healy grid (default ``0.25``). The deflection is long-wavelength, so a coarser grid (e.g. ``0.5`` or ``1.0``) is markedly cheaper with little change to the result — the cheapest speed-up for the serial global solve.
        k. ``maxIter``: maximum Anderson/Picard iterations for the **varying-Te** (``temap``) iterative solve (default ``50``).
        l. ``tol``: relative convergence tolerance ``|dw|/|w|`` for that iteration (default ``5.e-4``).
        m. ``relax``: under-relaxation factor (0–1) for the varying-Te update (default ``1.0``; lower damps strong Te contrasts).

        .. note::

          The varying-Te global solve warm-starts each step from the previous step's converged deflection, so after the first step it typically needs noticeably fewer iterations. With a constant ``thick`` (no ``temap``) the global solve is a single spectral pass and ``maxIter``/``tol``/``relax`` are unused.

        .. note::

          For flat/2D (``'fem'``) models the per-side boundary conditions ``bcN/bcE/bcS/bcW`` (default ``0Slope0Shear``) can be:

            - ``0Slope0Shear``: zero slope and zero shear at the edge (``∂w/∂n = 0``, ``∂³w/∂n³ = 0``) — the natural boundary; use it (or ``Mirror``) when the load is far from the edge.
            - ``Mirror``: reflective/symmetry boundary. For a thin plate this is identical to ``0Slope0Shear`` and uses the same (natural) treatment.
            - ``0Displacement0Slope``: clamped edge (``w = 0`` and ``∂w/∂n = 0``).

          ``0Moment0Shear`` and ``Periodic`` are not available with the ``'fem'`` solver.

        .. note::

          The ``'fem'`` solver caches its operator and factorisation and reuses them every step (serial direct LU / parallel MUMPS), so it is fast even on small meshes; only the surface load changes between steps. The boundary discretisation differs from the former gFlex solver by ~10 % only where the flexural wavelength approaches the domain size (the deflection reaches the edges); they agree closely otherwise.


In addition, it is possible to define a spatially-variable (and time-varying) lithospheric elastic thickness with the ``temap`` key below, where specific maps can be added through time during the run. This works with both ``'fem'`` and ``'global'`` methods (for ``'fem'`` the varying :math:`T_e` is still a single linear solve).

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