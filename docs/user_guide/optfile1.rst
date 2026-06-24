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

        .. code:: yaml

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

        .. code:: yaml

            climate:
                - start: 1000.
                  map: ['rain_map','r']
                - start: 20000.
                  uniform: 1.
                - start: 30000.
                  zscale: [0.0008, 0.5]


        The climatic forcing is defined in a similar fashion as the tectonic one. It is based on a sequence of events with each event starting at a given time (``start`` in years) and corresponding to a given precipitation condition. Each event must specify **one** of three rainfall sources:

        a. ``uniform`` — a single rainfall value (m/yr) applied across the whole mesh.
        b. ``map`` — a path-and-key pair pointing at a precipitation grid stored in an ``.npz`` file (e.g. ``['rain_map', 'r']`` where ``'r'`` is the key holding the per-vertex values, in m/yr).
        c. ``zscale`` — a two-element list ``[A, B]`` that makes the rainfall scale **linearly with elevation** at every node and every time step: :math:`\mathrm{P_i = A \cdot \eta_i + B}` (clamped to non-negative values), where :math:`\mathrm{\eta_i}` is the current surface elevation in metres and :math:`\mathrm{P_i}` is the resulting precipitation rate in m/yr. ``A`` is therefore in ``(m/yr) per m`` and ``B`` is the sea-level intercept in ``m/yr``.

        The precipitation field is re-evaluated at every time step (defined by ``dt``). For ``zscale`` this means the rain pattern follows the evolving topography — newly uplifted regions wet, eroded regions dry — without needing pre-baked maps for every snapshot.

.. important::

    When defining a precipitation grid, one needs to use the **npz** format and needs to specify the key corresponding to the precipitation variable in the file. In the above example this key is ``'r'``. The precipitation grid needs to define values for all vertices in the mesh.

.. note::

    The ``zscale`` option is a cheap proxy for orographic precipitation: it captures the first-order increase in rainfall with elevation but ignores wind direction, atmospheric dynamics, and rain-shadow effects. For the full orographic precipitation model — including airflow, condensation, and downwind drying — see the dedicated *Orographic precipitation* section below.

    Negative values of :math:`\mathrm{P_i}` (which can occur for low ``B`` and sub-sea-level nodes) are clamped to zero before being routed through the flow accumulation, so coastal depressions do not act as moisture sinks.


Evaporation (optional)
^^^^^^^^^^^^^^^^^^^^^^

Each climate event may additionally declare an **evaporation rate** that is subtracted from runoff before river flow accumulation and from lake-surface inflow before pit overflow. Evaporation is opt-in per event — any row that omits both keys contributes zero evaporation for the corresponding interval, so existing input files are unaffected.

.. grid:: 1
    :padding: 3

    .. grid-item-card::

        **Declaration example**:

        .. code:: yaml

            climate:
                - start: 0.
                  uniform: 1.0
                  evap_uniform: 0.3
                - start: 50000.
                  map: ['rain_map', 'r']
                  evap_map: ['evap_map', 'e']
                - start: 100000.
                  uniform: 0.5
                  # no evap_* key — evaporation is zero for this interval

        Two source types are supported:

        a. ``evap_uniform`` — a single evaporation rate (m/yr) applied at every land node.
        b. ``evap_map`` — a path-and-key pair pointing at a per-vertex evaporation grid stored in an ``.npz`` file, matching the same convention as the ``map`` rainfall source.

        Elevation-banded evaporation (an analogue of ``zscale``) is not yet supported.

.. note::

    Evaporation is applied as a **single rate**, used identically over channels and over lake surfaces. Cells with ``evap > rain`` produce zero runoff (the excess evaporative capacity is dropped — groundwater is out of scope). Lakes whose lake-surface evaporation budget exceeds their inflow simply do not form; ``waterFill`` stays at the bare-topography level in those depressions.

    The lake-evaporation budget per depression is computed using the **maximum-fill surface area** (the full extent the lake would occupy at spillover) rather than the actual fill area. This is intentionally conservative: it over-estimates lake evaporation at partial-fill levels, which biases the model towards drier basins — the physically appropriate direction when potential evaporation is the limiting term.

    Channels and lakes share the same per-cell rate; users who need a higher lake-surface rate (e.g. to capture open-water evaporation over endorheic basins) should encode the spatial pattern directly in ``evap_map``.

.. important::

    The ``Evap`` per-node field is written to the HDF5/XDMF output whenever the YAML declares any evaporation; it can be visualised in Paraview alongside ``Rain`` to verify the input pattern. The reduced flow accumulation ``FA``, the reduced lake levels in ``waterFill``, and the lake-erosion proxy ``fillFA`` all propagate the evaporation losses automatically.


Orographic precipitation definition
------------------------------------

.. warning::

    The orographic precipitation is only available for 2D grids and will not be working for global models.

goSPL implements the Linear Theory of Orographic Precipitation following `Smith & Barstad (2004) <https://journals.ametsoc.org/view/journals/atsc/61/12/1520-0469_2004_061_1377_altoop_2.0.co_2.xml>`_.

.. note::

    The model consists of two vertically-integrated steady-state advection equations describing:

    1. the cloud water density and
    2. the hydrometeor density,

    forced by the terrain-driven uplift ``Cw·(v·∇h)`` (condensation on windward slopes, evaporation on lee slopes). goSPL solves these **directly on the unstructured mesh in parallel** (a first-order upwind finite-volume discretisation), so no regular grid or FFT is required. The windward/lee rain shadow is reproduced. The stratified mountain-wave (airflow) term of the original spectral solution is dropped, so the ``latitude``, ``nm`` and ``hw`` parameters are inert; this term mainly refines the precipitation pattern for narrow or strongly stratified ranges.

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: yaml

            orography:
                wind_speed: 10.0
                wind_dir: 270
                env_lapse_rate: -4
                moist_lapse_rate: -7
                ref_density: 7.4e-3
                conv_time: 1000.
                fall_time: 1000.
                precip_base: 7.0
                precip_min: 0.01
                rainfall_frequency: 1

        This part of the input file defines the parameters for the orographic rain. All keys are optional and fall back to the defaults given below.

        **Wind** (sets the advection of moisture over the terrain):

        a. ``wind_speed``: uniform wind speed [m/s]; default 10. Together with the timescales below it sets the *advection length* ``wind_speed × conv/fall_time`` — how far moisture drifts downwind before it rains out. Faster wind pushes the precipitation maximum further onto the lee.
        b. ``wind_dir``: wind direction the flow comes **from** (meteorological convention) [degrees, ``0`` = north, ``90`` = east, ``180`` = south, ``270`` = west]; default 0. The wind therefore *blows toward* the opposite bearing — e.g. ``wind_dir: 270`` is a westerly that blows **eastward**, and ``wind_dir: 0`` is a northerly that blows **southward**. It fixes which slopes are windward (uplift, rain) and which are lee (subsidence, rain shadow).

        **Moisture sensitivity** (combine into the uplift coefficient ``Cw = ref_density × moist_lapse_rate / env_lapse_rate``, which scales the condensation source ``Cw·(v·∇h)`` and hence the overall orographic rain intensity):

        c. ``env_lapse_rate``: environmental lapse rate [°C/km]; default -4.0 (must be non-zero).
        d. ``moist_lapse_rate``: moist adiabatic lapse rate [°C/km]; default -7.0.
        e. ``ref_density``: reference saturation water-vapour density [kg/m³]; default 7.4e-3.

        **Microphysical timescales** (control where the rain falls relative to the crest and how sharp the rain shadow is):

        f. ``conv_time``: cloud-water → hydrometeor conversion time :math:`\\tau_c` [s]; default 1000.
        g. ``fall_time``: hydrometeor fallout time :math:`\\tau_f` [s]; default 1000. Short timescales make the rain fall on/near the windward crest with a crisp dry lee; long timescales smear it downwind.

        **Rate conversion**:

        h. ``precip_base``: uniform background (non-orographic) precipitation rate added everywhere [mm/h]; default 7.0.
        i. ``precip_min``: minimum precipitation floor [mm/h] (the dried lee shadow is clamped to this); default 0.01.
        j. ``rainfall_frequency``: number of 1-hour storms per day, used to convert the rate to m/yr; default 1.

.. note::

    The orographic precipitation is solved **directly on the unstructured mesh and in parallel** — there is no regular grid, no FFT and no ``regdx`` parameter. The stratified mountain-wave (airflow) term of the original spectral Smith-Barstad solution is dropped; the windward/lee rain shadow is retained, so the parameters ``latitude``, ``nm`` and ``hw`` of earlier versions are no longer used (and are silently ignored if present).

    The orographic uplift is forced by the surface the airflow follows, so the elevation is **clamped to sea level**: submarine bathymetry produces no orographic rain (the air flows over the flat sea surface), and the forcing only switches on where land rises above the sea.