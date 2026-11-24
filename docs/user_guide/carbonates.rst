.. _carbonates:

==============================
Carbonate platforms and reefs
==============================

goSPL can produce carbonate sediment in the sea and record it in the
stratigraphy as its own lithologies, next to the river (clastic) sediment.
Three opt-in blocks work together:

- ``carbonate:`` **benthic production** on the seafloor by one or more
  carbonate factories (reef builders, shallow and deeper producers). Part of
  each factory's production is built in place as reef framework; the rest is
  loose carbonate sand and mud spread by the marine transport.
- ``pelagic:`` **pelagic carbonate rain** (plankton) preserved above the
  carbonate compensation depth, and **hemipelagic clay** settling from the
  water column.
- ``strata: lithologies:`` (expert) the **lithology list** that both blocks
  use: each lithology has its own pile in the stratigraphy, compaction curve,
  erodibility and transport behaviour.

The same equations build carbonate platforms on coarse meshes and reef-rimmed
platforms, barriers and lagoons on fine ones: what changes with the resolution
is what the mesh resolves, and the two resolution-dependent choices you make
(the open-water length and the framework fractions). See
:ref:`carb_resolution`. The equations are on
:ref:`the carbonate page of the technical guide <carb>`.

.. contents:: On this page
   :local:
   :depth: 1


Requirements
------------

- ``strata: dual: True`` and a ``strat`` interval in the ``time`` block: the
  carbonate is recorded as lithologies of the dual-lithology stratigraphy.
- ``domain: seadepo: True`` (and ``nodep: False``) for the ``carbonate``
  block: the loose carbonate is spread by the marine transport.
- The ``pelagic`` block needs only the stratigraphy (it settles in place,
  unless ``resuspend`` is set).

goSPL stops with a message when one of these is missing. Without the blocks,
every run is unchanged (dual-lithology runs are bitwise identical to before).


Quick start
-----------

.. code:: yaml

    time:
        start: 0.
        end: 1.e6
        tout: 1.e4
        dt: 1.e3
        strat: 1.e4

    strata:
        dual: True

    carbonate:
        insolation: 400.      # three default factories (euphotic, oligophotic, aphotic)

This produces carbonate on every marine cell with the three factories of
CarboKitten (Hidding et al., 2026), builds part of it in place as ``reef``,
and spreads the rest as ``carb_sand`` and ``carb_mud``. Add a ``pelagic``
block for the open ocean (:ref:`carb_pelagic`).


How it works
------------

**Factories.** Each factory produces carbonate at the Bosscher and Schlager
(1992) rate, which is saturated in shallow water and falls off with depth as
light is absorbed:

.. code:: text

    rate = growth * tanh(insolation * exp(-extinction * depth) / saturation)

The three default factories are those of CarboKitten:

.. list-table::
   :header-rows: 1
   :widths: 18 14 14 14 14 12

   * - Factory
     - growth (m/Myr)
     - extinction (1/m)
     - saturation (W/m²)
     - framework
     - sand
   * - euphotic
     - 500
     - 0.8
     - 60
     - 0.5
     - 0.7
   * - oligophotic
     - 400
     - 0.1
     - 60
     - 0.1
     - 0.5
   * - aphotic
     - 100
     - 0.005
     - 60
     - 0
     - 0.3

The growth, extinction and saturation are CarboKitten's; the framework and
sand fractions are illustrative defaults to calibrate.

**Growth within a step.** The production of a step is integrated with the
depth it creates, so a cell never builds above sea level (minus
``min_depth``) and the result does not depend on the time step. Because
production increases as the water shallows, there is a critical depth
:math:`w^*` where production equals the rate of relative sea-level rise: a
platform shallower than :math:`w^*` catches up with the sea and keeps up with
it; a deeper one falls behind and drowns (give-up). goSPL reproduces this
exactly (benchmark in the technical guide).

**Framework and loose sediment.** A fraction ``framework`` of each factory's
production is built in place and recorded as the lithology ``reef`` (low
erodibility and porosity by default). The rest is loose: a fraction ``sand``
of it becomes ``carb_sand`` (it behaves like coarse sediment) and the rest
``carb_mud`` (like fine sediment). The loose carbonate joins the marine
transport of the same step, with the river sediment: the clinoform cascade,
the marine diffusion, and later slope failure or erosion by gravity currents.
River gravity currents do not carry it.

**What controls production.** Besides the depth (through light), each factory
can be scaled by:

- a ``habitat`` map (0 to 1, one value per mesh vertex), for instance from
  the sea-surface temperature or the latitude on a global run;
- ``response`` tables on four controls (below), multiplied together.

**Open-water index.** Reefs grow best at the platform margin, facing open
water. The ``exposure`` control is the elevation minus its average over the
length ``exposure_length``: positive on margins, rims and highs, negative in
lagoons and behind barriers. A table that cuts production at negative
``exposure`` concentrates a factory on the margin.


Parameters
----------

**The** ``carbonate`` **block**

.. list-table::
   :header-rows: 1
   :widths: 20 12 14 54

   * - Key
     - Unit
     - Default
     - Meaning
   * - ``insolation``
     - W/m²
     - 400
     - Light at the sea surface :math:`I_0`, shared by all factories.
   * - ``min_depth``
     - m
     - 0
     - Production stops this far below sea level (the top of a reef flat).
   * - ``exposure_length``
     - m
     - none
     - Length of the open-water index; required when a factory uses the
       ``exposure`` control. A few cells long: the width of a lagoon or of a
       platform top.
   * - ``factories``
     - list
     - the three above
     - Carbonate factories (next table). Giving the list replaces the
       defaults.
   * - ``lithologies``
     - mapping
     - see below
     - Changes to the properties of ``carb_sand``, ``carb_mud`` and ``reef``.

**A factory**

.. list-table::
   :header-rows: 1
   :widths: 20 12 14 54

   * - Key
     - Unit
     - Default
     - Meaning
   * - ``name``
     - text
     - factory<k>
     - Label (unique).
   * - ``growth``
     - m/Myr
     - 0
     - Maximum production rate :math:`g`.
   * - ``extinction``
     - 1/m
     - 0.1
     - Light extinction coefficient :math:`k`: high for light-dependent
       (euphotic) producers.
   * - ``saturation``
     - W/m²
     - 60
     - Light at which production saturates :math:`I_k`.
   * - ``framework``
     - 0 to 1
     - 0
     - Fraction built in place (``reef``). Corals and other frame builders:
       high; pelagic and most deeper producers: 0.
   * - ``sand``
     - 0 to 1
     - 0.5
     - Fraction of the loose part that is sand (``carb_sand``); the rest is
       mud (``carb_mud``).
   * - ``habitat``
     - [file, key]
     - none
     - Map (npz archive, one value per mesh vertex, clipped at 0) multiplying
       the production.
   * - ``response``
     - list
     - none
     - Tables ``{control: ..., x: [...], f: [...]}``: the production is
       multiplied by ``f`` interpolated at the control value (``x`` increasing,
       ``f`` >= 0, constant beyond the ends).

**Response controls**

.. list-table::
   :header-rows: 1
   :widths: 16 12 72

   * - ``control``
     - Unit
     - Value used
   * - ``depth``
     - m
     - Water depth, inside the integration (an extra depth limit on top of the
       light law, e.g. a narrow depth window for a reef factory).
   * - ``slope``
     - m/m
     - Seafloor gradient (steep fore-reef or slope).
   * - ``sedrate``
     - m/yr
     - Clastic (non-carbonate) marine deposition rate of the previous step:
       turbidity and burial near river mouths. ``x: [0., 1.e-3], f: [1., 0.]``
       stops production where 1 mm/yr of river sediment accumulates.
   * - ``exposure``
     - m
     - Open-water index (needs ``exposure_length``): positive on margins.
       ``x: [-5., 5.], f: [0.2, 1.]`` favours the margin over the lagoon.
   * - ``wave``
     - 0 to 1
     - Relative wave energy from the ``waves`` block (below): 1 in the open
       ocean, 0 on land and in the lee of islands and headlands.
       ``x: [0., 0.5], f: [0.1, 1.]`` keeps reef builders on exposed margins.
   * - ``temperature``
     - °C
     - Sea-surface temperature from the ``sst`` block (below). For reef
       builders ``x: [16., 18., 28., 32.], f: [0., 1., 1., 0.]``: no growth
       below 16 °C or above 32 °C, full growth between 18 and 28 °C.

**Carbonate lithologies** (``carbonate: lithologies:``)

.. list-table::
   :header-rows: 1
   :widths: 14 12 12 12 14 14 22

   * - Lithology
     - size
     - phi0
     - z0 (m)
     - k_factor
     - diff_factor
     - Origin
   * - ``carb_sand``
     - coarse
     - 0.45
     - 2500
     - 1
     - 1
     - loose sand
   * - ``carb_mud``
     - fine
     - 0.65
     - 1500
     - 1
     - 1
     - loose mud, pelagic rain
   * - ``reef``
     - coarse
     - 0.2
     - 5000
     - 0.1
     - 0.1
     - framework (only if a factory has ``framework`` > 0)

``phi0`` and ``z0`` set the porosity-depth curve
:math:`\phi = \phi_0 e^{-z/z_0}` used for compaction; ``k_factor`` multiplies
the river erodibility of an exposed layer and ``diff_factor`` the hillslope
diffusivity, in proportion to the lithology's share of the surface layer;
``size`` decides how pits, the marine spreading and the gravity currents treat
it (coarse settles first, fine travels further). The values are orders of
magnitude from published carbonate compaction trends (grainstone, lime mud,
boundstone), to adapt to your system::

    carbonate:
        lithologies:
            reef: {phi0: 0.1, k_factor: 0.05}


.. _carb_waves:

Wave energy and shadows
-----------------------

.. code:: yaml

    waves:
        climate:                       # directions the waves come FROM
          - {from: 90., weight: 0.7}   # degrees: 0 = north, 90 = east
          - {from: 180., weight: 0.3}  # or from: zonal, or from_map: [file, key]
        fetch: 3.e5                    # growth to full energy over open water (m)
        decay: 1.e5                    # attenuation length in deep water (m)
        shoal_depth: 20.               # waves feel the bottom shallower than this (m)
        shoal_length: 2.e3             # attenuation length at the shoreline (m)

The ``waves`` block computes a relative wave energy (0 to 1, output
``waveE``): for each direction of the climate the energy enters through the
domain edges facing the waves, travels across the sea, is lost in shallow
water (more and more as the depth drops below ``shoal_depth``) and stops on
land, so islands, headlands and barriers cast a shadow behind them. The
directions are combined with their weights (normalised to 1). It is the
``wave`` response control of the factories, and the energy of the coastal
erosion planned in phase C4.

*Directions.* ``from`` is one bearing everywhere; ``from: zonal`` gives the
global wind belts (trade winds from the east within 30° of the equator,
westerlies to 60°, polar easterlies beyond; global meshes only);
``from_map: [file, key]`` gives a bearing per mesh vertex (degrees), for any
other pattern.

*Fetch.* With ``fetch`` the waves grow to full energy over that length of open
water, and grow again behind islands (fetch-limited waves). A global mesh has
no domain edges to bring the waves in, so ``fetch`` is required there (a few
hundred km); on a flat mesh the edges facing the waves bring the open-ocean
energy, and ``fetch`` is optional.

Choose ``decay`` much longer than the domain if the open ocean should keep its
energy, and ``shoal_length`` of the width of the surf zone to the reef flat.
The waves follow the wind directions: refraction (turning towards the shore in
shallow water) is planned as an option.

.. _carb_sst:

Sea-surface temperature
-----------------------

.. code:: yaml

    sst:
        events:                         # global: maps through time
          - {start: 0., map: ['sst_0Ma', 't']}
          - {start: 5.e6, map: ['sst_5Ma', 't']}
        curve: 'sst_anomaly.csv'        # optional: time (yr), SST (°C)

The ``sst`` block gives the sea-surface temperature used by the
``temperature`` response of the factories and by the ``temperature`` table of
the pelagic rain. Events hold until the next one starts; a map has one value
per mesh vertex (°C). The curve (two columns, time in years and SST in °C,
comma-separated for a ``.csv`` file) is interpolated linearly in time and
constant beyond its ends. With events and a curve, the curve is an anomaly
added to the event values (for instance glacial cooling on top of a present-day
map); with a curve alone it is the SST everywhere, which suits a regional
model::

    sst:
        curve: 'sst_regional.csv'

Each event can also be uniform (``{start: 0., uniform: 26.}``). The SST is
written as the output ``sst``.

.. _carb_pelagic:

Pelagic carbonate and hemipelagic clay
--------------------------------------

.. code:: yaml

    pelagic:
      - start: 0.
        carbonate: 5.         # pelagic carbonate rain (m/Myr) or [file, key]
        lysocline: 3500.      # depth (m) where dissolution starts
        ccd: 4500.            # carbonate compensation depth (m)
        clay: 1.              # hemipelagic clay (m/Myr) or [file, key]
        clay_decay: 5.e4      # optional: clay * exp(-d/L) from the coast (m)
        resuspend: False      # optional: hand the rain to the marine diffusion
        temperature: {x: [5., 25.], f: [0.2, 1.]}   # optional, needs `sst`
      - start: 2.e6           # a later event replaces the earlier one
        carbonate: ['prod2Ma', 'p']
        lysocline: 3800.
        ccd: 4800.
        clay: 1.

.. list-table::
   :header-rows: 1
   :widths: 18 14 14 54

   * - Key
     - Unit
     - Default
     - Meaning
   * - ``start``
     - yr
     - required
     - The event holds from this time until the next event starts.
   * - ``carbonate``
     - m/Myr or [file, key]
     - 0
     - Pelagic carbonate rain (uniform, or a map such as an export
       production field converted to a thickness rate).
   * - ``lysocline``
     - m
     - required with carbonate
     - Depth below which the rain starts to dissolve.
   * - ``ccd``
     - m
     - = lysocline
     - Carbonate compensation depth: no carbonate preserved below it; linear
       between the two.
   * - ``clay``
     - m/Myr or [file, key]
     - 0
     - Hemipelagic clay rain (recorded as ``fine``).
   * - ``clay_decay``
     - m
     - none
     - The clay rate becomes the coastal rate times :math:`e^{-d/L}`, ``d``
       the distance from land measured through the sea (around headlands and
       islands). No land at all gives no clay.
   * - ``resuspend``
     - true/false
     - false
     - Hand the settled rain to the marine diffusion of the same step instead
       of laying it in place: it gathers into the lows of the seafloor.
   * - ``temperature``
     - table
     - none
     - ``{x: [...], f: [...]}``: the carbonate rain is multiplied by ``f``
       interpolated at the sea-surface temperature (``sst`` block).

The rain settles on marine cells only and never builds above sea level (the
excess is reported as bypassed). Thin deposits (below 0.1 mm per step) are
carried over until they are thick enough to be recorded. The ``pelagic`` block
can be used with or without the ``carbonate`` block.


.. _carb_resolution:

Setting it up for your resolution
---------------------------------

.. list-table::
   :header-rows: 1
   :widths: 14 30 56

   * - Mesh
     - What forms
     - How to set it up
   * - Global, 10-50 km cells
     - Carbonate shelves and platforms, their growth and drowning with sea
       level; pelagic carbonate and clay in the ocean basins. A reef is
       smaller than a cell: ``reef`` is the framework share of a cell's
       carbonate (a margin facies), not a resolved reef.
     - an ``sst`` block with SST maps through time and a ``temperature``
       response (no reef outside the warm belt), or ``habitat`` maps; a
       ``pelagic`` block with a productivity map and
       a lysocline and CCD per event; ``clay_decay`` of tens to hundreds of
       km; usually no ``exposure`` (cells are wider than lagoons); a ``waves`` block
       with ``from: zonal`` and a ``fetch`` for exposed and sheltered coasts; time steps
       of 10-100 kyr are fine (the production does not depend on ``dt``).
   * - Regional, 0.5-5 km cells
     - Platforms with a margin, an interior and a slope, carbonate aprons,
       mixed carbonate-clastic shelves, drowning and backstepping. Reefs are
       a margin facies a few cells wide.
     - ``exposure`` with ``exposure_length`` of 3-10 cells, or a ``waves``
       climate with a ``wave`` response, to concentrate the euphotic factory on
       exposed margins; ``sedrate`` to stop production near
       deltas; an ``sst`` curve with a ``temperature`` response for climatic
       changes through time; framework 0.3-0.6 for the euphotic factory; ``dt`` of 1-10 kyr.
   * - Local, 50-200 m cells
     - Rimmed platforms, barriers and lagoons, patch reefs several cells
       wide; reef terraces with coastal erosion (planned, phase C4).
     - ``exposure_length`` of a lagoon width (a few hundred metres to a few
       km); a ``depth`` window for the reef factory (for instance full
       production to 10 m, none below 30 m); framework up to 0.8 for reef
       builders; ``min_depth`` for the reef flat; ``dt`` of 0.1-1 kyr.

Structures narrower than about five cells are not resolved at any
resolution, and growth is vertical (as for all goSPL deposition): lateral
progradation comes from the loose sediment. There is no cellular automaton,
so the ecological patchiness of carbonate factories is not represented.

**Global example**

.. code:: yaml

    strata:
        dual: True
    carbonate:
        insolation: 400.
        factories:
          - {name: shallow, growth: 400., extinction: 0.1, saturation: 60.,
             framework: 0.3, sand: 0.6, habitat: ['sst_habitat', 'h']}
          - {name: deep, growth: 100., extinction: 0.005, saturation: 60.,
             framework: 0., sand: 0.3}
    pelagic:
      - start: 0.
        carbonate: ['export_prod', 'p']
        lysocline: 3500.
        ccd: 4500.
        clay: 1.
        clay_decay: 3.e5

**Regional example**

.. code:: yaml

    carbonate:
        exposure_length: 1.e4         # 10 km on a 1 km mesh
        factories:
          - name: euphotic
            growth: 500.
            extinction: 0.8
            saturation: 60.
            framework: 0.5
            sand: 0.7
            response:
              - {control: exposure, x: [-5., 5.], f: [0.2, 1.]}
              - {control: sedrate, x: [0., 1.e-3], f: [1., 0.]}
          - {name: oligophotic, growth: 400., extinction: 0.1, framework: 0.1}
          - {name: aphotic, growth: 100., extinction: 0.005}

**Local example**

.. code:: yaml

    carbonate:
        min_depth: 0.5
        exposure_length: 1500.
        factories:
          - name: reef
            growth: 800.
            extinction: 0.8
            saturation: 60.
            framework: 0.8
            sand: 0.8
            response:
              - {control: depth, x: [0., 10., 30.], f: [1., 1., 0.]}
              - {control: exposure, x: [-2., 2.], f: [0.1, 1.]}
          - {name: lagoon, growth: 300., extinction: 0.2, framework: 0., sand: 0.4}


Outputs
-------

- ``carbProd`` (m/yr): carbonate production rate of the step (all factories).
- ``pelDep`` (m/yr): settled (or resuspended) pelagic thickness rate.
- ``sst`` (°C): sea-surface temperature of the step (with an ``sst`` block).
- ``waveE`` (0 to 1): relative wave energy (with a ``waves`` block, when a
  factory uses it).
- Stratigraphy: one pile per lithology, ``stratH_<name>`` and ``phi_<name>``
  (``stratHf``/``phiF`` for the dual-lithology fine fraction), in the
  ``stratal`` files.
- ``--summary``: a ``carbonate`` budget per step (produced, framework, loose,
  m³) and a ``pelagic`` budget (carbonate rain, settled, dissolved, clay,
  bypassed).

The section and well tools (``gospl.analyse.stratasection``) currently draw
the coarse/fine composition; to plot the carbonate lithologies, read the
``stratH_<name>`` datasets.


With the other options
----------------------

- *Gravity currents.* River gravity currents carry the river load only; the
  carbonate on the slope can still fail (``failure_slope``) or be eroded by the
  currents, and is then carried with them.
- *Provenance.* The provenance classes describe the clastic sources. Loose
  carbonate takes the river composition of the cell where it is produced, and
  the bedrock class offshore.
- *Lakes.* No carbonate is produced in lakes.
- *Flexure, compaction, plate motion.* The carbonate piles load the
  lithosphere, compact with their own curves and move with the stratigraphy
  like any lithology.


Cautions
--------

- *Calibration.* The modules have been checked against exact solutions and
  for conservation and parallel consistency, but not yet calibrated on real
  platforms: start from the defaults, compare with your observations
  (platform thickness and width, margin position, drowning events) and adjust
  the growth rates and framework fractions first.
- *Response tables are multiplied*: several restrictive tables can switch a
  factory off almost everywhere. Check ``carbProd``.
- *The open-water index depends on the length, not the cells*: a length
  shorter than two or three cells gives a noisy index.
- *Budget.* With ``-v`` each step prints the produced, framework and loose
  volumes and the pelagic budget; the loose carbonate leaving through open
  boundaries is the only part not recorded in the stratigraphy.
