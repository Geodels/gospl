.. _gw_user:

========================================
Groundwater, duricrust and geochemistry
========================================

Adding an optional ``groundwater`` section turns on a **water table** and a
generic **duricrust**, near-surface hydrology plus chemical armoring of the
erodibility. It is fully **opt-in**: with no ``groundwater`` block the model is
byte-identical to a run without it.

The physical picture: rainfall that does not run off **infiltrates** and feeds a
water table, solved each step as an implicit Dupuit–Boussinesq head on the mesh
(draining to the sea, lakes and rivers). Where the water-table depth sits in a
shallow **capillary fringe**, an indurated crust (calcrete / silcrete / ferricrete,
treated generically) precipitates and **hardens the surface**, so crusted cells
erode more slowly, producing relief inversion and, when stratigraphy is on,
stacked duricrusts that are buried and later re-exposed. The duricrust runs
**soil-independent** by default and **couples to** :ref:`soil <surfproc>`
production when soil is tracked. A ``geochem`` sub-block adds a conservative
solute geochemistry (dissolution, transport, precipitation in the crust and
export down the rivers to the ocean). The equations are in the technical guide
(:ref:`groundwater`).

.. contents:: On this page
   :local:
   :depth: 1


Declaration
-----------

.. grid:: 1
    :padding: 3

    .. grid-item-card::

        **Declaration example**:

        .. code:: yaml

            groundwater:
                Ksat: 3.65e4          # hydraulic conductivity K_h (m/yr)
                specific_yield: 0.1   # S (drainable porosity)
                aquifer_base: 50.0    # z_bed depth below surface (m)
                bedrock_depth: 0.0    # permeable rock below lHbed (from_soil only)
                min_sat_thickness: 1.0
                infiltration: 0.3     # fraction of (rain − evap) recharging
                conserve_baseflow: True
                picard_its: 3
                seepage_passes: 4
                duricrust:
                    form_rate: 1.0e-4     # k_form (m/yr at Φ=Ψ=1)
                    max_thickness: 5.0    # duriH_max (m)
                    fringe_depth: 3.0     # d0, fringe centre below surface (m)
                    fringe_width: 2.0     # w, Gaussian half-width (m)
                    discharge_gate: False # crust only in gw discharge zones
                    supply_exp: 1.0       # p on (rain − evap) in the proxy
                    weather_Ea: 0.0       # Arrhenius activation energy (0 ⇒ off)
                    armor_max: 0.9        # max fractional K reduction (0..1)
                    armor_diffusion: False  # also armor hillslope Cd
                    break_rate: 1.0       # k_break per unit incision
                    decay_rate: 1.0e-6    # k_decay disequilibrium (1/yr)
                    weathering:
                        mode: proxy       # proxy | rate | prodsoil
                        C_eq: 1.0
                        Dw: 1.0
                        path_length: 20.0
                        weatherability: 1.0


Water table
-----------

The **water-table (hydrology) keys** are:

a. ``Ksat``: saturated hydraulic conductivity ``K_h`` (m/yr), one value for the whole domain.
b. ``specific_yield``: drainable porosity ``S`` (the storage coefficient linking recharge to head change), default ``0.1``.
c. ``aquifer_base``: depth of the impermeable base ``z_bed`` below the surface (m): a **scalar**, a per-vertex **map** ``[file, key]``, or the string ``from_soil`` (tie the base to the bedrock elevation ``z_bed = lHbed − bedrock_depth``: requires soil tracking; in a depositional basin the base deepens to the bottom of the porous sediment fill). Default ``50.0``.
d. ``bedrock_depth``: permeable weathered/fractured-rock thickness below ``lHbed`` (m), used **only** with ``aquifer_base: from_soil`` (default ``0``).
e. ``min_sat_thickness``: floor ``b_min`` on the saturated thickness so the transmissivity stays positive near the base (m, default ``1.0``).
f. ``infiltration``: fraction ``f_infil`` of ``max(0, rain − evap)`` that recharges the aquifer; a scalar or a per-vertex map ``[file, key]`` (default ``0.3``).
g. ``conserve_baseflow``: return the seepage discharge to the river network so total river discharge stays ``≈ rain − evap`` (default ``True``). The infiltrated recharge leaves surface runoff and is **re-injected as baseflow** at the seepage nodes (rivers become baseflow-fed); also writes the ``baseflow`` output.
h. ``lake_exchange``: opt-in lake ↔ aquifer **volume** coupling (default ``False``). When on, the signed across-bed groundwater flux debits/credits each lake's fill budget, a lake ringed by a higher water table gains groundwater, one ringed by a lower table leaks. Off ⇒ lakes are fixed-head only (unchanged). With evaporation, the exchange enters the lake's water balance: it decides with the inflow whether the lake spills or stays closed, and the level of a closed lake (see :ref:`flow_evaporation`).
i. ``subglacial_recharge``: fraction of the glacial meltwater (``iceMeltRiverL``) that infiltrates the aquifer where the ice melts out (default ``0``: under ice the rain path is gated off; this is the one recharge path allowed there).
j. ``fine_infil_factor``: multiplier on ``f_infil`` for the fine end-member (dual lithology); ``< 1`` makes clay/fine surfaces infiltrate less than coarse/sand (default ``1``: no lithology dependence).
k. ``infil_slope_ref``: reference slope for a ``f/(1 + slope/infil_slope_ref)`` reduction of infiltration on steep terrain (default ``0``: off; slope is the steepest-descent gradient).
l. ``picard_its`` / ``seepage_passes``: inner iteration counts for the unconfined non-linearity ``T(h)`` and the seepage free-boundary discovery (defaults ``3`` / ``4``).

.. important::

    ``aquifer_base: from_soil`` needs :ref:`soil production <surfproc>` tracked (it
    reads the bedrock elevation ``lHbed``). Without soil it falls back to the
    surface as the base with a warning.


Duricrust
---------

The nested **duricrust keys** (omit the ``duricrust:`` block for a water table with no crust) are:

a. ``form_rate``: crust formation rate ``k_form`` (m/yr at full favourability and supply).
b. ``max_thickness``: maximum crust thickness ``duriH_max`` (m); the induration degree is ``duriF = duriH/max_thickness``.
c. ``fringe_depth`` / ``fringe_width``: centre ``d0`` and Gaussian half-width ``w`` (m) of the capillary-fringe favourability band ``Φ`` on the water-table depth.
d. ``discharge_gate``: restrict crust formation to groundwater **discharge zones** (default ``False``). When ``True`` the favourability ``Φ`` is multiplied by ``G = (−∇·q)⁺ / ((−∇·q)⁺ + R)``, the fraction of the cell's upward discharge imported by **lateral convergence** of the groundwater flux (vs local recharge ``R``): ``G→1`` at a convergent valley floor / seepage face, ``G→0`` at a divergent recharge rise. This is the **absolute-accumulation** (lateral, valley/footslope ferricrete) style, as opposed to the default **relative-accumulation** style that indurates wherever the water table is shallow (the plateau/*bowal* cuirasses that blanket flat uplands). Off ⇒ weight ``1`` everywhere (unchanged).
e. ``supply_exp``: exponent ``p`` on ``(rain − evap)`` in the default climate proxy supply.
f. ``weather_Ea``: Arrhenius activation energy (J/mol) for an optional temperature scaling of the supply (``0`` ⇒ off; reuses the soil ``tempMap`` when present).
g. ``armor_max``: maximum fractional erodibility reduction (0–1); a fully indurated cell (``duriF = 1``) has its ``K`` multiplied by ``1 − armor_max`` (e.g. ``0.9`` ⇒ 10× more resistant).
h. ``armor_diffusion``: also armor the hillslope diffusivity ``Cd`` by the same factor (default ``False``).
i. ``break_rate`` / ``decay_rate``: breakdown per unit surface incision ``k_break`` and the slow disequilibrium decay ``k_decay`` (1/yr) away from the fringe.

The optional ``weathering:`` sub-block selects the **solute supply** ``Ψ`` feeding
formation: ``mode: proxy`` (default, climate/temperature stand-in), ``rate`` (an
explicit Maher–Chamberlain chemical-weathering rate driven by the recharge, keys
``C_eq``/``Dw``/``path_length``/``weatherability``), or ``prodsoil`` (reuse the
soil production rate). ``rate`` and ``prodsoil`` fall back to the proxy when soil
is off.


Solute geochemistry
-------------------

A further optional ``geochem:`` sub-block turns on **conservative solute
geochemistry**. One or more lumped **tracers** (a "species" = a solute
class such as carbonate, silica or iron) are, each step: **dissolved** on
subaerial land, **transported** down the groundwater flux, **precipitated** where
the water table sits in the capillary fringe (feeding the crust ``duriH``) and
the remainder **exported** to the rivers, a closed mass budget
(``dissolved = precipitated + exported``) with a dissolved flux to the ocean.
When ``geochem:`` is on the transported precipitation **replaces** the local
proxy/rate supply as the crust source, so the duricrust ``form_rate`` is **not
used**; the per-species ``precip_rate`` (below) governs crust growth instead.
Coupled multi-species aqueous equilibrium (speciation, pH, activity) is
deliberately **out of scope**, the lumped multi-tracer model is the
scale-appropriate choice at km / My.

A single tracer ships by default; several tracers give calcrete / silcrete /
ferricrete typing (the dominant one per node is the ``crust_type`` output). Full
declaration with sensible values:

.. code:: yaml

    groundwater:
        # ... hydrology + duricrust keys ...
        geochem:
            conserve: True                # keep/report the closed solute budget
            weatherability_from: source_class   # (optional) lithology label for the tables below
            lithology: [litho, rock_class]      # (optional) standalone per-vertex lithology map
            river_load: True              # (optional) route the export down the rivers
            marine_coupling: True         # (optional) accumulate the coastal delivery
            species:
              - name: carbonate           # tracer label (names the per-species outputs)
                weatherability: 1.0        # relative dissolution rate (scalar, map, or by-class)
                precip_rate: 1.0           # fringe-precipitation RATE (1/yr); scale down for large dt
                solid_volume: 1.0          # crust volume per unit precipitated solute
                c_sat: 1.0                 # saturation threshold (reserved; see note)
                river_decay: 0.0           # in-transit river loss (0 = conservative)
              - name: silica
                weatherability: 0.5
                precip_rate: 0.8
                river_decay: 0.3

**Block-level** ``geochem:`` keys:

.. list-table::
    :header-rows: 1
    :widths: 22 58 20

    * - Key
      - Meaning
      - Default
    * - ``species``
      - list of tracers (one mapping per species; keys below). Empty / absent ⇒ one default tracer.
      - one tracer
    * - ``conserve``
      - keep the per-step mass budget (``dissolved = precipitated + exported``) and report it.
      - ``True``
    * - ``weatherability_from``
      - lithology label for the ``weatherability_by_class`` tables: ``source_class`` (static bedrock) or ``surface_class`` (dynamic top stratigraphic layer). Needs ``provenance:``.
      - *(none)*
    * - ``lithology``
      - ``[file, key]`` per-vertex **integer** lithology map (an alternative label for the by-class tables, independent of provenance).
      - *(none)*
    * - ``river_load``
      - route the exported solute **down the drainage network** to the coast (per species), the river dissolved load ``riverSolute``.
      - ``False``
    * - ``marine_coupling``
      - accumulate the delivered coastal flux into a per-species marine reservoir + write ``marineSoluteInput``. Needs ``river_load``.
      - ``False``

**Per-species** keys (inside each ``species:`` entry). ``weatherability``,
``precip_rate`` and ``solid_volume`` are **relative, dimensionless** tuning knobs
,  what matters is their ratio between species, not the absolute value:

.. list-table::
    :header-rows: 1
    :widths: 20 60 20

    * - Key
      - Meaning
      - Default
    * - ``name``
      - tracer label; names the per-species outputs (e.g. ``solute_carbonate``).
      - ``solute<i>``
    * - ``weatherability``
      - relative **dissolution** rate (scales the local weathering supply for this tracer). A scalar, a per-vertex map ``[file, key]``, or a per-lithology ``weatherability_by_class`` list. ``0`` ⇒ this species does not weather here.
      - ``1.0``
    * - ``weatherability_by_class``
      - list of per-lithology weatherabilities gathered by the ``lithology`` map or the provenance class (see the three forms below).
      - *(none)*
    * - ``precip_rate``
      - **fringe-precipitation rate** ``k_p`` (per year), how fast the transported solute precipitates into crust in the capillary fringe. Precipitation is **self-limiting** (``·(1−duriH/max_thickness)``): a crust at ``max_thickness`` rejects further solute, which is exported instead. ``k_p`` is a **rate**, so scale it with the timestep, at a My ``dt`` a value near 1 saturates the crust to ``max_thickness`` in one step; use a small value (e.g. ``1e-6``) for gradual growth.
      - ``1.0``
    * - ``solid_volume``
      - crust **volume produced per unit precipitated solute** (a molar-volume-like factor turning precipitated mass into ``duriH``).
      - ``1.0``
    * - ``source_pool``
      - dissolvable **source-rock reservoir** per unit cell area (the weatherable rock mass). Dissolution debits it, so a pool too small for the run's weathering rate **exhausts** and the tracer goes inert (a one-time ``[gw] geochem: source pool exhausting`` warning is printed). Set it **large** (a thick saprolite source) to keep weathering *rate-limited* over long, fast-weathering runs, or small to model a finite, exhaustible weathering front (e.g. a wet-phase iron pulse).
      - ``1.0e6``
    * - ``c_sat``
      - saturation threshold (**reserved**, precipitation is currently a linear fringe sink; a hard ``c > c_sat`` gate is a future nonlinear refinement).
      - ``1.0``
    * - ``river_decay``
      - first-order **in-transit loss** along the river (in-channel precipitation / uptake), per drainage step. ``0`` = conservative; a small value (``0.1``–``0.5``) for a reactive species.
      - ``0.0``

**Spatial (lithology-driven) weatherability.** The ``weatherability`` may vary in
space so *lithology* controls which species each region yields (mafic rock →
Fe / silica, a carbonate platform → carbonate). Three forms, on top of the scalar
default:

.. code:: yaml

        geochem:
            lithology: [litho, rock_class]         # (b) per-vertex integer lithology map
            # weatherability_from: source_class    # (c) OR reuse the provenance regions
            species:
              - {name: carbonate, weatherability_by_class: [1.0, 0.0]}   # per rock class
              - {name: silica,    weatherability: [litho, sil_wab]}      # (a) per-vertex map

* **(a)** ``weatherability: [file, key]``: a **per-vertex map** for that species
  (loaded like ``infiltration``).
* **(b)** ``weatherability_by_class: [...]`` with a standalone ``lithology:
  [file, key]`` integer map, a **per-(class, species) table** gathered by that
  per-vertex lithology label (independent of provenance).
* **(c)** the same ``weatherability_by_class`` with ``weatherability_from:``
  a provenance class, ``source_class`` (**static** bedrock) or ``surface_class``
  (**dynamic**, the dominant provenance of the top stratigraphic layer,
  re-derived each step, so the weatherability tracks the rock actually exposed as
  erosion exhumes deeper layers or deposition buries the surface). Needs
  ``provenance:``; no extra input.

The label for (b)/(c) is the ``lithology:`` map when given, otherwise the
provenance class.


River dissolved load and the ocean
----------------------------------

**River dissolved load** (``river_load: True``). The exported (seepage / baseflow)
solute is routed **down the drainage network** to the shoreline, the dominant
natural pathway of weathering products to the sea, **per species** on the flow
matrix. ``riverSolute`` (m³/yr) grows downstream, is delivered at the coast, and
is trapped in closed continental basins (evaporite behaviour). Two refinements:

* per-species ``river_decay``: a first-order **in-transit loss** (in-channel
  precipitation / uptake; ``0`` = conservative);
* ``marine_coupling: True``: accumulate the delivered coastal flux into a
  per-species marine reservoir + write the per-node ``marineSoluteInput``
  (where weathering solute enters the sea).

Extract the per-basin dissolved / species flux at river mouths with
``gospl-catchment`` (see :ref:`running`).


Outputs
-------

Water table: ``recharge``, ``wtable``, ``wtdepth``, and
``baseflow`` (with ``conserve_baseflow``). Duricrust: ``duricrust``,
``induration``, ``Karmor``; with stratigraphy, the per-layer archive
``stratDuri`` (degree) plus ``stratCrustType`` / ``stratCrustSource`` (which
crust and from where), all shown by ``gospl-strata-volume``. Geochemistry:
``solute`` (concentration), ``soluteflux`` (dissolved export), ``crust_type``
(dominant former, several tracers), ``crust_source`` (with ``provenance:``),
``riverSolute`` (with ``river_load``) and ``marineSoluteInput`` (with
``marine_coupling``). With several tracers every aggregate is **also written
per species** (``solute_<name>``, ``soluteflux_<name>``, ``crust_<name>``,
``riverSolute_<name>``). See the technical
`groundwater documentation <../tech_guide/groundwater.html>`_
for the formulation.
