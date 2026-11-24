.. _gravity:

==============================
Marine gravity currents
==============================

.. grid:: 1
    :padding: 3

    .. grid-item-card::

        **Declaration example** (the regional margin of ``dual_lithology``, 615 m
        cells, ``dt`` = 5 kyr):

        .. code:: yaml

            gravity:
                fraction: 0.3         # or a flood block (below)
                runout: 1.e4          # coarse run-out length (m)
                runout_fine: 2.e4     # fine run-out length (m)
                bypass_depth: 50.     # m below sea level
                bypass_slope: 0.01    # gradient
                receivers: 2          # coarse network: narrow belts
                flow_exp: 4.
                receivers_fine: 6     # fine network: broad lobes
                flow_exp_fine: 1.
                network_smooth: 3     # passes
                slope_smooth: 3       # passes
                deposit_spread: 1500. # m
                coarse_factor: 0.3

        Adds a second pathway for the part of the river load that travels as a
        gravity current (hyperpycnal flow or turbidity current) instead of being
        spread by the clinoform scheme. That part is routed down the seafloor,
        deposits gradually with a run-out length, and builds coarse belts at the
        foot of the slope and finer lobes on the basin floor. Without the block,
        the marine deposition is unchanged (bitwise).

**Parameters**

All keys are optional. Lengths are in metres, gradients are dimensionless
(m/m), passes are numbers of neighbour-averaging sweeps.

.. list-table::
   :header-rows: 1
   :widths: 18 10 14 58

   * - Key
     - Unit
     - Default
     - Meaning
   * - ``fraction``
     - 0 to 1
     - absent
     - Fixed fraction of the river load reaching the sea that becomes a
       gravity current. When absent, the ``flood`` trigger is used.
   * - ``flood: c_crit``
     - kg/m³
     - 40
     - Hyperpycnal threshold concentration (Mulder and Syvitski, 1995).
   * - ``flood: rating_exp``
     - -
     - 1.8
     - Exponent b of the rating curve Qs ~ Q^b (must be > 1).
   * - ``flood: sigma``
     - -
     - 1.2
     - Log-standard deviation of daily discharge; larger for flashy rivers.
   * - ``flood: rho_s``
     - kg/m³
     - 2650
     - Sediment density used to turn the load into a concentration.
   * - ``coarse_factor``
     - 0 to 1
     - 0.3
     - Dual lithology only: the coarse load is diverted with this share of the
       fraction (the fine load with the full fraction). Gravity currents carry
       mostly the suspended load; 1 diverts as much coarse as fine.
   * - ``runout``
     - m
     - 5.e4
     - Distance over which the coarse (or single-fraction) deposit decays by
       a factor *e*. Must exceed the cell size.
   * - ``runout_fine``
     - m
     - 4 × ``runout``
     - Same for the fine fraction (dual lithology).
   * - ``bypass_depth``
     - m below sea level
     - 50
     - No deposition (or erosion) shallower than this; full rate at twice
       this depth. Also the accommodation limit: a gravity deposit never
       builds a cell above it. 0 turns both off.
   * - ``bypass_slope``
     - gradient
     - off
     - Deposition decreases on steeper seafloor (one half at this value,
       6 % at twice it), so flows cross the foresets.
   * - ``network``
     - ``raw`` / ``smoothed``
     - ``raw``
     - Surface the currents follow: the actual seafloor with its depressions
       filled (``raw``), or the clinoform network (``smoothed``, which
       ignores the five keys below).
   * - ``receivers``
     - 1 to 12
     - 6
     - Downslope neighbours sharing the coarse (or single) flux.
   * - ``flow_exp``
     - -
     - 2
     - Flux shared in proportion to slope^flow_exp: higher concentrates the
       flow on the steepest path.
   * - ``receivers_fine``, ``flow_exp_fine``
     - as above
     - coarse values
     - Network of the fine fraction (dual lithology).
   * - ``network_smooth``
     - passes
     - 0
     - Smoothing of the routing surface before the depressions are filled.
   * - ``slope_smooth``
     - passes
     - 0
     - Smoothing of the gradient used by ``bypass_slope`` and the erosion.
   * - ``deposit_spread``
     - m
     - 2 × mean cell size
     - Lateral spreading length of each step's deposit (volume conserving);
       0 turns it off.
   * - ``failure_slope``
     - gradient
     - absent (off)
     - Slope failure: sediment deposited during the run fails where the
       seafloor (deeper than ``bypass_depth``) is steeper than this, down to
       this gradient, and joins the gravity currents.
   * - ``max_failure``
     - m per step
     - 10
     - Cap of the failed thickness in a cell and step.
   * - ``erosion_length``
     - m
     - absent (off)
     - Expert option, for meshes of 200 m cells or finer (see the cautions).
       Slope erosion: the fraction of the flux picked up in a cell is
       (cell size / ``erosion_length``) × (S/S_e)^4 / (1 + (S/S_e)^4), at
       most 1. Must be much larger than the cell size (tens of cells).
   * - ``erosion_slope``
     - gradient
     - 0.02
     - S_e: gradient at which erosion reaches half its rate.
   * - ``max_erosion``
     - m per step
     - 10
     - Cap of the erosion in a cell and step.
   * - ``erosion_exp``, ``erosion_flux``
     - -, m³/yr
     - 1, 1.e4
     - Non-linear pick-up, multiplied by (Q / ``erosion_flux``)^(a - 1).

**Setting it up**

1. *Trigger.* Start with a fixed ``fraction`` (0.1 to 0.3) to see where the
   currents go, then switch to ``flood`` for a process-based fraction. The
   flood fraction depends on how concentrated the rivers are: with the default
   ``sigma`` it was below 1 % on the regional margin of the examples (dilute
   rivers) and about a fifth of the load on the global example, whose rivers
   carry large loads for their discharge. Raise ``sigma`` (about 2) for small
   steep catchments with flashy floods.
2. *Composition* (dual lithology). Gravity currents carry mostly the suspended
   load: the default ``coarse_factor`` of 0.3 diverts 0.3 x the fraction of the
   coarse load. With 1, the coarse load is diverted as fully as the fine one and
   builds thick coarse lobes at the foot of the slope.
3. *Run-out lengths.* Choose them from the distance between the shelf edge and
   the basin floor you want to fill: most of the coarse deposit lies within a
   few run-out lengths of the slope base, and the fine one reaches further (on
   the margin above, a 10 km coarse run-out left about 80 % of the coarse
   gravity deposit within 15 km of the slope base).
4. *Where they deposit.* Keep ``bypass_depth`` near the storm wave base (about
   50 m). On a regional mesh, read the seafloor gradients of your initial
   surface and set ``bypass_slope`` between the shelf and the foreset
   gradients, so the currents cross the foresets and deposit at their toe.
5. *Belts and lobes* (regional, dual lithology). Use a concentrated coarse
   network (``receivers: 2``, ``flow_exp: 4``) and a spreading fine one
   (``receivers_fine: 6``, ``flow_exp_fine: 1``), smooth the routing surface
   and the slope over a few cells (``network_smooth: 3``, ``slope_smooth: 3``).
   The deposit is spread over two cell sizes by default (``deposit_spread``).
6. *Global runs.* Use ``network: smoothed`` and leave the regional options
   (``bypass_slope``, the networks, the smoothing, the erosion) off; set the
   run-outs to several cell sizes (tens to hundreds of kilometres). This
   configuration has only been run on the test meshes so far.
7. *Slope failure* (optional). ``failure_slope`` remobilises sediment that
   builds steeper than a critical gradient, typically the upper slope and the
   foresets. Choose it from the gradients of your mesh (smoothed if
   ``slope_smooth`` is set): on the regional margin of the examples (615 m
   cells, foresets up to about 0.03), ``failure_slope: 0.02`` failed about
   0.95 km³ per 5 kyr step, a quarter of the river supply to the currents.
8. *Check the budget.* Run with ``-v`` (or ``--summary``): each step prints the
   supplied, eroded (including failed), deposited and exported volumes; the closure must stay at
   the solver tolerance (below 1e-6). ``gospl-inspect`` flags a budget that does
   not close.

**Set-up by resolution**

.. list-table::
   :header-rows: 1
   :widths: 16 28 56

   * - Mesh
     - What the currents represent
     - Settings
   * - Global, 10-50 km cells
     - Sediment delivered beyond the shelf to the continental rise and the
       abyssal plain, as a budget; no channels or lobes.
     - ``network: smoothed``; ``fraction`` or ``flood``; run-outs of several
       cells (tens to hundreds of km); leave ``bypass_slope``, the networks,
       the smoothing, failure and erosion off.
   * - Regional, 0.5-5 km cells
     - Coarse belts at the foot of the slope, fine lobes on the basin floor,
       slope failure on the foresets.
     - ``network: raw``; ``bypass_slope`` between the shelf and foreset
       gradients of the mesh; concentrated coarse network (``receivers: 2``,
       ``flow_exp: 4``), spreading fine one (``receivers_fine: 6``,
       ``flow_exp_fine: 1``); ``network_smooth``/``slope_smooth`` of a few
       passes; optional ``failure_slope`` from the mesh gradients.
   * - Local, 200 m cells or finer
     - As regional, with channels several cells wide.
     - As regional; slope erosion (``erosion_length`` of tens of cells,
       small ``max_erosion``) becomes meaningful (expert option, finer meshes
       not tested yet).

**Cautions**

- *Resolution.* Gradients depend on the mesh (a coarse mesh smears the slope),
  so ``bypass_slope``, ``erosion_slope`` and the network options must be chosen
  for each mesh. The smoothing keys are in passes, so their width scales with
  the cell size; ``deposit_spread`` and the run-outs are lengths and do not.
- *Slope erosion is an expert option: do not use it on meshes coarser than
  about 200 m.* It needs the channels resolved by the mesh and the time step:
  submarine channels are a few hundred metres to a few kilometres wide and
  need several cells across. On the 615 m regional mesh with 5 kyr steps it
  only roughened the toe of the slope (scours of tens of metres), and goSPL
  prints a warning when it is used on cells coarser than 200 m (an estimate:
  finer meshes have not been tested). Keep ``max_erosion`` small. ``erosion_length``
  must be many cell sizes: the pick-up compounds along a path, and with a
  length shorter than a cell the flux grows so fast that the erosion solve
  fails (``[flow] routing repair failed`` in the log) and no erosion is
  applied that step. ``erosion_exp``
  has only been tested up to 2.
- *Single-cell pathways.* A concentrated network (few receivers, high
  exponent) with ``deposit_spread: 0`` builds ridges one cell wide.
- *Slope failure* depends on the mesh gradients like ``bypass_slope``; it adds
  some cell-scale roughness on the slope (on the regional margin, the
  cell-scale roughness of the gravity deposit rose from 0.50 to 0.63 m), and
  only sediment deposited during the run can fail.
- *Domain walls.* Cells along a domain edge are half cells; their deposited
  fraction, set by the cell size, is smaller, so flows running along a wall
  deposit a little less.
- *Open boundaries export.* Currents reaching an open or fixed edge leave the
  domain with their load (reported as exported). Place the basin floor inside
  the domain, or use walls, if the distal lobes matter.
- *Global meshes.* Small steep coastal catchments (the most hyperpycnal
  rivers) are merged into larger outlets, so their gravity currents are
  underestimated. The seafloor depressions are filled on one processor (the
  whole surface is gathered on it and filled, then each processor receives
  its own nodes), as for the clinoform
  network that goSPL builds anyway. ``network: smoothed`` reuses that fill;
  ``raw`` adds a second one every step, which costs time and memory on meshes
  of millions of nodes: prefer ``smoothed`` there.
- *Calibration.* The module has been tested on one regional margin and on the
  test meshes. Compare the deposit with the clinoform wedge without the block
  (the gravity part moves sediment from the foresets to the basin floor) and
  check that the result does not depend on ``dt``.

**Outputs and cost**

The cumulative gravity deposit (net of any slope erosion) is written as
``gravDep`` (m). With a ``provenance`` block, the source classes travel with the
currents. With ``--summary``, each step's record has a ``gravity`` budget. On the
regional margin above (52,805 nodes, 4 processors) a step takes 1.25 s with the
block and 0.96 s without, about 30 % more: the routing solves take most of it
(about 0.2 s), the network about 0.04 s and the spreading about 0.04 s. Design
and tests: ``docs/DESIGN_GRAVITY_CURRENTS.md``.

Theory and numerics: :ref:`the gravity-current page of the technical guide
<tech_gravity>`.
