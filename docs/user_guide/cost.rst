.. _option_cost:

==========================================
Cost of the options: time and memory
==========================================

Every process you switch on makes a step longer and holds more memory. This
page gives rough numbers to plan a run: start from a basic model, add the
processes you need, and add up their costs. They are measured, not
modelled, on one regional mesh; read them as orders of magnitude.

.. contents:: On this page
   :local:
   :depth: 1


How it was measured
-------------------

- **Mesh**: the 200 x 100 km margin of the goSPL-examples
  ``mixed_carbonate_clastic`` example, 47,672 nodes (~700 m cells), with
  highlands, two large rivers, a shelf and a slope.
- **Base run**: linear stream power law with deposition (``spl: K, d``),
  linear hillslope and marine diffusion, marine deposition, a sea-level
  curve and a rainfall map, multiple flow directions; no stratigraphy.
- **Each option** is added alone to the base (or to its prerequisite, see the
  *added to* column), with the input blocks of the examples.
- **10 steps** of 5 kyr; the time of a step is the median over the steps
  after the first (the first sets up the solvers). The *own phases* column
  is the time spent in the option's own part of the step (the phases of the
  ``--summary`` record: ``groundwater``, ``carbonate``, ...).
- **Machine**: Apple M4 Max, macOS, goSPL 2026.7.14 (branch with the gravity
  currents, carbonates and karst), PETSc 3.21; one and four MPI ranks.

Repeated runs of the base agree to about 3 %. The script is
``scripts/option_cost.py`` (see :ref:`cost_measure` to run it on your mesh);
the raw numbers are in ``docs/user_guide/scaling/option_cost_np1.json`` and
``option_cost_np4.json``.


Time per step
-------------

The base step takes 2.6 s on one rank and 0.71 s on four. Values are
percentages of the base step; *own phases* is on four ranks.

.. list-table::
   :header-rows: 1
   :widths: 30 14 10 10 10 26

   * - Option
     - added to
     - 1 rank
     - 4 ranks
     - own phases
     - notes
   * - base
     -
     - 100
     - 100
     -
     - mostly the marine deposition (about two thirds) and the flow accumulation (about a third)
   * - **Erosion and hillslope laws**
     -
     -
     -
     -
     -
   * - non-linear SPL (``n: 1.5``)
     - base
     - 50
     - 58
     - 0
     - the slope exponent changes the relief, and the marine step does less work
   * - transport-limited SPL (``G: 0.5``)
     - base
     - 104
     - 124
     - 25
     -
   * - non-linear hillslope (``hillslopeSc``)
     - base
     - 101
     - 111
     - 0
     -
   * - soil / regolith (``soil:``)
     - base
     - 179
     - 173
     - 78
     - soil-aware SPL and soil diffusion (SNES, TS)
   * - **Climate, ice and the lithosphere**
     -
     -
     -
     -
     -
   * - flexure (``flexure: fem``)
     - base
     - 99
     - 100
     - 1
     - cached factorised operator
   * - orographic rain (``orography:``)
     - base
     - 45
     - 53
     - 2
     - lowers the rain here, so the flow and marine steps do less work
   * - diagnostic glacial model (``ice:``)
     - base
     - 95
     - 105
     - 8
     -
   * - **Stratigraphy**
     -
     -
     -
     -
     -
   * - stratigraphy (``strat`` interval)
     - base
     - 97
     - 103
     - 2
     - one layer per step; memory grows with the layers (below)
   * - dual lithology (``strata: dual``)
     - stratigraphy
     - 106
     - 115
     - 0
     - a second (mud) routing in the cascades
   * - provenance tracers (3 classes)
     - dual
     - 145
     - 159
     - 0
     - one routed column per class
   * - **Groundwater and geochemistry**
     -
     -
     -
     -
     -
   * - groundwater (water table)
     - base
     - 109
     - 134
     - 29
     - AMG head solve
   * - + duricrust
     - groundwater
     - 109
     - 130
     - 28
     -
   * - + geochemistry (2 species, river load)
     - duricrust
     - 570
     - 170
     - 95
     - direct solve per species: see the note below
   * - **Marine processes**
     -
     -
     -
     -
     -
   * - gravity currents (``gravity:``)
     - dual
     - 118
     - 127
     - (marine)
     - +12 % over dual
   * - pelagic rain (``pelagic:``)
     - dual
     - 148
     - 167
     - 18
     - +42 to +52 % over dual (carbonate and clay lithologies added)
   * - carbonate production (2 factories)
     - dual
     - 174
     - 182
     - 2
     - +68 % over dual: three more lithologies in the marine routing
   * - + open-water index (``exposure``)
     - carbonate
     - 189
     - 192
     - 4
     - one screened-Poisson solve
   * - + turbidity control
     - carbonate
     - 173
     - 188
     - 11
     -
   * - + waves (``waves:``, one direction)
     - carbonate
     - 174
     - 187
     - 9
     - one upwind solve per wave direction
   * - + coastal erosion
     - carbonate
     - 174
     - 188
     - 2
     -
   * - + karst (``karst:``)
     - carbonate
     - 172
     - 182
     - 1
     -

**How to read it.** The *own phases* column is the cost of the process
itself: it adds up. The whole-step columns also contain how the option
changes the landscape the other processes work on, which can make a run
cheaper (orographic rain and the non-linear SPL lower the marine work on
this margin) or dearer; that part depends on your model. The main costs:

- **What multiplies the routed columns**: every lithology and every
  provenance class is routed through the cascades and spread by the marine
  diffusion with the bulk sediment. Dual lithology adds one column,
  carbonates three (``carb_sand``, ``carb_mud``, ``reef``), provenance one per
  class: each adds about 10 to 20 % of the step.
- **Implicit solves of their own**: soil (+75 %), the groundwater head
  (+30 %) and its geochemistry (another +40 % on four ranks, one direct
  solve per species).
- **Cheap**: flexure (the operator is factorised once), ice, stratigraphy
  itself, the carbonate factories, waves, coastal erosion and karst (a few
  per cent each).

**The geochemistry on one rank.** Its transport is a direct factorisation
per species: with MPI it uses MUMPS, on one rank PETSc's own LU, which is
much slower on this matrix (570 % of the base step against 170 % on four
ranks). Asking for MUMPS on one rank brings it to 114 %:

.. code:: bash

    PETSC_OPTIONS="-gw_solute_pc_factor_mat_solver_type mumps" gospl -i input.yml

**Set-up time** (building the model) was under a second for every option
on this mesh; on large meshes it is dominated by the mesh distribution, see
:ref:`hpc_scaling`.


Memory
------

Each rank first holds about **560 MB** of libraries (PETSc, MPI, NumPy,
SciPy, VTK), whatever the mesh. The base model then adds, per node it holds,
a few kB (the mesh, the flow and diffusion matrices and their solvers: 5 to
8 kB per node, see the memory note in the developer guide), which is about
250 MB per rank here at four ranks.

On this mesh the measured peak memory per rank changed by -25 % to +20 %
between options, within the scatter of the peaks, which are set by
transient solver work (the marine diffusion, the pit filling) that depends
on the landscape. Only the geochemistry raised the peak systematically
(+19 % on one rank, its factorisations). So the options' own memory is small
next to the base model, **except the stratigraphy**, which allocates every
layer of the run at the start:

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Stratigraphy
     - bytes per node and per layer
   * - layers (thickness, elevation, porosity, erodibility)
     - 32
   * - dual lithology (mud pile and its porosity)
     - +16
   * - each further lithology (carbonates: three)
     - +16 each
   * - provenance tracers
     - +8 per class (+1 class with carbonates or pelagic rain)
   * - duricrust record, crust type, crust source (groundwater)
     - +8 each

The number of layers is the run length over the ``strat`` interval. A
1-million-node model with 100 layers, dual lithology and carbonates holds
1e6 x 100 x (48 + 48) B = 9.6 GB of stratigraphy, shared among the ranks:
more than the rest of the model. With many layers, choose ``strat`` for the
resolution you need in the record, not the time step.

Other per-node state of the options is a few arrays each (8 B per node per
array): negligible next to the base model.


Building up a run
-----------------

Rough times relative to the base step, from the table (four ranks):

- erosion, hillslope and marine deposition: **1x**
- + stratigraphy and dual lithology: **~1.15x**
- + provenance (3 classes): **~1.6x**
- + carbonates, waves, coastal erosion and karst (on dual): **~1.9x**
- + groundwater, duricrust and geochemistry: **+0.7x**
- + soil: **+0.7x**

A run with everything on is then about three times the base step, and its
memory is set by the base model and the stratigraphy.


.. _cost_measure:

Measuring it on your mesh
-------------------------

The numbers depend on the mesh, the resolution, the rank count and the
landscape. ``scripts/option_cost.py`` measures them for your own input:

.. code:: bash

    python scripts/option_cost.py my_input.yml -n 4 --steps 10 \
        --workdir cost_runs --json cost.json
    python scripts/option_cost.py --report cost_runs     # the table again

It keeps the core blocks of your input (domain, time, spl, diffusion, sea,
climate), runs the base and each option for ``--steps`` steps with
``mpirun -n``, and prints the median step time, the option's own phases and
the peak memory per rank, relative to the base. ``--only`` runs a subset
(``--only flexure,gw_geochem``); edit the ``VARIANTS`` list in the script to
add your own blocks. For the memory by stage of the set-up and the step, use
``scripts/mem_profile.py``; for the time by sub-stage,
``scripts/stage_profile.py``.
