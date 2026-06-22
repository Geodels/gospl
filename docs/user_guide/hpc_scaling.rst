.. _hpc_scaling:

==============================
HPC scaling & performance
==============================

This page reports **honest** strong-scaling measurements for goSPL on
production, global-scale meshes. "Strong scaling" means the *same* problem is
solved on an increasing number of MPI ranks (cores): ideally the wall-clock time
halves when the core count doubles. We report both the **speedup** and the
**parallel efficiency**, and we are explicit about where — and why — efficiency
falls away.

All runs exercise the full model (flow accumulation, river incision/deposition,
continental and marine sediment routing, hillslope diffusion, flexural isostasy,
tectonic forcing) so the numbers reflect a realistic workload, not a single
kernel.


Global mesh resolutions
-----------------------

The benchmarks use three icosahedral-derived global meshes (variable resolution;
the edge lengths below are in km):

.. list-table::
   :header-rows: 1
   :widths: 18 16 16 16 22

   * - Cell width (km)
     - Edge min (km)
     - Edge max (km)
     - Edge mean (km)
     - Nodes
   * - 5
     - 1.1
     - 4.5
     - 2.8
     - 23,632,811
   * - 8
     - 1.8
     - 7.2
     - 4.6
     - 9,236,387
   * - 10
     - 2.2
     - 8.9
     - 5.7
     - 5,912,778


How the benchmark is run
------------------------

* **Hardware**: NCI *Gadi* (48-core nodes; counts above 48 span whole nodes).
  The harness lives in ``scripts/scaling/`` (see its ``README.md``); each rank
  count is a separate, optionally chained job (``submit_sweep.sh`` → ``gadi.pbs``).
* **What is timed**: a fixed number of model steps with the wall-clock phase
  profiler enabled and **I/O disabled** (``--io off``), so the numbers reflect
  *compute + communication*, not HDF5 bandwidth.
* **Baseline**: a global mesh is **memory-infeasible at very low core counts**
  (the 10 km mesh needs ~6.6 GB on the heaviest rank and cannot start below
  ~16 ranks; finer meshes need more). There is therefore **no 1-CPU run** to
  normalise against. We follow the standard convention and **baseline to the
  smallest feasible rank count** :math:`p_0` in each sweep:

  .. math::

     S(p) = \frac{T(p_0)}{T(p)}, \qquad E(p) = S(p)\,\frac{p_0}{p}, \qquad
     \text{ideal: } S = p/p_0,\ E = 1.

  When reading the figure, the dashed line is the *ideal* :math:`p/p_0`, and
  efficiency ``1.0`` is perfect scaling relative to :math:`p_0` — **not** to a
  single core.


Results
-------

.. figure:: ../images/scaling_hpc.png
   :width: 100%
   :align: center
   :alt: goSPL strong-scaling speedup and parallel efficiency vs MPI ranks.

   Strong-scaling **speedup** (a) and **parallel efficiency** (b) versus MPI
   rank count, for the global meshes above. Markers are measured points; dashed
   lines are the ideal (each sweep baselined to its own smallest rank count, so
   the 10 km and 8 km curves have different ideal references). *(The 5 km curve
   is being collected and will be added here.)*

**10 km mesh (5.9 M nodes), baseline = 16 ranks:**

.. list-table::
   :header-rows: 1
   :widths: 12 16 14 16 20

   * - Ranks
     - Wall (s)
     - Speedup
     - Efficiency
     - RSS / rank (GB)
   * - 16
     - 521.3
     - 1.00
     - 1.00
     - 6.6
   * - 24
     - 385.6
     - 1.35
     - 0.90
     - 6.6
   * - 48
     - 187.9
     - 2.77
     - 0.92
     - 6.6
   * - 96
     - 98.6
     - 5.29
     - 0.88
     - 6.6
   * - 144
     - 69.9
     - 7.46
     - 0.83
     - 6.7
   * - 192
     - 54.9
     - 9.49
     - 0.79
     - 6.6
   * - 240
     - 46.9
     - 11.11
     - 0.74
     - 6.6
   * - 288
     - 44.4
     - 11.73
     - 0.65
     - 6.6

**8 km mesh (9.2 M nodes), baseline = 24 ranks:**

.. list-table::
   :header-rows: 1
   :widths: 12 16 14 16 20

   * - Ranks
     - Wall (s)
     - Speedup
     - Efficiency
     - RSS / rank (GB)
   * - 24
     - 800.7
     - 1.00
     - 1.00
     - 9.9
   * - 48
     - 366.4
     - 2.19
     - 1.09
     - 9.7
   * - 96
     - 194.7
     - 4.11
     - 1.03
     - 9.8
   * - 144
     - 150.0
     - 5.34
     - 0.89
     - 9.8
   * - 192
     - 116.3
     - 6.89
     - 0.86
     - 9.8
   * - 240
     - 96.0
     - 8.34
     - 0.83
     - 9.8
   * - 288
     - 71.4
     - 11.22
     - 0.93
     - 9.8


Reading the results
-------------------

* **Near-ideal to ~100 ranks.** Efficiency holds around 0.88–0.92 out to 96
  ranks (a 5.3× speedup), i.e. the compute-heavy phases — sediment routing
  (``sed``), marine deposition (``sea``), flow accumulation, erosion — all
  parallelise well and are perfectly load-balanced (per-phase imbalance ≈ 1.00).
* **Gentle decline beyond that, by Amdahl's law.** Efficiency eases to 0.74 at
  240 ranks and 0.65 at 288. This is **not** a load-balance problem; it is the
  growing weight of the small *serial* sections as the parallel work shrinks:

  - **flexure** (global spherical-harmonic solve, ``pyshtools``) is a flat
    ~4 s floor at every rank count;
  - the **pit-filling spillover graph** (Barnes-2016 master solve, serial on
    rank 0 by design — its cost scales with the partition *perimeter*) is a flat
    ~8 s floor and grows slightly with rank count.

  Together these set the practical ceiling.
* **Practical sweet spot ≈ 240 ranks** for the 10 km mesh: 11.1× speedup at 74 %
  efficiency. Going to 288 buys almost nothing (11.1× → 11.7×, efficiency
  0.74 → 0.65) — past here you are mostly paying the serial floors.
* **Memory is the lower bound, not a per-rank growth.** Peak RSS on the heaviest
  rank stays flat with rank count — ~6.6 GB (10 km) / ~9.9 GB (8 km) — because
  the dominant arrays are replicated, mesh-sized (``mpoints``) globals that do
  not decompose. That flat floor is what sets the *smallest* feasible rank count
  (the 8 km mesh cannot start below ~24 ranks); the aggregate footprint grows
  roughly linearly with rank count.
* **The 8 km mesh scales superlinearly to ~96 ranks** (efficiency 1.09 at 48,
  1.03 at 96 — *above* the ideal line). This is genuine, not a measurement
  artefact: at the 24-rank baseline the per-rank working set is large and spills
  cache / saturates memory bandwidth, so the baseline is *slow*; adding ranks
  shrinks each rank's slice until it fits, and per-core throughput rises. The
  effect fades once the slices are cache-resident (≥144 ranks), after which the
  same Amdahl decline takes over (0.89 → 0.83 to 240, with 288 a favourable
  outlier at 0.93). We report efficiency **uncapped** — superlinear points are
  shown as measured, not clipped to 1.0.

.. note::

   Efficiency is reported relative to each mesh's smallest feasible run
   (:math:`p_0` = 16 ranks at 10 km, 24 at 8 km), **not** a single core. We
   deliberately do not extrapolate a virtual 1-core time to inflate the speedup.
   One consequence: because the baseline itself is memory-bandwidth bound,
   efficiency can legitimately exceed 1.0 at intermediate rank counts (see the
   8 km mesh) — a real super-linear effect we leave un-clipped rather than hide.


Reproducing
-----------

The summary CSVs and the figure script live in
``docs/user_guide/scaling/``. To regenerate the figure (e.g. after adding the
8 km / 5 km sweeps as ``scaling_8km.csv`` / ``scaling_5km.csv``)::

    python docs/user_guide/scaling/make_scaling_figure.py

To run a sweep and produce a CSV, see ``scripts/scaling/README.md`` (the
``submit_sweep.sh`` → ``analyze_scaling.py`` → ``plot_scaling.py`` workflow).
