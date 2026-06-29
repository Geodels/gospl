.. _running:

==================
Running goSPL
==================

Once the :ref:`input file <inputfile>` is ready, a goSPL simulation can be
launched either from the command line or from a Python script — both run the
same forward model.

Command line
------------

Installing goSPL provides a ``gospl`` console command:

.. code-block:: bash

    # serial
    gospl -i input.yml

    # parallel (recommended) — choose any number of processes
    mpirun -np 8 gospl -i input.yml -v

Options:

- ``-i`` / ``--input`` *(required)* — the YAML input file.
- ``-v`` / ``--verbose`` — print per-step model progress (phase-grouped). The
  first line reports the goSPL version.
- ``--log`` — write the PETSc solver log summary.
- ``--profile`` — record a per-phase wall-clock profile (``profile.json``).
- ``--version`` — print the goSPL version and exit.

``--i`` / ``--v`` are accepted as aliases of ``--input`` / ``--verbose``. The
equivalent ``python -m gospl -i input.yml -v`` is available without a
console-script install. Note that **fast mode** is not a flag — set
``domain: fast: true`` in the input file.

Python script
-------------

The command above is a thin wrapper around the :class:`~gospl.model.Model`
class, so a run can equally be driven from Python (handy for parameter sweeps,
notebooks or coupling):

.. code-block:: python

    from gospl.model import Model as sim

    # Read the input file  (filename, verbose, showlog)
    model = sim(args.input, args.verbose, args.log)

    # Run the forward model
    model.runProcesses()

    # Clean up the PETSc objects
    model.destroy()

Launch it under MPI the same way as any mpi4py program::

    mpirun -np 8 python run_model.py -i input.yml

The number of processes used to **run** the model is independent of the number
used later for post-processing.


Post-processing & utility commands
----------------------------------

Installing goSPL also provides several command-line tools (each is equally
``python -m gospl.<module>`` if you haven't reinstalled to register the
console scripts). They read the model output in the ``h5`` directory; the
:ref:`output fields <outputs>` page describes what each variable means.

Stratigraphic volume for ParaView — ``gospl-strata-volume``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Stacks the recorded layers into a 3-D **wedge (triangular-prism) volume** and
writes an XDMF (per-partition HDF5 + one ``.xdmf``) for ParaView::

    gospl-strata-volume --h5dir myrun/h5 --field lithology     # or provenance
    mpirun -np 4 gospl-strata-volume --h5dir myrun/h5 --field provenance \
        --tout 100000                                          # ParaView time in yr

Key options: ``--outdir``/``--out`` (default ``strata``), ``--field
lithology|provenance``, ``--tout``/``--tstart`` (label time in years),
``--steps``, ``--first-layer``. Runs serially or under ``mpirun`` (partitions
split across ranks, independent of the run's processor count).

Per-cell fields attached to every wedge: ``thickness``, ``elevation``
(wedge mid-height) and ``layer`` (the recorded stratigraphic-layer index)
always; ``--field lithology`` adds coarse/fine thickness,
fine fraction and per-fraction porosity; ``--field provenance`` adds the
per-class volume fraction (``src_classN``), the ``dominant`` source class and
the per-layer ``porosity`` (plus ``phiFine`` for a dual-lithology run). An
eroded / pinched-out layer (zero thickness, hence no source) inherits the
``dominant`` source and composition of the cell directly below it in the
column, so it never renders as a no-source cell.

Publication sections, wells & Wheeler — ``gospl-section``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Vector (PDF/SVG) stratigraphic figures via matplotlib::

    # vertical cross-section coloured by lithology, with surface + basement lines
    gospl-section --h5dir myrun/h5 --mesh input/mesh.npz:v:c --kind cross \
        --along x --color-by lithology --out section.pdf

    gospl-section ... --kind slice   --z -50          # map at elevation z
    gospl-section ... --kind well    --xy 4.2e5,3.1e5  # synthetic borehole
    gospl-section ... --kind wheeler --color-by thickness \
        --strat-dt 1e4                                 # chronostratigraphic chart

    # depositional facies from the deposition water depth, exaggerated x10,
    # custom figure size, an interface line every 2 layers
    gospl-section ... --kind cross --color-by facies --vexag 10 \
        --figsize 12,3 --layer-lines 2

``--color-by`` ∈ {deposition, thickness, lithology, coarse, porosity, age,
provenance, **facies**}; ``--along x|y`` or ``--path x0,y0;x1,y1;...``. The
cross-section y-axis shows **true elevation**; ``--vexag`` applies a real
vertical exaggeration (data aspect, labels stay true); ``--figsize W,H`` sets
the figure size; ``--layer-lines N`` overlays a thin interface line every N
layers. **``facies``** classifies each layer by the **water depth at
deposition** (``sea_level − stratZ``): fluvial/deltaic plain (subaerial),
shoreface (0–20 m), distal offshore (20–50 m), upper slope (50–75 m), lower
slope (>75 m) — all tunable via ``--facies-depths`` (bin edges) and
``--facies-colors``; ``--color-by facies`` works for both ``cross`` and
``wheeler`` (the Wheeler shoreline overlay is the subaerial↔marine boundary,
consistent with the facies). ``--figsize`` applies to every ``--kind``. The
**time axis** of the Wheeler / ``age`` colouring is the **real simulation
time**: the per-layer interval is derived from the step's stratal display time
(the ``.xmf`` ``Time``) so the surface layer maps to that time —
``--strat-dt YR`` overrides it. The sea-level datum (section line / Wheeler
shoreline trajectory / facies depth reference) is **read from the simulation**
for the step by default (the step's ``.xmf`` ``sea`` constant);
``--sea-level FLOAT`` overrides it. ``--legend-loc`` positions the legend
(e.g. ``'upper left'``), ``--title-fontsize`` sets the title size, and
``--xlim MIN,MAX`` / ``--ylim MIN,MAX`` clip the cross-section / Wheeler to a
distance / elevation (or time) window (pinning a limit leaves the aspect auto,
so ``--vexag`` then has no effect). ``--figsize`` is honoured in the saved file
by default (pass
``--tight`` to crop to content instead); for the cross-section, ``--vexag``
exaggerates via the data aspect *without* overriding ``--figsize``. Horizontal
**distance is labelled in km and time in ky** (the data stay in m / yr). The
synthetic **well** draws its colour bar **horizontally at the base** (few ticks,
small font, so it stays readable for a tall narrow ``figsize`` like ``1,8``);
:func:`~gospl.analyse.stratasection.well_panel` draws **several wells on one
figure** with a shared colour scale and a single colour bar (titled ``well 1``,
``well 2``, … by default, or pass ``labels``).
The same functions are importable for
notebooks (with extra keyword args — ``figsize``, ``layer_lines``,
``facies_depths``, ``facies_colors``, ``facies_labels``, ``legend_loc``,
``title_fontsize``, ``xlim``, ``ylim``, well ``cbar_orientation`` / ``vmin`` /
``vmax`` / ``colorbar`` / ``ylim``, panel ``labels`` / ``vmin`` / ``vmax`` /
``ylim``)
(:func:`~gospl.analyse.stratasection.cross_section`,
:func:`~gospl.analyse.stratasection.horizontal_slice`,
:func:`~gospl.analyse.stratasection.synthetic_well`,
:func:`~gospl.analyse.stratasection.wheeler`,
:func:`~gospl.analyse.stratasection.well_panel`).

Regular-grid NetCDF for PyGMT / ArcGIS — ``gospl-grid``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rasterises the surface to a CF-NetCDF grid with drainage **basins**, **chi** and
drainage area; auto lon/lat for global meshes::

    gospl-grid --h5dir myrun/h5 --mesh input/mesh.npz:v:c --out surface.nc \
        --spacing 1000 --mn 0.5

Base level for catchments/chi defaults to the run's sea level (``--base-level``
overrides); ``--latlim`` crops the polar caps (default 89.9). The sea level
used is written into the NetCDF as the global attribute ``sea_level`` and a
scalar ``sea_level`` variable (and returned in the result dict as
``base_level``), so the grid is self-describing. The priority-flood-**filled**
elevation is also written (``filled``). Every NetCDF variable carries its
``units`` and a ``long_name`` definition (CF-style), so the file is
self-describing. The companion notebook API extracts /
plots per-basin river longitudinal profiles
(:func:`~gospl.analyse.gridexport.basin_rivers`,
:func:`~gospl.analyse.gridexport.plot_long_profile`,
:func:`~gospl.analyse.gridexport.plot_basin_map`). ``plot_long_profile`` plots
the **raw** elevation by default (``which='elev'`` — keeps real lakes /
depressions + interpolation roughness as small peaks); ``which='filled'`` plots
the hydrologically-conditioned elevation, which is **strictly monotonic**
upstream. ``plot_basin_map`` overlays the **sea-level coastline** (the
``elev == sea_level`` contour; defaults to the run's sea level) and takes a
``figsize``.

Per-basin outflow fluxes — ``gospl-catchment``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For every drainage **basin** of a gridded surface, extract the cell of maximum
water discharge and the cell of maximum sediment load — i.e. each basin's
**river-mouth (outflow) point** and its flux. It reads the ``gospl-grid``
NetCDF directly (it uses the same variable names — ``FA``, ``sedLoad``,
``basin``, ``lon``/``lat`` — falling back to the legacy ``flowDischarge`` /
``sedimentLoad`` / ``basinID`` names, so old files still work). The water flux
defaults to ``FA``; pass ``--flow-var fillFA`` to use the depression-filled
accumulation instead (routes the trunk discharge *through* lakes)::

    gospl-catchment -i index.csv -o flowsed

The index CSV has two columns ``time,netcdf`` (one row per step)::

    time,netcdf
    1,results/surface1.nc
    5,results/surface5.nc
    10,results/surface10.nc

For each step it writes ``flowsed/flow{time}.csv`` and ``flowsed/sed{time}.csv``
(columns ``basin,lon,lat,val``; ``val`` in m³/yr). Basins with ``<= --min-cells``
cells (default 10) are skipped. The per-basin maximum is a single vectorised
``lexsort`` (grouped arg-max), so a global 0.1° grid is processed in ~1 s per
step **serially** — the old MPI fan-out is no longer needed. From a notebook:
:func:`~gospl.analyse.catchment.catchment_flux` (batch → CSVs) and
:func:`~gospl.analyse.catchment.basin_outflow` (one file → two DataFrames).

Sediment provenance — ``gospl-provenance``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Source-to-sink attribution (per-basin / per-pixel source mix, transport
distance, optional Cu-fertility layer) → HDF5 + XDMF + CSV. See
:ref:`provenance <provenance>` for the inputs and ``--cu-weights``.

ELA maps from paleo-temperature — ``gospl-ela``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A **pre-processing** helper that derives the glacier-geometry maps
(``hela``/``hice``) for the ``ice`` block from a per-vertex temperature map by
inverting the lapse rate::

    gospl-ela --temperature t2m.npz --t-key t2m --lapse 0.0065 --t-ela -5 \
        --elevation mesh.npz --z-key z --band 800 --out glaciers_0Ma

Then list the output in the ``ice.glaciers`` time series (``hela: ['glaciers_0Ma',
'hela']`` etc.; see :ref:`surfproc`). ``--reference surface|sealevel`` selects
whether the temperature map is at the surface or reduced to sea level.

.. note::

   The console commands (``gospl``, ``gospl-strata-volume``, ``gospl-section``,
   ``gospl-grid``, ``gospl-catchment``, ``gospl-provenance``, ``gospl-ela``)
   appear after installing
   goSPL. Until then, use the equivalent ``python -m gospl.<module>`` form. The
   analysis tools need the optional extras: ``pip install gospl[analysis]``
   (numba, geopandas, netCDF4, matplotlib).

From a Jupyter notebook (import the API)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Every tool is also importable, so it can be driven and the figures shown inline
in a notebook (no CLI needed):

.. code-block:: python

    # --- stratigraphic sections / wells / Wheeler (matplotlib, inline) ---
    from gospl.analyse.stratasection import (
        load_strata, cross_section, horizontal_slice, synthetic_well, wheeler)

    data = load_strata("myrun/h5", "input/mesh.npz", step=20)
    cross_section(data, kind="x", color_by="lithology", vexag=20, sea_level=0)
    synthetic_well(data, 4.2e5, 3.1e5, color_by="porosity")
    wheeler(data, kind="x", color_by="thickness", sea_level=0, dt=1e4)
    horizontal_slice(data, z=-50.0, color_by="age")

    # --- gridded NetCDF + per-basin river profiles (PyGMT/ArcGIS + matplotlib) ---
    from gospl.analyse.gridexport import (
        grid_export, to_netcdf, basin_rivers, plot_long_profile, plot_basin_map)

    g = grid_export("myrun/h5", "input/mesh.npz", step=20, spacing=1000.)
    to_netcdf(g, "surface.nc")                       # for PyGMT / ArcGIS
    riv = basin_rivers(g, basin_id=g["main_basin"], area_threshold=5e6)
    plot_long_profile(riv); plot_basin_map(g, riv)

    # --- ELA maps from temperature (pre-processing) ---
    import numpy as np
    from gospl.tools.ela_from_temperature import derive_ela
    hela, hice = derive_ela(np.load("t2m.npz")["t2m"], 0.0065, -5.0, 800,
                            reference="surface", elevation=np.load("mesh.npz")["z"])

Each plotting call returns the matplotlib ``Axes``/``Figure`` (or PyGMT-ready
NetCDF) so you can tune it for a paper. The wedge volume (``gospl-strata-volume``)
is written for ParaView rather than inline rendering.
