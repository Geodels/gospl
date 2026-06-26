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
        --along x --color-by lithology --sea-level 0 --out section.pdf

    gospl-section ... --kind slice   --z -50          # map at elevation z
    gospl-section ... --kind well    --xy 4.2e5,3.1e5  # synthetic borehole
    gospl-section ... --kind wheeler --color-by thickness --sea-level 0 \
        --strat-dt 1e4                                 # chronostratigraphic chart

``--color-by`` ∈ {deposition, thickness, lithology, coarse, porosity, age,
provenance}; ``--along x|y`` or ``--path x0,y0;x1,y1;...``; ``--vexag`` vertical
exaggeration. The same functions are importable for notebooks
(:func:`~gospl.analyse.stratasection.cross_section`,
:func:`~gospl.analyse.stratasection.horizontal_slice`,
:func:`~gospl.analyse.stratasection.synthetic_well`,
:func:`~gospl.analyse.stratasection.wheeler`).

Regular-grid NetCDF for PyGMT / ArcGIS — ``gospl-grid``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rasterises the surface to a CF-NetCDF grid with drainage **basins**, **chi** and
drainage area; auto lon/lat for global meshes::

    gospl-grid --h5dir myrun/h5 --mesh input/mesh.npz:v:c --out surface.nc \
        --spacing 1000 --mn 0.5

Base level for catchments/chi defaults to the run's sea level (``--base-level``
overrides); ``--latlim`` crops the polar caps (default 89.9). The companion
notebook API extracts/plots per-basin river longitudinal profiles
(:func:`~gospl.analyse.gridexport.basin_rivers`,
:func:`~gospl.analyse.gridexport.plot_long_profile`,
:func:`~gospl.analyse.gridexport.plot_basin_map`).

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
   ``gospl-grid``, ``gospl-provenance``, ``gospl-ela``) appear after installing
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
