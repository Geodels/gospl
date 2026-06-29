.. _analyse_ref:

==========================
Post-processing tools
==========================

The :mod:`gospl.analyse` package and a couple of helpers in :mod:`gospl.tools`
provide the **post-processing / analysis** scripts that turn a finished goSPL
run (XDMF + per-partition HDF5) into publication products: regular-grid NetCDFs,
per-basin river profiles and outflow fluxes, stratigraphic sections / wells /
Wheeler diagrams, a 3-D stratigraphic wedge and source-to-sink provenance.

Each is importable as a function API **and** runnable as a console command (see
:ref:`running` for the command-line usage). They need the optional extras
``pip install gospl[analysis]`` (``numba``, ``geopandas``, ``netCDF4``,
``matplotlib``); the heavy dependencies are imported lazily so importing the
modules themselves is cheap.

.. contents::
    :local:
    :depth: 1


Regular-grid NetCDF + river profiles — ``gospl-grid``
-----------------------------------------------------

.. automodule:: analyse.gridexport

.. rubric:: Functions

.. autofunction:: analyse.gridexport.grid_export
.. autofunction:: analyse.gridexport.to_netcdf
.. autofunction:: analyse.gridexport.basin_rivers
.. autofunction:: analyse.gridexport.plot_long_profile
.. autofunction:: analyse.gridexport.plot_basin_map


Per-basin outflow fluxes — ``gospl-catchment``
----------------------------------------------

.. automodule:: analyse.catchment

.. rubric:: Functions

.. autofunction:: analyse.catchment.basin_outflow
.. autofunction:: analyse.catchment.catchment_flux


Stratigraphic sections, wells & Wheeler — ``gospl-section``
-----------------------------------------------------------

.. automodule:: analyse.stratasection

.. rubric:: Functions

.. autofunction:: analyse.stratasection.load_strata
.. autofunction:: analyse.stratasection.cross_section
.. autofunction:: analyse.stratasection.horizontal_slice
.. autofunction:: analyse.stratasection.synthetic_well
.. autofunction:: analyse.stratasection.well_panel
.. autofunction:: analyse.stratasection.wheeler


3-D stratigraphic wedge for ParaView — ``gospl-strata-volume``
--------------------------------------------------------------

.. automodule:: analyse.stratamesh

.. rubric:: Functions

.. autofunction:: analyse.stratamesh.discover_steps
.. autofunction:: analyse.stratamesh.detect_first_layer
.. autofunction:: analyse.stratamesh.build_partition
.. autofunction:: analyse.stratamesh.write_partition_h5
.. autofunction:: analyse.stratamesh.write_xdmf


Sediment provenance — ``gospl-provenance``
------------------------------------------

.. automodule:: analyse.provenance

.. rubric:: Classes

.. autoclass:: analyse.provenance.ProvenanceTracker
    :members:

.. autoclass:: analyse.provenance.GosplOutput
    :members:

.. rubric:: Functions

.. autofunction:: analyse.provenance.build_neighbours
.. autofunction:: analyse.provenance.steepest_receivers
.. autofunction:: analyse.provenance.downhill_edges
.. autofunction:: analyse.provenance.classes_from_shapefile
.. autofunction:: analyse.provenance.write_xdmf


ELA maps from paleo-temperature — ``gospl-ela``
-----------------------------------------------

.. automodule:: tools.ela_from_temperature

.. rubric:: Functions

.. autofunction:: tools.ela_from_temperature.derive_ela
