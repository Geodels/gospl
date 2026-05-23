.. _voro_ref:

================
Class VoroBuild
================

Voronoi helper used by :class:`mesher.unstructuredmesh.UnstMesh` to build the Centroidal Voronoi Tessellation underlying the finite-volume discretisation. Most users do not interact with it directly — :meth:`~mesher.unstructuredmesh.UnstMesh._meshStructure` drives it.

.. autoclass:: mesher.meshfunc.VoroBuild

   .. rubric:: Initialise

   .. autosummary::

      ~VoroBuild.initVoronoi

   .. rubric:: Public Methods

   .. autosummary::

      ~VoroBuild.angles
      ~VoroBuild.cell_barycenters
      ~VoroBuild.cell_centroids
      ~VoroBuild.cell_circumcenters
      ~VoroBuild.cell_partitions
      ~VoroBuild.ce_ratios_per_interior_edge
      ~VoroBuild.circumradius
      ~VoroBuild.compute_curl
      ~VoroBuild.control_volume_centroids
      ~VoroBuild.control_volumes
      ~VoroBuild.create_edges
      ~VoroBuild.edge_gid_to_edge_list
      ~VoroBuild.edge_lengths
      ~VoroBuild.edges_cells
      ~VoroBuild.face_partitions
      ~VoroBuild.inradius
      ~VoroBuild.is_boundary_facet
      ~VoroBuild.is_boundary_node
      ~VoroBuild.is_interior_node
      ~VoroBuild.mark_boundary
      ~VoroBuild.signed_cell_areas
      ~VoroBuild.surface_areas
      ~VoroBuild.triangle_quality

Public functions
---------------------

.. automethod:: mesher.meshfunc.VoroBuild.initVoronoi
.. automethod:: mesher.meshfunc.VoroBuild.angles
.. automethod:: mesher.meshfunc.VoroBuild.cell_barycenters
.. automethod:: mesher.meshfunc.VoroBuild.cell_centroids
.. automethod:: mesher.meshfunc.VoroBuild.cell_circumcenters
.. automethod:: mesher.meshfunc.VoroBuild.cell_partitions
.. automethod:: mesher.meshfunc.VoroBuild.ce_ratios_per_interior_edge
.. automethod:: mesher.meshfunc.VoroBuild.circumradius
.. automethod:: mesher.meshfunc.VoroBuild.compute_curl
.. automethod:: mesher.meshfunc.VoroBuild.control_volume_centroids
.. automethod:: mesher.meshfunc.VoroBuild.control_volumes
.. automethod:: mesher.meshfunc.VoroBuild.create_edges
.. automethod:: mesher.meshfunc.VoroBuild.edge_gid_to_edge_list
.. automethod:: mesher.meshfunc.VoroBuild.edge_lengths
.. automethod:: mesher.meshfunc.VoroBuild.edges_cells
.. automethod:: mesher.meshfunc.VoroBuild.face_partitions
.. automethod:: mesher.meshfunc.VoroBuild.inradius
.. automethod:: mesher.meshfunc.VoroBuild.is_boundary_facet
.. automethod:: mesher.meshfunc.VoroBuild.is_boundary_node
.. automethod:: mesher.meshfunc.VoroBuild.is_interior_node
.. automethod:: mesher.meshfunc.VoroBuild.mark_boundary
.. automethod:: mesher.meshfunc.VoroBuild.signed_cell_areas
.. automethod:: mesher.meshfunc.VoroBuild.surface_areas
.. automethod:: mesher.meshfunc.VoroBuild.triangle_quality
