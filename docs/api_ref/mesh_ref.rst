.. _mesh_ref:

==============
Class UnstMesh
==============

.. autoclass:: mesher.unstructuredmesh.UnstMesh

   .. rubric:: Initialise

   .. autosummary::

      ~UnstMesh.__init__

   .. rubric:: Public Methods

   .. autosummary::

      ~UnstMesh.applyForces
      ~UnstMesh.applyTectonics
      ~UnstMesh.destroy_DMPlex

   .. rubric:: Private Methods

   .. autosummary::

      ~UnstMesh._buildMesh
      ~UnstMesh._generateVTKmesh
      ~UnstMesh._get_boundary
      ~UnstMesh._meshfrom_cell_list
      ~UnstMesh._meshStructure
      ~UnstMesh._readErosionDeposition
      ~UnstMesh._set_DMPlex_boundary_points
      ~UnstMesh._updateRain
      ~UnstMesh._updateEroFactor
      ~UnstMesh._xyz2lonlat

Public functions
---------------------

.. automethod:: mesher.unstructuredmesh.UnstMesh.applyForces
.. automethod:: mesher.unstructuredmesh.UnstMesh.applyTectonics
.. automethod:: mesher.unstructuredmesh.UnstMesh.destroy_DMPlex


Private functions
---------------------

.. automethod:: mesher.unstructuredmesh.UnstMesh._buildMesh
.. automethod:: mesher.unstructuredmesh.UnstMesh._generateVTKmesh
.. automethod:: mesher.unstructuredmesh.UnstMesh._get_boundary
.. automethod:: mesher.unstructuredmesh.UnstMesh._meshfrom_cell_list
.. automethod:: mesher.unstructuredmesh.UnstMesh._meshStructure
.. automethod:: mesher.unstructuredmesh.UnstMesh._readErosionDeposition
.. automethod:: mesher.unstructuredmesh.UnstMesh._set_DMPlex_boundary_points
.. automethod:: mesher.unstructuredmesh.UnstMesh._updateRain
.. automethod:: mesher.unstructuredmesh.UnstMesh._updateEroFactor
.. automethod:: mesher.unstructuredmesh.UnstMesh._xyz2lonlat
