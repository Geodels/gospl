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
      ~UnstMesh.initExtForce
      ~UnstMesh.reInitialiseElev
      ~UnstMesh.updatePaleomap

   .. rubric:: Private Methods

   .. autosummary::

      ~UnstMesh._buildMesh
      ~UnstMesh._generateVTKmesh
      ~UnstMesh._localStrat
      ~UnstMesh._meshAdvectorSphere
      ~UnstMesh._meshUpliftSubsidence
      ~UnstMesh._meshfrom_cell_list
      ~UnstMesh._readErosionDeposition
      ~UnstMesh._readStratLayers
      ~UnstMesh._stratalRecord
      ~UnstMesh._updateDispStrata
      ~UnstMesh._updateRain
      ~UnstMesh._updateTectonics
      ~UnstMesh._xyz2lonlat

Public functions
---------------------

.. automethod:: mesher.unstructuredmesh.UnstMesh.applyForces
.. automethod:: mesher.unstructuredmesh.UnstMesh.applyTectonics
.. automethod:: mesher.unstructuredmesh.UnstMesh.destroy_DMPlex
.. automethod:: mesher.unstructuredmesh.UnstMesh.initExtForce
.. automethod:: mesher.unstructuredmesh.UnstMesh.reInitialiseElev
.. automethod:: mesher.unstructuredmesh.UnstMesh.updatePaleomap


Private functions
---------------------

.. automethod:: mesher.unstructuredmesh.UnstMesh._buildMesh
.. automethod:: mesher.unstructuredmesh.UnstMesh._generateVTKmesh
.. automethod:: mesher.unstructuredmesh.UnstMesh._localStrat
.. automethod:: mesher.unstructuredmesh.UnstMesh._meshAdvectorSphere
.. automethod:: mesher.unstructuredmesh.UnstMesh._meshUpliftSubsidence
.. automethod:: mesher.unstructuredmesh.UnstMesh._meshfrom_cell_list
.. automethod:: mesher.unstructuredmesh.UnstMesh._readErosionDeposition
.. automethod:: mesher.unstructuredmesh.UnstMesh._readStratLayers
.. automethod:: mesher.unstructuredmesh.UnstMesh._stratalRecord
.. automethod:: mesher.unstructuredmesh.UnstMesh._updateDispStrata
.. automethod:: mesher.unstructuredmesh.UnstMesh._updateRain
.. automethod:: mesher.unstructuredmesh.UnstMesh._updateTectonics
.. automethod:: mesher.unstructuredmesh.UnstMesh._xyz2lonlat
