.. _plate_ref:


=================
Class Tectonics
=================

.. autoclass:: mesher.tectonics.Tectonics

   .. rubric:: Initialise

   .. autosummary::

      ~Tectonics.__init__

   .. rubric:: Public Methods

   .. autosummary::

      ~Tectonics.getTectonics
      ~Tectonics.updatePaleoZ

   .. rubric:: Private Methods

   .. autosummary::

      ~Tectonics._buildAdvecMat
      ~Tectonics._advectorIIOE2
      ~Tectonics._varAdvector
      ~Tectonics._readAdvectionData
      ~Tectonics._advectPlates
      ~Tectonics._advectStrati
      ~Tectonics._advectStratiScheme
      ~Tectonics._fetchOwnerRows
      ~Tectonics._ownerPetscIndex
      ~Tectonics._cylinderDisplace
      ~Tectonics._cylinderVelocity

Public functions
---------------------

.. automethod:: mesher.tectonics.Tectonics.getTectonics
.. automethod:: mesher.tectonics.Tectonics.updatePaleoZ


Private functions
---------------------

.. automethod:: mesher.tectonics.Tectonics._buildAdvecMat
.. automethod:: mesher.tectonics.Tectonics._advectorIIOE2
.. automethod:: mesher.tectonics.Tectonics._varAdvector
.. automethod:: mesher.tectonics.Tectonics._readAdvectionData
.. automethod:: mesher.tectonics.Tectonics._advectPlates
.. automethod:: mesher.tectonics.Tectonics._advectStrati
.. automethod:: mesher.tectonics.Tectonics._advectStratiScheme
.. automethod:: mesher.tectonics.Tectonics._fetchOwnerRows
.. automethod:: mesher.tectonics.Tectonics._ownerPetscIndex
.. automethod:: mesher.tectonics.Tectonics._cylinderDisplace
.. automethod:: mesher.tectonics.Tectonics._cylinderVelocity
