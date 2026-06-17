.. _soilspl_ref:


==============
Class soilSPL
==============

.. autoclass:: eroder.soilSPL.soilSPL

   .. rubric:: Initialise

   .. autosummary::

      ~soilSPL.__init__

   .. rubric:: Public Methods

   .. autosummary::

      ~soilSPL.updateSoilThickness
      ~soilSPL.erodepSPLsoil
      ~soilSPL.diffuseSoil

   .. rubric:: Private Methods

   .. autosummary::

      ~soilSPL._form_residual_soil
      ~soilSPL._monitorsoil
      ~soilSPL._build_soil_snes
      ~soilSPL._solveSoil
      ~soilSPL._getEroDepRateSoil
      ~soilSPL._evalFunctionSoil
      ~soilSPL._evalJacobianSoil
      ~soilSPL._evalSolutionSoil

Public functions
---------------------

.. automethod:: eroder.soilSPL.soilSPL.updateSoilThickness
.. automethod:: eroder.soilSPL.soilSPL.erodepSPLsoil
.. automethod:: eroder.soilSPL.soilSPL.diffuseSoil


Private functions
---------------------

.. automethod:: eroder.soilSPL.soilSPL._form_residual_soil
.. automethod:: eroder.soilSPL.soilSPL._monitorsoil
.. automethod:: eroder.soilSPL.soilSPL._build_soil_snes
.. automethod:: eroder.soilSPL.soilSPL._solveSoil
.. automethod:: eroder.soilSPL.soilSPL._getEroDepRateSoil
.. automethod:: eroder.soilSPL.soilSPL._evalFunctionSoil
.. automethod:: eroder.soilSPL.soilSPL._evalJacobianSoil
.. automethod:: eroder.soilSPL.soilSPL._evalSolutionSoil
