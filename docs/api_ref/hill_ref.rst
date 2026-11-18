.. _hill_ref:


==============
Class hillSLP
==============

.. autoclass:: sed.hillslope.hillSLP

   .. rubric:: Initialise

   .. autosummary::

      ~hillSLP.__init__

   .. rubric:: Public Methods

   .. autosummary::

      ~hillSLP.getHillslope
      ~hillSLP.hillslopeBeforeRouting

   .. rubric:: Private Methods

   .. autosummary::

      ~hillSLP._hillSlope
      ~hillSLP._buildDiffMat
      ~hillSLP._channelMask
      ~hillSLP._channelDiffusion
      ~hillSLP._channelComposition
      ~hillSLP._checkChannelOptions
      ~hillSLP._assembleDiffMat
      ~hillSLP._assembleDiffMatCSR
      ~hillSLP._donorDiffCoeffs
      ~hillSLP._makeDiffusionKSP
      ~hillSLP._solveSmooth
      ~hillSLP._diff_nl_monitor
      ~hillSLP._form_residual_nl_hillslope
      ~hillSLP._hillSlopeNL
      ~hillSLP._hillslopeFineFraction
      ~hillSLP._evalFunctionMardDiff
      ~hillSLP._evalJacobianMardDiff
      ~hillSLP._evalSolutionMardDiff
      ~hillSLP._diffuseOcean
      ~hillSLP._diffuseImplicit
      ~hillSLP._diffuseImplicitPicard

Public functions
---------------------

.. automethod:: sed.hillslope.hillSLP.getHillslope
.. automethod:: sed.hillslope.hillSLP.hillslopeBeforeRouting


Private functions
---------------------

.. automethod:: sed.hillslope.hillSLP._hillSlope
.. automethod:: sed.hillslope.hillSLP._buildDiffMat
.. automethod:: sed.hillslope.hillSLP._channelMask
.. automethod:: sed.hillslope.hillSLP._channelDiffusion
.. automethod:: sed.hillslope.hillSLP._channelComposition
.. automethod:: sed.hillslope.hillSLP._checkChannelOptions
.. automethod:: sed.hillslope.hillSLP._assembleDiffMat
.. automethod:: sed.hillslope.hillSLP._assembleDiffMatCSR
.. automethod:: sed.hillslope.hillSLP._donorDiffCoeffs
.. automethod:: sed.hillslope.hillSLP._makeDiffusionKSP
.. automethod:: sed.hillslope.hillSLP._solveSmooth
.. automethod:: sed.hillslope.hillSLP._diff_nl_monitor
.. automethod:: sed.hillslope.hillSLP._form_residual_nl_hillslope
.. automethod:: sed.hillslope.hillSLP._hillSlopeNL
.. automethod:: sed.hillslope.hillSLP._hillslopeFineFraction
.. automethod:: sed.hillslope.hillSLP._evalFunctionMardDiff
.. automethod:: sed.hillslope.hillSLP._evalJacobianMardDiff
.. automethod:: sed.hillslope.hillSLP._evalSolutionMardDiff
.. automethod:: sed.hillslope.hillSLP._diffuseOcean
.. automethod:: sed.hillslope.hillSLP._diffuseImplicit
.. automethod:: sed.hillslope.hillSLP._diffuseImplicitPicard
