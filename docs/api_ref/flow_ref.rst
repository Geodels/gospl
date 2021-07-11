.. _flow_ref:


==============
Class FAMesh
==============

.. autoclass:: flow.flowplex.FAMesh

   .. rubric:: Initialise

   .. autosummary::

      ~FAMesh.__init__

   .. rubric:: Public Methods

   .. autosummary::

      ~FAMesh.flowAccumulation
      ~FAMesh.matrixFlow
      ~FAMesh.riverIncision

   .. rubric:: Private Methods

   .. autosummary::

      ~FAMesh._buildFlowDirection
      ~FAMesh._depressionlessSurface
      ~FAMesh._distanceCoasts
      ~FAMesh._distributeFlowExcess
      ~FAMesh._flowSurface
      ~FAMesh._getErosionRate
      ~FAMesh._matrix_build
      ~FAMesh._matrix_build_diag
      ~FAMesh._solve_KSP

Public functions
---------------------

.. automethod:: flow.flowplex.FAMesh.flowAccumulation
.. automethod:: flow.flowplex.FAMesh.matrixFlow
.. automethod:: flow.flowplex.FAMesh.riverIncision


Private functions
---------------------

.. automethod:: flow.flowplex.FAMesh._buildFlowDirection
.. automethod:: flow.flowplex.FAMesh._depressionlessSurface
.. automethod:: flow.flowplex.FAMesh._distanceCoasts
.. automethod:: flow.flowplex.FAMesh._distributeFlowExcess
.. automethod:: flow.flowplex.FAMesh._flowSurface
.. automethod:: flow.flowplex.FAMesh._getErosionRate
.. automethod:: flow.flowplex.FAMesh._matrix_build
.. automethod:: flow.flowplex.FAMesh._matrix_build_diag
.. automethod:: flow.flowplex.FAMesh._solve_KSP
