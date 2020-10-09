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
      ~FAMesh.riverIncision

   .. rubric:: Private Methods

   .. autosummary::

      ~FAMesh._distanceCoasts
      ~FAMesh._buildFlowDirection
      ~FAMesh._depressionlessSurface
      ~FAMesh._getErosionRate
      ~FAMesh._matrixFA
      ~FAMesh._matrix_build
      ~FAMesh._matrix_build_diag
      ~FAMesh._pitInformation
      ~FAMesh._solve_KSP

Public functions
---------------------

.. automethod:: flow.flowplex.FAMesh.flowAccumulation
.. automethod:: flow.flowplex.FAMesh.riverIncision


Private functions
---------------------

.. automethod:: flow.flowplex.FAMesh._distanceCoasts
.. automethod:: flow.flowplex.FAMesh._buildFlowDirection
.. automethod:: flow.flowplex.FAMesh._depressionlessSurface
.. automethod:: flow.flowplex.FAMesh._getErosionRate
.. automethod:: flow.flowplex.FAMesh._matrixFA
.. automethod:: flow.flowplex.FAMesh._matrix_build
.. automethod:: flow.flowplex.FAMesh._matrix_build_diag
.. automethod:: flow.flowplex.FAMesh._pitInformation
.. automethod:: flow.flowplex.FAMesh._solve_KSP
