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
      ~FAMesh._distanceCoasts
      ~FAMesh._distributeDownstream
      ~FAMesh._getErosionRate
      ~FAMesh._globalCoastsTree
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
.. automethod:: flow.flowplex.FAMesh._distanceCoasts
.. automethod:: flow.flowplex.FAMesh._distributeDownstream
.. automethod:: flow.flowplex.FAMesh._getErosionRate
.. automethod:: flow.flowplex.FAMesh._globalCoastsTree
.. automethod:: flow.flowplex.FAMesh._matrix_build
.. automethod:: flow.flowplex.FAMesh._matrix_build_diag
.. automethod:: flow.flowplex.FAMesh._solve_KSP
