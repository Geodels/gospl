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

      ~FAMesh.matrixFlow
      ~FAMesh.flowAccumulation

   .. rubric:: Private Methods

   .. autosummary::

      ~FAMesh._matrix_build
      ~FAMesh._matrix_build_diag
      ~FAMesh._make_reasons
      ~FAMesh._solve_KSP2
      ~FAMesh._solve_KSP
      ~FAMesh._buildFlowDirection
      ~FAMesh._potentialLakeEvap
      ~FAMesh._lakeEvapCurve
      ~FAMesh._equilibriumLakeLevel
      ~FAMesh._distributeDownstream
      ~FAMesh._cascadeResetProgress
      ~FAMesh._cascadeStopReason
      ~FAMesh._losingStreamSolve
      ~FAMesh._undrainedCap
      ~FAMesh._pondRoutingRegion
      ~FAMesh._reportRoutingRegion
      ~FAMesh._takeRoutingPond
      ~FAMesh._trueResidualOK
      ~FAMesh._warmStartGuard
      ~FAMesh._solveIDAExact
      ~FAMesh._routingSurface
      ~FAMesh._runoffWork

Public functions
---------------------

.. automethod:: flow.flowplex.FAMesh.matrixFlow
.. automethod:: flow.flowplex.FAMesh.flowAccumulation


Private functions
---------------------

.. automethod:: flow.flowplex.FAMesh._matrix_build
.. automethod:: flow.flowplex.FAMesh._matrix_build_diag
.. automethod:: flow.flowplex.FAMesh._make_reasons
.. automethod:: flow.flowplex.FAMesh._solve_KSP2
.. automethod:: flow.flowplex.FAMesh._solve_KSP
.. automethod:: flow.flowplex.FAMesh._buildFlowDirection
.. automethod:: flow.flowplex.FAMesh._potentialLakeEvap
.. automethod:: flow.flowplex.FAMesh._lakeEvapCurve
.. automethod:: flow.flowplex.FAMesh._equilibriumLakeLevel
.. automethod:: flow.flowplex.FAMesh._distributeDownstream
.. automethod:: flow.flowplex.FAMesh._cascadeResetProgress
.. automethod:: flow.flowplex.FAMesh._cascadeStopReason
.. automethod:: flow.flowplex.FAMesh._losingStreamSolve
.. automethod:: flow.flowplex.FAMesh._undrainedCap
.. automethod:: flow.flowplex.FAMesh._pondRoutingRegion
.. automethod:: flow.flowplex.FAMesh._reportRoutingRegion
.. automethod:: flow.flowplex.FAMesh._takeRoutingPond
.. automethod:: flow.flowplex.FAMesh._trueResidualOK
.. automethod:: flow.flowplex.FAMesh._warmStartGuard
.. automethod:: flow.flowplex.FAMesh._solveIDAExact
.. automethod:: flow.flowplex.FAMesh._routingSurface
.. automethod:: flow.flowplex.FAMesh._runoffWork
