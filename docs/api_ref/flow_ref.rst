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
      ~FAMesh._distributeDownstream

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
.. automethod:: flow.flowplex.FAMesh._distributeDownstream
