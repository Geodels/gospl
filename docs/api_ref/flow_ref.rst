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

      ~FAMesh.erodepSPL
      ~FAMesh.flowAccumulation
      ~FAMesh.matrixFlow

   .. rubric:: Private Methods

   .. autosummary::

      ~FAMesh._buildFlowDirection
      ~FAMesh._coupledEDSystem
      ~FAMesh._distributeDownstream
      ~FAMesh._eroMats
      ~FAMesh._getEroDepRate
      ~FAMesh._matrix_build
      ~FAMesh._matrix_build_diag
      ~FAMesh._solve_KSP

Public functions
---------------------

.. automethod:: flow.flowplex.FAMesh.erodepSPL
.. automethod:: flow.flowplex.FAMesh.flowAccumulation
.. automethod:: flow.flowplex.FAMesh.matrixFlow


Private functions
---------------------

.. automethod:: flow.flowplex.FAMesh._buildFlowDirection
.. automethod:: flow.flowplex.FAMesh._coupledEDSystem
.. automethod:: flow.flowplex.FAMesh._distributeDownstream
.. automethod:: flow.flowplex.FAMesh._eroMats
.. automethod:: flow.flowplex.FAMesh._getEroDepRate
.. automethod:: flow.flowplex.FAMesh._matrix_build
.. automethod:: flow.flowplex.FAMesh._matrix_build_diag
.. automethod:: flow.flowplex.FAMesh._solve_KSP
