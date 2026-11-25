.. _utils_ref:

==========================
Run-time utilities
==========================

Helpers used by every run: the ``gospl`` command, the wall-clock profiler,
the elevation-spike probe (``-v`` runs), the PETSc garbage-collection
wrapper and the constants shared by the kernels.

.. contents:: On this page
    :local:
    :depth: 1


Command line (``gospl``)
------------------------

.. automodule:: cli

.. autofunction:: cli.main
.. autofunction:: cli._parse_args


Wall-clock profiler
-------------------

.. automodule:: tools.profiler

.. autoclass:: tools.profiler.Profiler
   :members: phase, start, stop, reduce, report
   :private-members: _accumulate


Elevation-spike probe
---------------------

.. automodule:: tools.zprobe

.. autofunction:: tools.zprobe.probeZ
.. autofunction:: tools.zprobe.resetZ
.. autofunction:: tools.zprobe.reportVolume
.. autofunction:: tools.zprobe._inputMeshId


PETSc garbage collection
------------------------

.. automodule:: tools.petscgc

.. autofunction:: tools.petscgc.safe_garbage_cleanup


Constants
---------

Sentinels, thresholds and floors shared by the kernels
(``gospl/tools/constants.py``), each with its role:

.. literalinclude:: ../../gospl/tools/constants.py
   :language: python
   :start-at: # Pre-fill value before MPI
