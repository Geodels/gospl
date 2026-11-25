.. _solvers:

==================================
Solvers: defaults and alternatives
==================================

Each process of goSPL is solved with a PETSc solver chosen for robustness and
speed on the meshes goSPL is used with. This page lists, for every process,
the default solver, the alternatives that have been **tested**, and how to
select them for a run.

.. contents:: On this page
   :local:
   :depth: 1


How to change a solver
----------------------

Two routes, depending on the process:

- **Input keys** for the solvers that have a choice in the YAML input (the
  soil, the non-linear stream power law, the marine diffusion): set them in
  the block of the process, for example ``diffusion: marineSolver: ts``.
- **PETSc options** for everything else. Every solver has its own options
  *prefix* (listed below); any PETSc option given with that prefix replaces
  the default for that solver only. Pass them in the ``PETSC_OPTIONS``
  environment variable (it reaches every MPI rank):

  .. code:: bash

      export PETSC_OPTIONS="-gw_solute_ksp_type preonly -gw_solute_pc_type lu \
          -gw_solute_pc_factor_mat_solver_type mumps"
      mpirun -np 8 gospl -i input.yml

  To check what a solver actually uses, add ``-<prefix>ksp_view`` (or
  ``-<prefix>snes_view``, ``-<prefix>ts_view``); ``-<prefix>ksp_converged_reason``
  prints the outcome of every solve. Some solvers (MUMPS, hypre) are external
  PETSc packages: goSPL uses them when PETSc was built with them (the conda,
  pip and Docker builds are) and falls back as noted below otherwise.

Unless stated otherwise, the timings below come from the 200 x 100 km
margin of the goSPL-examples ``mixed_carbonate_clastic`` example (47,672 and
189,497 nodes) on a laptop (Apple M4 Max), November 2026.


Flow accumulation and sediment routing
--------------------------------------

The implicit flow accumulation ``(I - W^T) q = b`` (water, sediment, solutes).

.. list-table::
   :header-rows: 1
   :widths: 22 30 48

   * - Solver
     - Default
     - Alternatives and settings
   * - main routing (``flowacc_``)
     - ``fgmres`` + ``bjacobi`` (ILU)
     - ``GOSPL_FLOW_KSP=richardson`` and ``GOSPL_FLOW_PC`` (environment
       variables) select the Krylov method and preconditioner; ``fgmres`` is
       the default because Richardson diverged partition-dependently on long
       drainage chains.
   * - exact rescue (``flowaccx_``)
     - block-Jacobi with exact LU blocks, used when the main solve fails
       (sticky)
     - ``GOSPL_FLOW_IDA_EXACT=1`` uses it from the start.
   * - bounded fallback (``flowaccfb_``)
     - ``richardson``, no preconditioner
     - none (last resort).

See :ref:`flow` for the pond and repair rules that surround these solves.


River incision (stream power law)
---------------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 30 48

   * - Process
     - Default
     - Alternatives and settings
   * - linear SPL, no deposition
     - direct accumulation solve (flow solver above)
     -
   * - linear SPL with deposition (``spl: G > 0``, ``spl_ed_``)
     - ``tfqmr`` with a Schur fieldsplit (``gasm`` blocks)
     - the fieldsplit and its blocks through ``spl_ed_`` options
       (``-spl_ed_pc_fieldsplit_type additive``,
       ``-spl_ed_fieldsplit_h_pc_type hypre``).
   * - non-linear SPL (``spl: n != 1``), no deposition
     - SNES ``nrichardson`` with the analytic Jacobian
     -
   * - non-linear SPL with deposition (``nlspled_``)
     - SNES ``qn`` (L-BFGS), ``ngmres`` + ``pcType`` fallback (``nlspledfb_``)
     - ``spl: solver: ngmres`` makes the multigrid-preconditioned solver the
       primary; ``spl: pcType``, ``rtol``, ``atol``, ``maxIter``. Not yet
       measured on the regional margin (the soil SNES, which had the same
       default, now uses ``ngmres``: see below).
   * - soil-aware SPL (``soil:``, ``soilspl_``)
     - SNES ``ngmres`` + ``nrichardson`` / CG / hypre BoomerAMG, rtol 1e-6;
       ``qn`` fallback (``soilsplfb_``, rtol x 100); adaptive sub-stepping if
       both fail
     - ``soil: solver: qn`` (the default until November 2026: it failed at
       nearly every step on the regional margin and on the global
       ``continental_geochem`` example, and the fallback finished at 1e-4;
       ``ngmres`` converged at 1e-6 with no fallback, in 58 % of the time at
       47k nodes). ``soil: pcType: gamg | bjacobi | asm``; ``rtol: 1.e-5``
       with ``hypre`` is a tested fast setting (user guide).


Hillslope and soil diffusion
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 30 48

   * - Process
     - Default
     - Alternatives and settings
   * - linear hillslope (``hilllin_``)
     - KSP Richardson + ``bjacobi``, operator cached and rebuilt when the
       coastline moves
     - ``hilllin_`` options.
   * - non-linear hillslope (``hillnl_``)
     - SNES ``qn``
     - bare ``ngmres`` aborted on fine meshes; ``hillnl_`` options.
   * - soil diffusion
     - TS ``rosw`` (Rosenbrock), each stage solved with GMRES + ``gasm``
       (atol 1e-3 m)
     - do not set the stage solve to ``preonly``: a Rosenbrock scheme needs
       each stage solved.


Marine and lake deposition
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 30 48

   * - Process
     - Default
     - Alternatives and settings
   * - marine and lake diffusion (``marpicard_``)
     - ``marineSolver: picard``: lagged-diffusivity backward Euler,
       ``picardSub: 20`` sub-steps, ``picardIts: 2`` Picard iterations, each
       a linear Richardson + ``bjacobi`` solve
     - ``marineSolver: ts``: adaptive non-linear TS ``rosw`` (the default until
       November 2026). Tested on 189,497 nodes: ``picard`` with 20 sub-steps
       took 68 % of the ``ts`` step and its deposits were within 1.4 % rms of
       ``ts``; with 10 sub-steps 51 % and 9 %; 30 sub-steps and 3 iterations
       114 % and 1.4 %. The error of ``picard`` grows as the cells shrink: on
       finer meshes scale ``picardSub`` by about ``(dx_old/dx_new)^2``. If
       ``ts`` gives up during a run, goSPL switches to ``picard``. ``picard``
       also limits each cell's export to its deposit, which ``ts`` does not.
   * - marine flow-direction smoothing (``marsmooth_``)
     - KSP Richardson + ``bjacobi``, cached
     - ``marsmooth_`` options.
   * - lithology and provenance columns in the marine diffusion
       (``provdiff_``)
     - the operator of the bulk diffusion, ``picardSub`` backward-Euler
       sub-steps, Richardson + ``bjacobi``
     - follows ``picardSub``.
   * - coupled marine deposition with ``Gmar > 0`` on a sphere
       (``marine_dep_``)
     - fieldsplit (rebuilt each call)
     - ``marine_dep_`` options.


Flexure and orography
---------------------

.. list-table::
   :header-rows: 1
   :widths: 22 30 48

   * - Process
     - Default
     - Alternatives and settings
   * - flat flexure (``flexure: method: fem``, ``flexfem_``)
     - direct LU (serial) / MUMPS (parallel), factorised once
     - an iterative solver for meshes too large to factorise:
       ``-flexfem_ksp_type gmres -flexfem_pc_type gamg``.
   * - global flexure (``method: global``)
     - spherical-harmonic solution on a regular grid (rank 0)
     -
   * - orographic rain (``oro_``)
     - KSP ``gmres`` + ``bjacobi``, rtol 1e-8, cached operators
     - ``oro_`` options.


Groundwater and geochemistry
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 30 48

   * - Process
     - Default
     - Alternatives and settings
   * - water-table head (``gw_``)
     - ``fgmres`` + hypre BoomerAMG
     - AMG is required (block-Jacobi stalls on the head operator).
   * - solute transport (``gw_solute_``)
     - ``fgmres`` + hypre BoomerAMG, rtol 1e-10, warm start; direct
       factorisation (``gw_solutelu_``, MUMPS) as automatic fallback of a
       failed solve, and as the default when PETSc has no hypre
     - the direct solve as primary:
       ``-gw_solute_ksp_type preonly -gw_solute_pc_type lu
       -gw_solute_pc_factor_mat_solver_type mumps`` (the default until
       November 2026). Tested on 189,497 nodes: the iterative solve took
       2.06 / 0.83 / 0.67 s per step at 1 / 4 / 8 ranks against 2.99 / 2.24 /
       2.41 s for MUMPS, with less memory with several ranks and the same
       mass budget (closure ~1e-6). In serial, ask MUMPS rather than PETSc's
       own LU (ten times slower on this operator).
