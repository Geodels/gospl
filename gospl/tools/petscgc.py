"""
Robust wrapper around PETSc's collective garbage collection.

``petsc4py.PETSc.garbage_cleanup()`` flushes PETSc's *delayed-destroy*
registry — objects created collectively but dropped (via Python refcount / GC)
in a rank-divergent order. It is a memory-hygiene optimisation, not a
correctness requirement, and on a healthy run it succeeds on every rank.

It is, however, a **collective** (it does an ``MPI_Allreduce`` to intersect the
per-rank garbage keys). If an *upstream* collective has already desynced the
ranks — e.g. a collective accidentally gated on a rank-local condition, or the
PETSc-3.21.x ``GarbageKeyAllReduceIntersect`` ``MPI_ERR_BUFFER`` edge case — the
cleanup fails **asymmetrically**: one rank raises while another blocks forever
inside the ``Allreduce``. Swallowing the error is therefore *unsafe* — it would
leave the blocked rank deadlocked (the run hangs and must be killed by hand).

So a failure here is treated as fatal: we print a clear diagnostic and
``MPI_Abort`` the whole job. Aborting from the rank that raised kills the
partner blocked in the collective, turning the hang into an immediate, clean
termination. In practice this should never fire — the known trigger (a
rank-0-gated ``Vec.norm``/``Mat.norm`` in the flow-KSP diagnostic) was fixed in
``flowplex._solve_KSP2`` — it is a defensive backstop for any future desync.
"""
import petsc4py


def safe_garbage_cleanup():
    """Run ``petsc4py.PETSc.garbage_cleanup()``; on failure, abort ALL ranks
    cleanly rather than risk a deadlock (see module docstring)."""
    try:
        petsc4py.PETSc.garbage_cleanup()
    except Exception as exc:
        from mpi4py import MPI

        if petsc4py.PETSc.COMM_WORLD.getRank() == 0:
            print(
                "[gospl] FATAL: PETSc garbage_cleanup() failed (%s). This is a "
                "collective call, so the failure means an upstream MPI desync "
                "(some collective taken by only a subset of ranks). Aborting "
                "all ranks to avoid a deadlock." % type(exc).__name__,
                flush=True,
            )
        MPI.COMM_WORLD.Abort(1)
