"""
goSPL — parallel landscape-evolution model.

Single PETSc initialization site. Until 2026-06 every submodule called
``petsc4py.init(sys.argv)`` at module import time (15 files); that pattern
was harmless (the call is idempotent) but obscured where PETSc state was
actually being created. The audit (REFACTOR_AUDIT.md §1.2) flagged this
as worth consolidating; AGENTS.md > MPI contract notes the same.

Because Python runs ``gospl/__init__.py`` before any submodule of the
``gospl`` package, this single call guarantees PETSc is initialised
before any module-level code in submodules executes (e.g.
``MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()`` in flowplex.py, which
runs at import time and requires PETSc to be live).

Do NOT re-introduce ``petsc4py.init`` in submodules. If you need PETSc
in a new file, just ``import petsc4py`` and use ``petsc4py.PETSc.X``
directly — init is guaranteed by the time your module loads.
"""

import sys

import petsc4py

petsc4py.init(sys.argv)

# --- version ---
try:
    from importlib.metadata import version, PackageNotFoundError
    __version__ = version("gospl")
except PackageNotFoundError:
    __version__ = "unknown"
