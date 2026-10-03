"""
scripts/lint_mpi_collectives.py: the static guard for the #1 parallel deadlock
class (AGENTS.md > MPI contract). Each case is a minimal reproduction of a bug
shape that has actually deadlocked goSPL at np>1, or of a pattern the lint must
NOT flag (otherwise it gets ignored).
"""

from __future__ import annotations

import importlib.util
import sys
import textwrap
from pathlib import Path

import pytest

pytestmark = pytest.mark.tools

_LINT = Path(__file__).resolve().parents[1] / "scripts" / "lint_mpi_collectives.py"


@pytest.fixture(scope="module")
def lint():
    spec = importlib.util.spec_from_file_location("lint_mpi_collectives", _LINT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod        # @dataclass resolves its module here
    spec.loader.exec_module(mod)
    return mod


def _findings(lint, tmp_path, src):
    f = tmp_path / "case.py"
    f.write_text(textwrap.dedent(src))
    return lint.lint_file(f, vec_names={"FAG", "tmp", "hGlobal"})


def test_flags_collective_norm_under_rank0(lint, tmp_path):
    # The 2026-06-24 flowplex._solve_KSP2 hang: Vec.norm/Mat.norm under rank 0.
    out = _findings(lint, tmp_path, """
        def diag(self, vector1, matrix):
            if MPIrank == 0:
                rhs = vector1.norm()
                m = matrix.norm()
                print(rhs, m)
    """)
    assert {f.collective for f in out} == {"vector1.norm", "matrix.norm"}
    assert all(f.kind == "guard" for f in out)


def test_flags_scatter_under_local_any(lint, tmp_path):
    # The sedplex closed-sink deadlock: localToGlobal under a rank-local .any().
    out = _findings(lint, tmp_path, """
        def dep(self):
            if self._closedDepo.any():
                self.dm.localToGlobal(self.tmpL, self.tmp)
    """)
    assert [f.collective for f in out] == ["localToGlobal"]


def test_reduced_flag_is_safe(lint, tmp_path):
    # The fix recipe: reduce the rank-local flag first.
    out = _findings(lint, tmp_path, """
        def dep(self):
            has_closed = MPI.COMM_WORLD.allreduce(self._closedDepo.any(), op=MPI.LOR)
            if has_closed:
                self.dm.localToGlobal(self.tmpL, self.tmp)
    """)
    assert out == []


def test_balanced_bcast_is_safe(lint, tmp_path):
    out = _findings(lint, tmp_path, """
        def mesh(self):
            if MPIrank == 0:
                MPIcomm.bcast(cells.shape, root=0)
            else:
                shape = MPIcomm.bcast(None, root=0)
    """)
    assert out == []


def test_flags_collective_after_rank_local_early_return(lint, tmp_path):
    out = _findings(lint, tmp_path, """
        def tilt(self, sel_node):
            if not sel_node.any():
                return None
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, buf, op=MPI.SUM)
    """)
    assert [(f.collective, f.kind) for f in out] == [("Allreduce", "early-exit")]


def test_vec_sum_only_on_known_vecs(lint, tmp_path):
    out = _findings(lint, tmp_path, """
        def f(self, arr):
            if arr.size > 0:
                a = arr.sum()            # numpy: fine
                b = self.FAG.sum()       # PETSc Vec: collective
    """)
    assert [f.collective for f in out] == ["FAG.sum"]


def test_ok_marker_needs_a_reason(lint, tmp_path):
    out = _findings(lint, tmp_path, """
        def f(self):
            if pit_select.any():  # mpi-lint: ok pit_select is global
                self.dm.localToGlobal(self.tmpL, self.tmp)
            if eV.any():  # mpi-lint: ok
                self.dm.globalToLocal(self.tmp, self.tmpL)
    """)
    assert [f.kind for f in out] == ["bare-ok"]


def test_codebase_is_clean(lint):
    """The package itself must lint clean (CI runs the same check)."""
    assert lint.main([]) == 0
