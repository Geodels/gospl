#!/usr/bin/env python3
"""Flag MPI/PETSc collectives reached under a rank-local condition.

This encodes the grep sweep in AGENTS.md > MPI contract ("THE #1 parallel
deadlock class"): a collective reached on some ranks but not others blocks
forever at np>1, while serial runs always pass. The lint walks every function
in ``gospl/`` and reports a collective call that sits

* inside an ``if`` / ``while`` / conditional expression whose test is
  rank-local (``MPIrank``, ``.any()``/``.all()``, ``len(...)``, ``.size``,
  ``.shape``, ``np.any/np.all/np.count_nonzero``), or
* after an early ``return`` / ``continue`` / ``break`` taken under such a test
  in the same function (reported as ``early-exit``).

A test is treated as GLOBAL (safe) when every name it reads was assigned from a
reduction (``allreduce``/``Allreduce``/``bcast``/``Bcast``/``reduce``, a
collective Vec reduction such as ``vec.sum()``), from ``self.<config>``
attributes, or from literals. ``MPIrank`` comparisons are always rank-local.

It is a heuristic, not a proof: it cannot see through helper functions, and it
does not know which per-pit arrays are global. A finding is therefore a
question, not a verdict. When a flagged site is correct, say why in place:

    if self.flatModel and pit_ids.any():  # mpi-lint: ok pit_ids is global (pitParams)

The marker may sit on the guard line or on the collective's line. The reason
is mandatory; a bare ``# mpi-lint: ok`` is itself reported.

Usage::

    python scripts/lint_mpi_collectives.py            # lint gospl/, exit 1 on findings
    python scripts/lint_mpi_collectives.py --json     # machine-readable
    python scripts/lint_mpi_collectives.py path.py    # specific files
"""

from __future__ import annotations

import argparse
import ast
import json
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

# Method names that are collective whatever object they are called on.
ALWAYS_COLLECTIVE = {
    "localToGlobal", "globalToLocal", "localToLocal",
    "Allreduce", "allreduce", "Bcast", "bcast", "Allgatherv", "allgather",
    "Allgather", "Gatherv", "gather", "Scatterv", "scatter", "Reduce",
    "Barrier", "barrier", "Alltoall", "alltoall", "Alltoallv",
    "garbage_cleanup", "safe_garbage_cleanup",
    "assemblyBegin", "assemblyEnd", "zeroRows", "zeroRowsLocal",
    "multTranspose", "setUp", "createNest", "distribute",
}
# Collective only on PETSc objects: we require a non-numpy receiver.
# (copy/axpy/scale are "logically collective" but do not communicate, so they
# cannot deadlock; they are deliberately left out to keep the signal clean.)
PETSC_COLLECTIVE = {"solve", "mult", "norm", "dot", "duplicate"}
# Collective Vec reductions -- flagged only when the receiver is a known Vec.
VEC_REDUCTIONS = {"sum", "max", "min", "norm", "dot"}
# Calls whose RESULT is identical on every rank (makes a name global).
REDUCERS = {"allreduce", "Allreduce", "bcast", "Bcast", "reduce", "allgather",
            "Allgather"} | VEC_REDUCTIONS
NUMPY_ROOTS = {"np", "numpy", "linalg", "la", "scipy", "sp", "spla", "math"}
RANK_NAMES = {"MPIrank", "rank", "myrank", "comm_rank"}
RANK_LOCAL_CALLS = {"any", "all", "count_nonzero", "nonzero", "argwhere", "where"}
RANK_LOCAL_ATTRS = {"size", "shape"}
OK_RE = re.compile(r"#\s*mpi-lint:\s*ok\b(.*)")


@dataclass
class Finding:
    file: str
    line: int
    function: str
    collective: str
    guard_line: int
    guard: str
    kind: str          # "guard" | "early-exit" | "bare-ok"

    def text(self) -> str:
        where = f"{self.file}:{self.line}"
        if self.kind == "bare-ok":
            return f"{where}: `# mpi-lint: ok` without a reason (in {self.function})"
        verb = "inside" if self.kind == "guard" else "after an early exit under"
        return (f"{where}: collective `{self.collective}` {verb} rank-local "
                f"condition `{self.guard}` (line {self.guard_line}, in {self.function})")


def _root_name(node: ast.AST) -> str | None:
    while isinstance(node, (ast.Attribute, ast.Subscript, ast.Call)):
        node = node.value if not isinstance(node, ast.Call) else node.func
    return node.id if isinstance(node, ast.Name) else None


def _dotted(node: ast.AST) -> str:
    try:
        return ast.unparse(node)
    except Exception:  # pragma: no cover
        return "?"


def collect_vec_names(files: list[Path]) -> set[str]:
    """Attribute names assigned a PETSc Vec anywhere in the package."""
    names: set[str] = set()
    vec_makers = ("createGlobalVec", "createLocalVec", "createGlobalVector",
                  "createLocalVector", "duplicate", "createWithArray",
                  "createMPI", "createSeq", "getVecs", "createVecs")
    for f in files:
        try:
            tree = ast.parse(f.read_text())
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign) and isinstance(node.value, ast.Call):
                fn = node.value.func
                if isinstance(fn, ast.Attribute) and fn.attr in vec_makers:
                    for t in node.targets:
                        if isinstance(t, ast.Attribute):
                            names.add(t.attr)
                        elif isinstance(t, ast.Name):
                            names.add(t.id)
    return names


class FunctionLinter:
    def __init__(self, fn: ast.AST, path: str, lines: list[str], vec_names: set[str]):
        self.fn = fn
        self.path = path
        self.lines = lines
        self.vec_names = vec_names
        self.global_names: set[str] = set()
        self.local_vecs: set[str] = {
            t.id for n in ast.walk(fn) if isinstance(n, ast.Assign)
            and isinstance(n.value, ast.Call) and isinstance(n.value.func, ast.Attribute)
            and n.value.func.attr in ("createGlobalVec", "createLocalVec",
                                      "createGlobalVector", "createLocalVector",
                                      "duplicate", "createWithArray")
            for t in n.targets if isinstance(t, ast.Name)}
        self.findings: list[Finding] = []
        self.qualname = getattr(fn, "name", "<module>")

    # -- classification -------------------------------------------------
    def collective_name(self, call: ast.Call) -> str | None:
        fn = call.func
        if isinstance(fn, ast.Name):
            return fn.id if fn.id in {"safe_garbage_cleanup"} else None
        if not isinstance(fn, ast.Attribute):
            return None
        attr = fn.attr
        root = _root_name(fn.value)
        if attr in ALWAYS_COLLECTIVE:
            return attr
        if root in NUMPY_ROOTS:
            return None
        recv = fn.value
        recv_name = recv.attr if isinstance(recv, ast.Attribute) else (
            recv.id if isinstance(recv, ast.Name) else None)
        # A bare local name only counts as a Vec if THIS function made it one
        # (a local numpy array often shares a name with a Vec attribute).
        is_vec = (recv_name in self.vec_names and isinstance(recv, ast.Attribute)) \
            or (isinstance(recv, ast.Name) and recv_name in self.local_vecs)
        if attr in VEC_REDUCTIONS and is_vec:
            return f"{recv_name}.{attr}"
        # numpy arrays have no .norm/.mult/.solve/.duplicate, so on any
        # non-numpy receiver these are PETSc (the 2026-06-24 hang was
        # `vector1.norm()` on a function argument).
        if attr in PETSC_COLLECTIVE and attr != "dot":
            return f"{_dotted(recv)}.{attr}"
        if attr == "dot" and is_vec:
            return f"{recv_name}.dot"
        return None

    def is_rank_local(self, test: ast.AST) -> bool:
        """True if any part of the test is rank-local and not shown global."""
        for node in ast.walk(test):
            if isinstance(node, ast.Name) and node.id in RANK_NAMES:
                return True
            if isinstance(node, ast.Call):
                fn = node.func
                name = fn.attr if isinstance(fn, ast.Attribute) else (
                    fn.id if isinstance(fn, ast.Name) else None)
                if name in RANK_LOCAL_CALLS or name == "len":
                    target = fn.value if isinstance(fn, ast.Attribute) else (
                        node.args[0] if node.args else None)
                    if target is not None and not self._is_global_expr(target):
                        return True
            if isinstance(node, ast.Attribute) and node.attr in RANK_LOCAL_ATTRS:
                if not self._is_global_expr(node.value):
                    return True
        return False

    def _is_global_expr(self, node: ast.AST) -> bool:
        names = {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}
        names.discard("self")
        if not names:
            return False  # e.g. self.foo.any(): an attribute, assume rank-local
        return names <= self.global_names

    def _mark_globals(self, node: ast.AST) -> None:
        # MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, buf, ...) makes `buf` global.
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            c = node.value
            if isinstance(c.func, ast.Attribute) and c.func.attr in {"Allreduce", "Bcast"}:
                for a in c.args[:2]:
                    if isinstance(a, ast.Name):
                        self.global_names.add(a.id)
            return
        if isinstance(node, (ast.Assign, ast.AnnAssign)) and node.value is not None:
            value = node.value
            calls = [c for c in ast.walk(value) if isinstance(c, ast.Call)]
            is_reduce = any(
                isinstance(c.func, ast.Attribute) and (
                    c.func.attr in {"allreduce", "Allreduce", "bcast", "Bcast",
                                    "allgather", "Allgather"} or
                    (c.func.attr in VEC_REDUCTIONS and self.collective_name(c)))
                for c in calls)
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for t in targets:
                for n in ast.walk(t):
                    if isinstance(n, ast.Name):
                        if is_reduce:
                            self.global_names.add(n.id)
                        else:
                            self.global_names.discard(n.id)

    # -- traversal -------------------------------------------------------
    def ok_reason(self, lineno: int) -> str | None:
        m = OK_RE.search(self.lines[lineno - 1]) if 0 < lineno <= len(self.lines) else None
        return None if m is None else m.group(1).strip()

    def report(self, call: ast.Call, name: str, guard: ast.AST, kind: str):
        for ln in (call.lineno, guard.lineno):
            reason = self.ok_reason(ln)
            if reason is not None:
                if not reason:
                    self.findings.append(Finding(self.path, ln, self.qualname,
                                                 name, guard.lineno, "", "bare-ok"))
                return
        g = getattr(guard, "test", guard)
        self.findings.append(Finding(self.path, call.lineno, self.qualname, name,
                                     guard.lineno, _dotted(g)[:80], kind))

    def walk(self, stmts, guards: list[ast.AST], exited: list[ast.AST]):
        for st in stmts:
            self._mark_globals(st)
            if isinstance(st, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                continue  # linted separately
            # collectives in this statement's own expressions (not nested bodies)
            for call in self._own_calls(st):
                name = self.collective_name(call)
                if name is None:
                    continue
                if guards:
                    self.report(call, name, guards[-1], "guard")
                elif exited:
                    self.report(call, name, exited[-1], "early-exit")
            if isinstance(st, (ast.If, ast.While)):
                local = self.is_rank_local(st.test) and not self._balanced(st)
                inner = guards + [st] if local else guards
                self.walk(st.body, inner, exited)
                self.walk(st.orelse, inner, exited)
                if local and not guards and self._exits(st.body):
                    exited = exited + [st]
            elif isinstance(st, (ast.For, ast.AsyncFor)):
                self.walk(st.body, guards, exited)
                self.walk(st.orelse, guards, exited)
            elif isinstance(st, (ast.With, ast.AsyncWith)):
                self.walk(st.body, guards, exited)
            elif isinstance(st, ast.Try):
                self.walk(st.body, guards, exited)
                for h in st.handlers:
                    self.walk(h.body, guards, exited)
                self.walk(st.orelse, guards, exited)
                self.walk(st.finalbody, guards, exited)

    def _balanced(self, st: ast.AST) -> bool:
        """`if MPIrank == 0: bcast(x) else: bcast(None)`: both branches make
        the same collective calls, so every rank takes part (AGENTS.md: the
        bcast send/receive pattern is fine)."""
        if not isinstance(st, ast.If) or not st.orelse:
            return False

        def names(stmts):
            out = []
            for s_ in stmts:
                for c in ast.walk(s_):
                    if isinstance(c, ast.Call):
                        n = self.collective_name(c)
                        if n and n.split(".")[-1] in {"bcast", "Bcast", "Barrier",
                                                       "barrier", "gather", "Gatherv",
                                                       "scatter", "Scatterv"}:
                            out.append(n)
            return sorted(out)

        b = names(st.body)
        return bool(b) and b == names(st.orelse)

    @staticmethod
    def _exits(body) -> bool:
        return any(isinstance(s, (ast.Return, ast.Continue, ast.Break)) for s in body)

    def _own_calls(self, st: ast.stmt):
        """Calls in a statement, excluding nested statement bodies."""
        if isinstance(st, (ast.If, ast.While)):
            roots = [st.test]
        elif isinstance(st, (ast.For, ast.AsyncFor)):
            roots = [st.iter]
        elif isinstance(st, (ast.With, ast.AsyncWith)):
            roots = [i.context_expr for i in st.items]
        elif isinstance(st, ast.Try):
            roots = []
        else:
            roots = [st]
        for r in roots:
            for node in ast.walk(r):
                if isinstance(node, ast.Call):
                    yield node
                elif isinstance(node, ast.IfExp) and self.is_rank_local(node.test):
                    for c in ast.walk(node.body):
                        if isinstance(c, ast.Call) and self.collective_name(c):
                            self.report(c, self.collective_name(c), node, "guard")
                    for c in ast.walk(node.orelse):
                        if isinstance(c, ast.Call) and self.collective_name(c):
                            self.report(c, self.collective_name(c), node, "guard")

    def run(self) -> list[Finding]:
        self.walk(self.fn.body, [], [])
        # de-duplicate (IfExp branches can be visited twice)
        seen, out = set(), []
        for f in self.findings:
            key = (f.line, f.collective, f.kind)
            if key not in seen:
                seen.add(key)
                out.append(f)
        return out


def lint_file(path: Path, vec_names: set[str]) -> list[Finding]:
    src = path.read_text()
    tree = ast.parse(src)
    lines = src.splitlines()
    rel = str(path.relative_to(REPO)) if path.is_relative_to(REPO) else str(path)
    out: list[Finding] = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            out += FunctionLinter(node, rel, lines, vec_names).run()
    return sorted(out, key=lambda f: f.line)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("paths", nargs="*", help="files to lint (default: gospl/**/*.py)")
    ap.add_argument("--json", action="store_true", help="emit JSON findings")
    args = ap.parse_args(argv)

    pkg = sorted((REPO / "gospl").rglob("*.py"))
    files = [Path(p).resolve() for p in args.paths] if args.paths else pkg
    vec_names = collect_vec_names(pkg)
    findings = [f for p in files for f in lint_file(p, vec_names)]

    if args.json:
        print(json.dumps([asdict(f) for f in findings], indent=1))
    else:
        for f in findings:
            print(f.text())
        n = len(findings)
        print(f"\n{n} finding(s) in {len(files)} file(s)." if n else
              f"OK: no rank-local collectives in {len(files)} file(s).",
              file=sys.stderr)
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
