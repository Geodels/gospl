#!/usr/bin/env python3
"""PostToolUse hook (Edit/Write/MultiEdit): fast checks on the file just edited.

* ``gospl/**/*.py``, ``tests/**/*.py``, ``scripts/**/*.py``, ``benchmarks/**/*.py``:
  byte-compile + pyflakes ERRORS only (undefined names, syntax errors -- the
  package already carries ~75 harmless unused-import warnings, which are not
  reported), and for ``gospl/`` the MPI-collective lint on that file
  (AGENTS.md > MPI contract). A finding blocks with exit 2 so the agent sees it.
* ``fortran/functions.F90`` / ``functions.pyf``: a reminder that both must change
  together and that ``scripts/fcheck.sh`` must run (AGENTS.md checklist 12).
* ``meson.build`` / ``conda/meta.yaml`` / ``docs/conf.py``: the four-place
  version-bump reminder (docs/dev/RELEASE.md).

Reads the hook JSON on stdin; never modifies anything. Each check is a few
hundred ms at most.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ERROR_MARKERS = ("undefined name", "referenced before assignment", "invalid syntax",
                 "SyntaxError", "unexpected indent", "unindent does not match",
                 "is not defined", "EOL while scanning", "unterminated")


def context(msg: str) -> None:
    print(json.dumps({"hookSpecificOutput": {"hookEventName": "PostToolUse",
                                             "additionalContext": msg}}))
    sys.exit(0)


def main() -> None:
    try:
        data = json.load(sys.stdin)
    except Exception:
        sys.exit(0)
    fp = (data.get("tool_input") or {}).get("file_path") or ""
    if not fp:
        sys.exit(0)
    path = Path(fp).resolve()
    try:
        rel = path.relative_to(REPO)
    except ValueError:
        sys.exit(0)
    top = rel.parts[0] if rel.parts else ""
    name = rel.name

    if name in ("functions.F90", "functions.pyf") and top == "fortran":
        other = "functions.pyf" if name.endswith(".F90") else "functions.F90"
        context(f"Fortran edit: change fortran/{other} too if a signature changed "
                "(f2py drift is silent), bound-check every write into a (npoints, 12) "
                "table, then run scripts/fcheck.sh before committing "
                "(fortran/AGENTS.md).")

    if str(rel) in ("meson.build", "conda/meta.yaml", "docs/conf.py",
                    "docs/_static/version_switch.json"):
        context("If this is a version bump: meson.build:4, conda/meta.yaml:2, "
                "docs/conf.py and docs/_static/version_switch.json move together, "
                "no leading zeros (docs/dev/RELEASE.md).")

    if path.suffix != ".py" or top not in ("gospl", "tests", "scripts", "benchmarks"):
        sys.exit(0)

    problems = []
    py = sys.executable
    r = subprocess.run([py, "-m", "pyflakes", str(path)], capture_output=True, text=True)
    if r.returncode not in (0, 1) and "No module named pyflakes" in r.stderr:
        r = subprocess.run([py, "-m", "py_compile", str(path)], capture_output=True,
                           text=True)
        if r.returncode:
            problems.append(r.stderr.strip())
    else:
        for line in (r.stdout + r.stderr).splitlines():
            if any(m in line for m in ERROR_MARKERS):
                problems.append(line)

    if top == "gospl" and not problems:
        lint = REPO / "scripts" / "lint_mpi_collectives.py"
        r = subprocess.run([py, str(lint), str(path)], capture_output=True, text=True,
                           cwd=REPO)
        if r.returncode:
            problems += [l for l in r.stdout.splitlines() if l.strip()]
            problems.append("(AGENTS.md > MPI contract: reduce a rank-local flag "
                            "with allreduce(..., op=MPI.LOR) before it gates a "
                            "collective, or justify with '# mpi-lint: ok <reason>')")

    if problems:
        sys.stderr.write("post-edit check failed for %s:\n%s\n" % (rel, "\n".join(problems)))
        sys.exit(2)
    sys.exit(0)


if __name__ == "__main__":
    main()
