#!/usr/bin/env python
"""
List the methods of goSPL's model classes that are missing from the API
reference (``docs/api_ref/*_ref.rst``).

Every method of a class documented with ``.. autoclass::`` must appear in the
page's ``.. autosummary::`` list and have an ``.. automethod::`` entry (the
lists are written by hand, AGENTS.md > Docs). Exit code 1 when some are
missing, so ``tests/test_contracts.py`` fails on a method added without its
API entry.

    python scripts/check_api_docs.py
"""

import ast
import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# classes documented with :members: (whole classes) or not model mixins
SKIP_CLASSES = {"_IndentWriter", "GridBuilder", "ProvenanceTracker", "GosplOutput",
                "Output", "Profiler", "RunSummary"}
SKIP_METHODS = {"_phase"}


def missing():
    pages = {p: open(p).read() for p in glob.glob(os.path.join(ROOT, "docs", "api_ref", "*_ref.rst"))}
    out = []
    for f in sorted(glob.glob(os.path.join(ROOT, "gospl", "**", "*.py"), recursive=True)):
        tree = ast.parse(open(f).read())
        for c in [n for n in tree.body if isinstance(n, ast.ClassDef)]:
            if c.name in SKIP_CLASSES:
                continue
            page = next((s for s in pages.values()
                         if re.search(r"autoclass:: [\w.]*\b%s\b" % c.name, s)), None)
            if page is None:
                continue
            for m in c.body:
                if not isinstance(m, ast.FunctionDef) or m.name.startswith("__") or m.name in SKIP_METHODS:
                    continue
                if any(getattr(d, "id", getattr(d, "attr", "")) in ("property", "setter")
                       for d in m.decorator_list):
                    continue
                listed = re.search(r"~%s\.%s\s*$" % (c.name, m.name), page, re.M)
                auto = re.search(r"automethod:: [\w.]*%s\.%s\s*$" % (c.name, m.name), page, re.M)
                if not (listed and auto):
                    out.append("%s %s.%s" % (os.path.relpath(f, ROOT), c.name, m.name))
    return out


if __name__ == "__main__":
    miss = missing()
    for m in miss:
        print("missing from the API reference:", m)
    sys.exit(1 if miss else 0)
