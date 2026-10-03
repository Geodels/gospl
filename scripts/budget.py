#!/usr/bin/env python3
"""Volume budget of a goSPL output directory (eroded / deposited / net per
output step). Thin wrapper over ``gospl-inspect <outdir> --budget``; see
gospl/analyse/runinspect.py. For the exact in-model per-step budget, run with
``gospl --summary run.jsonl`` and read its ``volume`` records.

    python scripts/budget.py <output-dir> [--json] [--mesh mesh.npz]
"""
import sys

from gospl.analyse.runinspect import main

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] + ["--budget"]))
