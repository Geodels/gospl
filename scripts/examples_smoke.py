#!/usr/bin/env python3
"""Smoke-run every goSPL-examples input for a few time steps.

The examples repository (github.com/Geodels/goSPL-examples) is the main way users
meet the model, and it has broken several times because of model changes it was
never run against (the NESW ``bc`` order, the soil solver default, open-edge
draining). This finds every input YAML in a checkout, runs each one for
``--steps`` time steps under ``mpirun -n --np`` in a private scratch directory
(symlinks to the example's files, so nothing is written into the checkout), with
``--summary`` on, and reports per example: passed / FAILED / missing-data,
wall time, and any anomaly in its run summary (NaN/inf, solver failures, ponded
regions, stalled cascades).

Usage::

    python scripts/examples_smoke.py ~/Workspace/goSPL-examples
    python scripts/examples_smoke.py <dir> --only Local-examples --np 2 --steps 2
    python scripts/examples_smoke.py <dir> --list
    python scripts/examples_smoke.py <dir> --json report.json

An input is any ``*.yml`` that is not ``environment*.yml`` and does not sit next
to an ``h5/`` directory (goSPL copies the input into its output directory).
"missing-data" means the run stopped on a missing input file -- typically data a
notebook or another example generates (e.g. continental_geochem reads
continental_flux's output) -- and does not count as a failure. Exit status: 1 if
any example FAILED, else 0.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

MISSING_RE = re.compile(r"FileNotFoundError|No such file or directory|"
                        r"does not exist|not present", re.I)


def find_inputs(root: Path, only: str | None):
    out = []
    for y in sorted(root.rglob("*.yml")):
        rel = y.relative_to(root)
        if any(p.startswith(".") for p in rel.parts) or y.name.startswith("environment"):
            continue
        if (y.parent / "h5").is_dir():
            continue                          # an output copy of an input
        if only and only not in str(rel):
            continue
        out.append(y)
    return out


def _child_env():
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("OMPI_", "PMIX_", "PRTE_", "OPAL_"))}
    if "OPAL_PREFIX" in os.environ:
        env["OPAL_PREFIX"] = os.environ["OPAL_PREFIX"]
    env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env.setdefault(k, "1")
    return env


def prepare(yml: Path, scratch: Path, steps: int) -> Path:
    """Mirror the example dir into `scratch` and write a shortened input."""
    from ruamel.yaml import YAML

    yaml = YAML()                              # round-trip: keeps comments/order
    data = yaml.load(yml.read_text())
    outdir = str((data.get("output") or {}).get("dir", "output"))
    for entry in yml.parent.iterdir():
        if entry.name == outdir or entry.name.startswith(outdir + "_") or \
                (entry.is_dir() and (entry / "h5").is_dir()):
            continue                          # never mirror output dirs
        try:
            (scratch / entry.name).symlink_to(entry.resolve())
        except OSError:
            pass
    t = data["time"]
    start, dt = float(t["start"]), float(t["dt"])
    t["end"] = start + steps * dt
    t["tout"] = dt
    for k in ("rstep",):
        t.pop(k, None)
    data.setdefault("output", {})
    data["output"]["dir"] = "smoke_out"
    data["output"]["makedir"] = False
    smoke = scratch / ("smoke_" + yml.name)
    with open(smoke, "w") as f:
        yaml.dump(data, f)
    return smoke


def run_one(yml: Path, root: Path, args) -> dict:
    rel = str(yml.relative_to(root))
    rec = {"input": rel, "status": None, "wall_s": None, "issues": []}
    with tempfile.TemporaryDirectory(prefix="gospl-smoke-") as tmp:
        scratch = Path(tmp)
        try:
            smoke = prepare(yml, scratch, args.steps)
        except Exception as exc:              # unparsable / unusual input
            rec.update(status="FAILED", error="prepare: %r" % exc)
            return rec
        summary = scratch / "run.jsonl"
        launcher = ["mpirun", "-n", str(args.np)] if shutil.which("mpirun") else []
        cmd = launcher + [sys.executable, "-m", "gospl", "-i", smoke.name,
                          "--summary", str(summary)]
        t0 = time.time()
        try:
            res = subprocess.run(cmd, cwd=scratch, env=_child_env(), capture_output=True,
                                 text=True, timeout=args.timeout)
            rc, out = res.returncode, res.stdout + res.stderr
        except subprocess.TimeoutExpired as exc:
            rc = -9
            out = (exc.stdout or b"").decode(errors="replace") if isinstance(exc.stdout, bytes) \
                else (exc.stdout or "")
            out += "\n[smoke] TIMEOUT after %.0f s" % args.timeout
        rec["wall_s"] = round(time.time() - t0, 1)
        if rc == 0 and summary.exists():
            rec["status"] = "passed"
            try:
                from gospl.analyse.runinspect import inspect_summary

                rep = inspect_summary(str(summary))
                rec["steps"] = rep["steps"]
                rec["issues"] = rep["issues"]
                if any(i["issue"] == "nonfinite" for i in rep["issues"]):
                    rec["status"] = "FAILED"
                    rec["error"] = "non-finite values in the run summary"
            except Exception as exc:
                rec["issues"] = [{"issue": "summary-unreadable", "detail": repr(exc)}]
        elif MISSING_RE.search(out):
            rec["status"] = "missing-data"
            rec["error"] = _tail_error(out)
        else:
            rec["status"] = "FAILED"
            rec["error"] = _tail_error(out)
    return rec


def _tail_error(out: str, n=12) -> str:
    lines = [l for l in out.splitlines() if l.strip()]
    err = [i for i, l in enumerate(lines) if "Error" in l or "error" in l or "TIMEOUT" in l]
    if err:
        i = err[-1]
        return "\n".join(lines[max(0, i - n + 1): i + 1])
    return "\n".join(lines[-n:])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("examples", type=Path, help="goSPL-examples checkout")
    ap.add_argument("--only", help="substring filter on the input path (e.g. Local-examples)")
    ap.add_argument("--np", type=int, default=2, help="MPI ranks per run (default 2)")
    ap.add_argument("--steps", type=int, default=2, help="time steps per run (default 2)")
    ap.add_argument("--timeout", type=float, default=900.0, help="per-example timeout (s)")
    ap.add_argument("--list", action="store_true", help="list the inputs and exit")
    ap.add_argument("--json", type=Path, help="write the report to this file")
    args = ap.parse_args(argv)

    root = args.examples.resolve()
    inputs = find_inputs(root, args.only)
    if args.list:
        print("\n".join(str(y.relative_to(root)) for y in inputs))
        return 0
    if not inputs:
        raise SystemExit("no example inputs under %s" % root)

    report = []
    for y in inputs:
        rec = run_one(y, root, args)
        report.append(rec)
        flag = {"passed": "ok", "missing-data": "skip", "FAILED": "FAIL"}[rec["status"]]
        extra = ("  (%d issue(s))" % len(rec["issues"])) if rec.get("issues") else ""
        print("[%4s] %-55s %7ss%s" % (flag, rec["input"], rec["wall_s"], extra), flush=True)
        if rec["status"] == "FAILED":
            print("       " + rec.get("error", "").replace("\n", "\n       "), flush=True)

    n = {s: sum(r["status"] == s for r in report) for s in ("passed", "missing-data", "FAILED")}
    print("\n%d passed, %d missing-data, %d FAILED (np=%d, %d step(s) each)" % (
        n["passed"], n["missing-data"], n["FAILED"], args.np, args.steps))
    if args.json:
        args.json.write_text(json.dumps({"np": args.np, "steps": args.steps,
                                         "results": report}, indent=1))
    return 1 if n["FAILED"] else 0


if __name__ == "__main__":
    sys.exit(main())
