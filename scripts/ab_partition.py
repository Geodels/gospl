#!/usr/bin/env python3
"""Run one goSPL input at two rank counts and compare the fields node by node.

The tool for AGENTS.md > MPI contract "THE #2 partition-dependence class" and
the drainage-invariance work in gospl/flow/AGENTS.md: a correct parallel model
gives the same answer (to solver round-off) whatever the decomposition. This
runs the model at ``np=1`` (or ``--ref-np``) and at ``-n N``, gathers every
field to rank 0 in INPUT-MESH order (``locIDs``, the partition-invariant id --
never PETSc's per-partition numbering), and reports, per field, the largest
difference, where it is, and how many nodes differ.

Usage::

    python scripts/ab_partition.py input.yml -n 4
    python scripts/ab_partition.py input.yml -n 8 --ref-np 2 --steps 3
    python scripts/ab_partition.py input.yml -n 4 --field iceHL --field soilH --json
    python scripts/ab_partition.py input.yml -n 4 --keep ab_out   # keep the npz dumps

``--steps K`` runs only the first K time steps (default: 1). Each run executes
in its own scratch directory whose entries are symlinks to the input's
directory, so model output never lands next to your input. The comparison
ignores nodes on the domain boundary only if ``--interior`` is given.

Exit status: 0 if every field agrees within ``--rtol``/``--atol``, 1 otherwise,
so it can gate a script. A difference is a lead, not a verdict. The default
tolerances are strict on purpose, and a healthy model does NOT pass them
node-by-node: KSP round-off flips near-tie flow routing at partition seams.
Reference baseline (tests/fixtures/minimal.yml, np=1 vs np=2, 2 steps, macOS):

    elev    max rel 4e-6   rel L2 8e-7      (round-off level: healthy)
    FA      max rel 0.3    rel L2 8e-2      (near-tie re-routing at ~175 nodes)
    cumED   max rel 3e-2   rel L2 2e-2

What signals a real partition bug is a jump well above such a baseline (orders
of magnitude, or a localised spike along the partition boundary), a non-finite
count that differs between the runs, or a difference that GROWS with --steps.
Measure the baseline for your input first, then compare against it.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

# Fields gathered by default: (name, attribute, kind). kind "vec" = a PETSc
# LOCAL Vec (read with getArray), "arr" = a per-local-node numpy array.
DEFAULT_FIELDS = [
    ("elev", "hLocal", "vec"),
    ("cumED", "cumEDLocal", "vec"),
    ("FA", "FAL", "vec"),
    ("EDrate", "EbLocal", "vec"),
    ("sedLoad", "vSedLocal", "vec"),
    ("rain", "rainVal", "arr"),
]

DRIVER = r'''
import json, sys
import numpy as np
from mpi4py import MPI
from gospl.model import Model

yml, out, steps, extra = sys.argv[1], sys.argv[2], int(sys.argv[3]), json.loads(sys.argv[4])
model = Model(yml, verbose=False, showlog=False)
try:
    if steps > 0:
        model.tEnd = min(model.tEnd, model.tStart + (steps - 1) * model.dt)
    model.runProcesses()
    comm = MPI.COMM_WORLD
    owned = model.inIDs == 1
    gid = np.asarray(model.locIDs)[owned]
    fields = {}
    for name, attr, kind in extra:
        obj = getattr(model, attr, None)
        if obj is None:
            continue
        val = obj.getArray().copy() if kind == "vec" else np.asarray(obj, dtype=float)
        if val.ndim != 1 or val.shape[0] != model.lpoints:
            continue
        fields[name] = val[owned]
    parts = comm.gather((gid, fields), root=0)
    if comm.rank == 0:
        n = model.mpoints
        res = {"_coords": np.asarray(model.mCoords)}
        names = set().union(*(p[1].keys() for p in parts))
        for name in names:
            arr = np.full(n, np.nan)
            for g, f in parts:
                if name in f:
                    arr[g] = f[name]
            res[name] = arr
        np.savez(out, **res)
finally:
    model.destroy()
'''


def _child_env():
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("OMPI_", "PMIX_", "PRTE_", "OPAL_"))}
    if "OPAL_PREFIX" in os.environ:
        env["OPAL_PREFIX"] = os.environ["OPAL_PREFIX"]
    env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    env.setdefault("OMP_NUM_THREADS", "1")
    return env


def _output_dir_name(yml: Path) -> str | None:
    try:
        from ruamel.yaml import YAML
        data = YAML(typ="safe").load(yml.read_text()) or {}
        return str((data.get("output") or {}).get("dir", "output")).strip("/")
    except Exception:
        return "output"


def _scratch_copy(yml: Path, root: Path, tag: str) -> Path:
    """A directory mirroring the input's directory through symlinks."""
    d = root / tag
    d.mkdir()
    outdir = _output_dir_name(yml)
    for entry in yml.parent.iterdir():
        # never mirror the input's own output dir (or its _N siblings): the
        # model rmtree's + recreates it at step 0, which must happen HERE.
        if outdir and (entry.name == outdir or entry.name.startswith(outdir + "_")):
            continue
        try:
            (d / entry.name).symlink_to(entry.resolve())
        except OSError:
            pass
    # the YAML itself is copied (goSPL copies the input into its output dir)
    (d / yml.name).unlink(missing_ok=True)
    shutil.copy2(yml, d / yml.name)
    return d


def run(yml: Path, np_: int, steps: int, fields, root: Path, timeout: float) -> Path:
    d = _scratch_copy(yml, root, f"np{np_}")
    driver = root / "ab_driver.py"
    driver.write_text(DRIVER)
    out = d / f"ab_np{np_}.npz"
    launcher = ["mpirun", "-n", str(np_)] if np_ > 1 or shutil.which("mpirun") else []
    cmd = launcher + [sys.executable, str(driver), yml.name, str(out), str(steps),
                      json.dumps(fields)]
    res = subprocess.run(cmd, cwd=d, env=_child_env(), capture_output=True,
                         text=True, timeout=timeout)
    if res.returncode != 0 or not out.exists():
        sys.stderr.write(res.stdout[-4000:] + res.stderr[-4000:])
        raise SystemExit(f"ab_partition: np={np_} run failed (rc={res.returncode})")
    return out


def compare(a: dict, b: dict, rtol: float, atol: float, interior=None):
    coords = a["_coords"]
    report = []
    for name in sorted(k for k in a.files if not k.startswith("_")):
        if name not in b.files:
            continue
        x, y = a[name], b[name]
        ok = np.isfinite(x) & np.isfinite(y)
        if interior is not None:
            ok &= interior
        d = np.abs(x - y)
        d[~ok] = 0.0
        scale = max(np.abs(x[ok]).max(initial=0.0), 1e-300)
        bad = d > (atol + rtol * np.abs(x))
        bad &= ok
        k = int(np.argmax(d))
        report.append({
            "field": name,
            "max_abs_diff": float(d.max(initial=0.0)),
            "max_rel_diff": float(d.max(initial=0.0) / scale),
            "rel_l2": float(np.linalg.norm(d[ok]) / max(np.linalg.norm(x[ok]), 1e-300)),
            "n_differ": int(bad.sum()),
            "n_nodes": int(ok.sum()),
            "worst_node": k,                      # input-mesh index (npz row)
            "worst_xyz": [float(v) for v in coords[k]],
            "ref_value": float(x[k]),
            "value": float(y[k]),
            "nonfinite_ref": int((~np.isfinite(x)).sum()),
            "nonfinite": int((~np.isfinite(y)).sum()),
            "pass": bool(bad.sum() == 0 and np.isfinite(x).sum() == np.isfinite(y).sum()),
        })
    return report


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("input", type=Path, help="goSPL YAML input")
    ap.add_argument("-n", "--np", type=int, default=2, help="rank count to test (default 2)")
    ap.add_argument("--ref-np", type=int, default=1, help="reference rank count (default 1)")
    ap.add_argument("--steps", type=int, default=1, help="time steps to run (0 = full run)")
    ap.add_argument("--field", action="append", default=[], metavar="ATTR",
                    help="extra Model attribute to compare (local Vec or per-node array)")
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--atol", type=float, default=1e-8)
    ap.add_argument("--interior", action="store_true",
                    help="ignore nodes on the domain boundary (flat meshes)")
    ap.add_argument("--timeout", type=float, default=3600.0, help="per-run timeout (s)")
    ap.add_argument("--keep", type=Path, help="copy both npz dumps into this directory")
    ap.add_argument("--json", action="store_true", help="print the report as JSON")
    args = ap.parse_args(argv)

    yml = args.input.resolve()
    fields = [list(f) for f in DEFAULT_FIELDS]
    for attr in args.field:
        kind = "vec" if attr.endswith(("L", "Local")) else "arr"
        fields.append([attr, attr, kind])

    with tempfile.TemporaryDirectory(prefix="gospl-ab-") as tmp:
        root = Path(tmp)
        ref = run(yml, args.ref_np, args.steps, fields, root, args.timeout)
        tst = run(yml, args.np, args.steps, fields, root, args.timeout)
        a, b = np.load(ref), np.load(tst)
        interior = None
        if args.interior:
            c = a["_coords"]
            lo, hi = c.min(0), c.max(0)
            axes = (hi - lo) > 0.0            # a flat mesh has a constant z
            span = (hi - lo)[axes]
            cc = c[:, axes]
            interior = np.all((cc - lo[axes] > 1e-6 * span)
                              & (hi[axes] - cc > 1e-6 * span), axis=1)
        report = compare(a, b, args.rtol, args.atol, interior)
        if args.keep:
            args.keep.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ref, args.keep / ref.name)
            shutil.copy2(tst, args.keep / tst.name)

    if args.json:
        print(json.dumps({"input": str(yml), "ref_np": args.ref_np, "np": args.np,
                          "steps": args.steps, "fields": report}, indent=1))
    else:
        print(f"{yml.name}: np={args.ref_np} vs np={args.np}, {args.steps or 'all'} step(s)")
        print(f"{'field':10s} {'max|d|':>11s} {'max rel':>9s} {'rel L2':>9s} "
              f"{'#differ':>8s}  worst node (input-mesh id)")
        for r in report:
            flag = "" if r["pass"] else "  <-- DIFFERS"
            print(f"{r['field']:10s} {r['max_abs_diff']:11.3e} {r['max_rel_diff']:9.2e} "
                  f"{r['rel_l2']:9.2e} {r['n_differ']:8d}  {r['worst_node']}{flag}")
    return 0 if all(r["pass"] for r in report) else 1


if __name__ == "__main__":
    sys.exit(main())
