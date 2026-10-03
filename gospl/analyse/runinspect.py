"""
Inspect a goSPL output directory (or a ``--summary`` JSON-Lines file).

``gospl-inspect`` answers the first questions asked of any run, without a
notebook or ad-hoc ``h5py`` script, and prints them as text or ``--json`` (one
blob that can be pasted from an HPC system):

* which steps / times / fields / partitions the output holds;
* per-field min / max / mean and the NaN/inf count at a step;
* the volume budget of every output step (eroded, deposited, net, from the
  change in ``erodep`` times the barycentric cell area). Outputs are float32
  and barycentric areas only approximate the model's Voronoi ``larea``, so a
  closed domain closes to ~1e-3 here; the exact in-model per-step budget is in
  the ``--summary`` file (``volume``);
* the largest elevation change between consecutive outputs and WHERE it is
  (input-mesh node id and coordinates): the usual first sign of a spike;
* ``--compare OTHER``: field-by-field differences against a second run at the
  same step (a code change A/B, or np=1 vs np=N output);
* ``--summary run.jsonl``: anomalies in a run summary written with
  ``gospl --summary`` (non-finite values, solver failures, aborts or ponded
  regions, stalled cascades, elevation jumps, budget drift).

Node ids are INPUT-MESH ids when the mesh can be found (the YAML copied into the
output directory names it via ``domain: npdata``; pass ``--mesh`` otherwise);
without the mesh they index a de-duplicated reassembly and are marked as such.

Usage::

    gospl-inspect outdir                       # overview + budget + spikes
    gospl-inspect outdir --step 12 --json
    gospl-inspect outdir --compare other_outdir
    gospl-inspect --summary run.jsonl
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

import numpy as np


# ---------------------------------------------------------------------------
# Reassembly
# ---------------------------------------------------------------------------


def _find_mesh(outdir):
    """Locate the run's npdata mesh from the YAML goSPL copies into outdir."""
    ymls = sorted(glob.glob(os.path.join(outdir, "*.yml")) +
                  glob.glob(os.path.join(outdir, "*.yaml")))
    for y in ymls:
        try:
            from ruamel.yaml import YAML

            dom = (YAML(typ="safe").load(open(y)) or {}).get("domain", {})
            spec = dom.get("npdata")
            if not spec:
                continue
            name = spec[0] if isinstance(spec, (list, tuple)) else spec
            name = name if name.endswith(".npz") else name + ".npz"
            # the run's cwd was the YAML's original directory: outdir's parent
            for base in (os.path.dirname(os.path.abspath(outdir)), outdir, os.getcwd()):
                cand = os.path.join(base, name)
                if os.path.exists(cand):
                    key = spec[1] if isinstance(spec, (list, tuple)) and len(spec) > 1 else "v"
                    return cand, key
        except Exception:
            continue
    return None, None


class Output(object):
    """A goSPL output directory reassembled onto one global node ordering."""

    def __init__(self, outdir, mesh=None, file_base="gospl"):
        import h5py

        self.outdir = outdir
        self.h5dir = os.path.join(outdir, "h5")
        self.file_base = file_base
        tfiles = sorted(glob.glob(os.path.join(self.h5dir, "topology.p*.h5")),
                        key=lambda f: int(re.search(r"\.p(\d+)\.h5$", f).group(1)))
        if not tfiles:
            raise FileNotFoundError("no h5/topology.p*.h5 under %s" % outdir)
        self.nparts = len(tfiles)
        lcoords, lcells = [], []
        for tf in tfiles:
            with h5py.File(tf, "r") as f:
                lcoords.append(np.asarray(f["coords"], dtype=np.float64))
                lcells.append(np.asarray(f["cells"], dtype=np.int64) - 1)  # 1-based

        mesh_path, key = (mesh, "v") if mesh else _find_mesh(outdir)
        if mesh_path:
            gcoords = np.asarray(np.load(mesh_path)[key], dtype=np.float64)
            from scipy.spatial import cKDTree

            tree = cKDTree(gcoords)
            self.maps = [tree.query(c)[1].astype(np.int64) for c in lcoords]
            self.input_ids = True
        else:
            allc = np.concatenate(lcoords)
            # shared (ghost) nodes are written from the same float64 value in
            # every partition, so their float32 coordinates match exactly
            gcoords, inv = np.unique(allc.astype(np.float32), axis=0,
                                     return_inverse=True)
            inv = inv.ravel()
            off = np.cumsum([0] + [len(c) for c in lcoords])
            self.maps = [inv[off[i]:off[i + 1]] for i in range(len(lcoords))]
            real = np.zeros(gcoords.shape, dtype=np.float64)
            real[inv] = allc
            gcoords = real
            self.input_ids = False
        self.mesh_path = mesh_path
        self.coords = gcoords
        self.npoints = len(gcoords)

        tri = np.concatenate([m[c] for m, c in zip(self.maps, lcells)])
        tri = np.unique(np.sort(tri, axis=1), axis=0)
        self.cells = tri
        p = gcoords[tri]
        a = 0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)
        self.area = np.bincount(tri.ravel(), weights=np.repeat(a / 3.0, 3),
                                minlength=self.npoints)

        steps = set()
        for f in glob.glob(os.path.join(self.h5dir, "%s.*.p0.h5" % file_base)):
            m = re.search(r"\.(\d+)\.p0\.h5$", f)
            if m:
                steps.add(int(m.group(1)))
        self.steps = sorted(steps)

    def fields(self, step):
        import h5py

        with h5py.File(self._file(step, 0), "r") as f:
            return sorted(k for k in f.keys() if f[k].ndim <= 2)

    def _file(self, step, part):
        return os.path.join(self.h5dir, "%s.%d.p%d.h5" % (self.file_base, step, part))

    def field(self, step, name):
        import h5py

        out = np.full(self.npoints, np.nan)
        for part, m in enumerate(self.maps):
            with h5py.File(self._file(step, part), "r") as f:
                if name not in f:
                    return None
                v = np.asarray(f[name], dtype=np.float64)
            if v.ndim > 1:
                if v.shape[1] != 1:
                    return None          # multi-column field (e.g. per-class)
                v = v[:, 0]
            out[m] = v
        return out

    def time(self, step):
        xmf = os.path.join(self.outdir, "xmf", "%s%d.xmf" % (self.file_base, step))
        try:
            m = re.search(r'<Time[^>]*Value="([^"]+)"', open(xmf).read())
            return float(m.group(1)) if m else None
        except OSError:
            return None

    def where(self, i):
        return {"node": int(i), "input_mesh_id": bool(self.input_ids),
                "xyz": [float(v) for v in self.coords[i]]}


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------


def field_stats(out, step):
    res = {}
    for name in out.fields(step):
        v = out.field(step, name)
        if v is None:
            continue
        ok = np.isfinite(v)
        res[name] = {
            "min": float(v[ok].min()) if ok.any() else None,
            "max": float(v[ok].max()) if ok.any() else None,
            "mean": float(v[ok].mean()) if ok.any() else None,
            "nonfinite": int((~ok).sum()),
        }
    return res


def budget(out):
    """Per-output-step eroded / deposited / net volume from erodep (m^3)."""
    rows, prev = [], None
    for s in out.steps:
        ed = out.field(s, "erodep")
        if ed is None:
            return []
        if prev is not None:
            d = np.where(np.isfinite(ed - prev), ed - prev, 0.0) * out.area
            ero, dep = float(d[d < 0].sum()), float(d[d > 0].sum())
            act = dep - ero
            rows.append({"step": s, "t": out.time(s), "eroded": ero, "deposited": dep,
                         "net": ero + dep,
                         "net_rel": (ero + dep) / act if act > 0 else 0.0})
        prev = ed
    return rows


def spikes(out, top=1):
    """Largest |elev change| between consecutive outputs, with location."""
    rows, prev = [], None
    for s in out.steps:
        z = out.field(s, "elev")
        if z is None:
            return []
        if prev is not None:
            d = np.abs(z - prev)
            d[~np.isfinite(d)] = -1.0
            i = int(np.argmax(d))
            rows.append({"step": s, "t": out.time(s), "max_abs_dz": float(d[i]),
                         "z_before": float(prev[i]), "z_after": float(z[i]),
                         **out.where(i)})
        prev = z
    return rows


def compare(a, b, step):
    if a.npoints != b.npoints:
        raise SystemExit("meshes differ: %d vs %d nodes" % (a.npoints, b.npoints))
    res = {}
    for name in sorted(set(a.fields(step)) & set(b.fields(step))):
        x, y = a.field(step, name), b.field(step, name)
        if x is None or y is None:
            continue
        ok = np.isfinite(x) & np.isfinite(y)
        d = np.where(ok, np.abs(x - y), 0.0)
        i = int(np.argmax(d))
        nx = np.linalg.norm(x[ok])
        res[name] = {"max_abs_diff": float(d[i]),
                     "rel_l2": float(np.linalg.norm(d[ok]) / nx) if nx > 0 else 0.0,
                     "nonfinite": [int((~np.isfinite(x)).sum()), int((~np.isfinite(y)).sum())],
                     **a.where(i)}
    return res


def inspect_summary(path, jump=0.25):
    """Anomalies in a `gospl --summary` JSON-Lines file."""
    recs = [json.loads(l) for l in open(path) if l.strip()]
    head = next((r for r in recs if r.get("start")), {})
    end = next((r for r in recs if r.get("end")), None)
    steps = [r for r in recs if "step" in r and "elev" in r]
    issues = []
    prev = None
    for r in steps:
        t = r["t"]
        nf = r.get("nonfinite", {})
        if any(v > 0 for v in nf.values()):
            issues.append({"t": t, "issue": "nonfinite", "detail": nf})
        rescued = {e.get("solve") for e in r.get("events", [])
                   if e.get("kind") == "flow_ksp_exact_rescue"}
        for k, s in r.get("ksp", {}).items():
            # A primary failure the exact block solver then recovered is the
            # designed long-drainage-chain path, not an issue.
            if s.get("fails") and k not in rescued and not k.endswith("_exact"):
                issues.append({"t": t, "issue": "ksp_primary_failed", "solver": k,
                               "detail": s})
        for e in r.get("events", []):
            kind = e.get("kind")
            # "benign" = a small isolated pocket left to pond: the documented,
            # expected outcome (AGENTS.md > Flow accumulation), not an issue.
            if kind == "flow_ksp_fallback_failed" and e.get("outcome") == "benign":
                continue
            if kind == "flow_ksp_fallback_failed" or kind == "soil_substep" or \
               kind.endswith("_failed") or (kind == "flow_cascade" and
                                            e.get("outcome") in ("stall", "max_steps")):
                issues.append({"t": t, "issue": kind, "detail": e})
        if prev is not None:
            # A needle moves an EXTREME far more than the mean; uniform uplift
            # or subsidence moves both. Flag the former only.
            rng = max(prev["elev"]["max"] - prev["elev"]["min"], 1.0)
            dz = max(abs(r["elev"]["max"] - prev["elev"]["max"]),
                     abs(r["elev"]["min"] - prev["elev"]["min"]))
            dmean = abs(r["elev"]["mean"] - prev["elev"]["mean"])
            if dz > jump * rng and dz > 5.0 * dmean:
                issues.append({"t": t, "issue": "elev_extreme_jump",
                               "detail": {"before": prev["elev"], "after": r["elev"]}})
        vol = r.get("volume", {})
        act = vol.get("deposited", 0.0) - vol.get("eroded", 0.0)
        # Only a CLOSED domain must balance; on an open one a negative net is
        # sediment leaving through the boundary. A gain is never physical.
        net_rel = vol.get("net", 0.0) / act if act > 0 else 0.0
        if act > 0 and (net_rel > 0.05 or (head.get("closed") and abs(net_rel) > 0.05)):
            issues.append({"t": t, "issue": "volume_imbalance",
                           "detail": {**vol, "net_rel": net_rel}})
        prev = r
    wall = [r["wall_s"] for r in steps]
    phase_tot = {}
    for r in steps:
        for k, v in r.get("phases", {}).items():
            phase_tot[k] = phase_tot.get(k, 0.0) + v
    return {
        "file": path, "run": head, "steps": len(steps),
        "finished": end is not None,
        "t_last": steps[-1]["t"] if steps else None,
        "wall_total_s": end["wall_s"] if end else float(np.sum(wall)) if wall else 0.0,
        "wall_step_s": {"mean": float(np.mean(wall)), "max": float(np.max(wall))} if wall else {},
        "phase_total_s": dict(sorted(phase_tot.items(), key=lambda kv: -kv[1])[:12]),
        "issues": issues,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _fmt(v):
    return "n/a" if v is None else ("%.4g" % v)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="gospl-inspect",
                                 description="Inspect a goSPL output directory or run summary.")
    ap.add_argument("outdir", nargs="?", help="goSPL output directory (the YAML `output: dir`)")
    ap.add_argument("--step", type=int, help="output step for the field table (default: last)")
    ap.add_argument("--compare", metavar="OTHER", help="second output dir to diff against")
    ap.add_argument("--summary", metavar="JSONL", help="inspect a `gospl --summary` file")
    ap.add_argument("--mesh", help="input npz mesh (if it cannot be found from the YAML)")
    ap.add_argument("--budget", action="store_true", help="only print the volume budget")
    ap.add_argument("--json", action="store_true", help="print one JSON object")
    args = ap.parse_args(argv)

    if args.summary:
        rep = inspect_summary(args.summary)
        if args.json:
            print(json.dumps(rep, indent=1))
        else:
            print("%s: %d steps, last t=%s, %s, wall %.1f s" % (
                rep["file"], rep["steps"], rep["t_last"],
                "finished" if rep["finished"] else "NOT FINISHED", rep["wall_total_s"]))
            print("slowest phases (s):", ", ".join(
                "%s %.1f" % kv for kv in list(rep["phase_total_s"].items())[:6]))
            if not rep["issues"]:
                print("no issues found")
            for i in rep["issues"]:
                print("  t=%-12s %-26s %s" % (i["t"], i["issue"], json.dumps(i["detail"])[:160]))
        return 1 if any(i["issue"] in ("nonfinite",) for i in rep["issues"]) else 0

    if not args.outdir:
        ap.error("an output directory or --summary is required")
    out = Output(args.outdir, mesh=args.mesh)
    if not out.steps:
        raise SystemExit("no %s.<step>.p0.h5 outputs under %s/h5" % (out.file_base, args.outdir))
    step = args.step if args.step is not None else out.steps[-1]

    rep = {"outdir": args.outdir, "nparts": out.nparts, "npoints": out.npoints,
           "mesh": out.mesh_path, "input_mesh_ids": out.input_ids,
           "steps": out.steps, "times": [out.time(s) for s in out.steps],
           "step": step}
    rep["budget"] = budget(out)
    if not args.budget:
        rep["fields"] = field_stats(out, step)
        rep["spikes"] = spikes(out)
    if args.compare:
        rep["compare"] = compare(out, Output(args.compare, mesh=args.mesh), step)

    if args.json:
        print(json.dumps(rep, indent=1))
        return 0

    print("%s: %d node(s), %d partition(s), steps %s..%s (t=%s..%s)%s" % (
        args.outdir, out.npoints, out.nparts, out.steps[0], out.steps[-1],
        _fmt(rep["times"][0]), _fmt(rep["times"][-1]),
        "" if out.input_ids else "  [node ids: reassembly order, mesh not found]"))
    if "fields" in rep:
        print("\nfields at step %d:" % step)
        print("  %-14s %12s %12s %12s %9s" % ("name", "min", "max", "mean", "nonfinite"))
        for k, v in rep["fields"].items():
            flag = "  <-- NON-FINITE" if v["nonfinite"] else ""
            print("  %-14s %12s %12s %12s %9d%s" % (k, _fmt(v["min"]), _fmt(v["max"]),
                                                     _fmt(v["mean"]), v["nonfinite"], flag))
    if rep["budget"]:
        print("\nvolume budget between outputs (m^3, barycentric areas):")
        print("  %6s %12s %12s %12s %12s %9s" % ("step", "t", "eroded", "deposited", "net", "net/act"))
        for r in rep["budget"]:
            print("  %6d %12s %12.4e %12.4e %12.4e %9.2e" % (
                r["step"], _fmt(r["t"]), r["eroded"], r["deposited"], r["net"], r["net_rel"]))
    if rep.get("spikes"):
        worst = max(rep["spikes"], key=lambda r: r["max_abs_dz"])
        print("\nlargest elevation change between outputs: %.4g m at step %d, node %d "
              "(%s), z %.4g -> %.4g" % (worst["max_abs_dz"], worst["step"], worst["node"],
                                        ", ".join("%.6g" % c for c in worst["xyz"]),
                                        worst["z_before"], worst["z_after"]))
    if "compare" in rep:
        print("\ncompare vs %s at step %d:" % (args.compare, step))
        for k, v in rep["compare"].items():
            print("  %-14s max|d| %11.4e  rel L2 %9.2e  worst node %d" % (
                k, v["max_abs_diff"], v["rel_l2"], v["node"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
