#!/usr/bin/env python3
"""
Analyse a goSPL scaling sweep.

Reads the ``scaling_p<N>.json`` records produced by ``run_scaling.py`` (one per
MPI rank count), computes strong-scaling speedup S(p)=T(p0)/T(p) and parallel
efficiency E(p)=S(p)·p0/p relative to the smallest rank count in the sweep, and
emits:

  * a console table,
  * ``results/scaling/scaling_summary.csv``  (raw numbers),
  * ``results/scaling/scaling_summary.md``   (Markdown table for reports),
  * ``results/scaling/scaling_speedup.png``  (speedup + efficiency vs p),
  * ``results/scaling/scaling_phases.png``   (per-phase mean wall vs p, stacked),

The plots are skipped (with a notice) if matplotlib is unavailable, so the CSV
and Markdown always get written. This script is serial — run it on a login node
after the sweep finishes.

Usage::

    python analyze_scaling.py scaling_results/ -o results/scaling
"""

import os
import csv
import glob
import json
import argparse


def _load(indir):
    records = []
    for path in sorted(glob.glob(os.path.join(indir, "scaling_p*.json"))):
        with open(path) as fh:
            records.append(json.load(fh))
    records.sort(key=lambda r: r["nranks"])
    return records


def _phase_union(records):
    names, seen = [], set()
    for r in records:
        for n in r["phases"]:
            if n not in seen:
                seen.add(n)
                names.append(n)
    # Stable, biggest-first by total mean time across the sweep.
    weight = {
        n: sum(r["phases"].get(n, {}).get("mean", 0.0) for r in records) for n in names
    }
    names.sort(key=lambda n: weight[n], reverse=True)
    return names


def analyse(records):
    base = records[0]
    p0, t0 = base["nranks"], base["run_wall_s"]
    rows = []
    for r in records:
        p, t = r["nranks"], r["run_wall_s"]
        speedup = t0 / t if t > 0 else 0.0
        eff = speedup * p0 / p if p > 0 else 0.0
        rows.append(
            {
                "nranks": p,
                "run_wall_s": t,
                "speedup": speedup,
                "efficiency": eff,
                "rss_max_mb_per_rank": r.get("rss_max_mb_per_rank", float("nan")),
                "rss_sum_mb": r.get("rss_sum_mb", float("nan")),
            }
        )
    return rows, p0


def write_csv(rows, names, records, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        header = [
            "nranks",
            "run_wall_s",
            "speedup",
            "efficiency",
            "rss_max_mb_per_rank",
            "rss_sum_mb",
        ] + ["phase_%s_mean_s" % n for n in names] + [
            "phase_%s_imbal" % n for n in names
        ]
        w.writerow(header)
        for row, r in zip(rows, records):
            line = [
                row["nranks"],
                "%.4f" % row["run_wall_s"],
                "%.4f" % row["speedup"],
                "%.4f" % row["efficiency"],
                "%.1f" % row["rss_max_mb_per_rank"],
                "%.1f" % row["rss_sum_mb"],
            ]
            line += ["%.4f" % r["phases"].get(n, {}).get("mean", 0.0) for n in names]
            line += ["%.2f" % r["phases"].get(n, {}).get("imbalance", 0.0) for n in names]
            w.writerow(line)


def write_md(rows, p0, path):
    lines = [
        "# goSPL strong-scaling summary",
        "",
        "Baseline = %d rank(s). Speedup S(p)=T(base)/T(p); efficiency "
        "E(p)=S·base/p (1.00 = perfect)." % p0,
        "",
        "| ranks | wall (s) | speedup | efficiency | RSS/rank (MB) |",
        "|------:|---------:|--------:|-----------:|--------------:|",
    ]
    for row in rows:
        lines.append(
            "| %d | %.2f | %.2f | %.0f%% | %.0f |"
            % (
                row["nranks"],
                row["run_wall_s"],
                row["speedup"],
                100.0 * row["efficiency"],
                row["rss_max_mb_per_rank"],
            )
        )
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def plot(rows, names, records, outdir):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("[analyze] matplotlib/numpy not available — skipping plots")
        return

    ps = [r["nranks"] for r in rows]

    # --- speedup + efficiency -------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(7, 5))
    sp = [r["speedup"] for r in rows]
    ax1.plot(ps, sp, "o-", color="tab:blue", label="measured speedup")
    ax1.plot(
        ps,
        [p / ps[0] for p in ps],
        "k--",
        alpha=0.6,
        label="ideal (linear)",
    )
    ax1.set_xlabel("MPI ranks")
    ax1.set_ylabel("speedup", color="tab:blue")
    ax1.set_xscale("log", base=2)
    ax1.set_yscale("log", base=2)
    ax1.set_xticks(ps)
    ax1.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax1.grid(True, which="both", alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(ps, [100 * r["efficiency"] for r in rows], "s-", color="tab:red",
             label="efficiency")
    ax2.set_ylabel("parallel efficiency (%)", color="tab:red")
    ax2.set_ylim(0, 110)
    ax1.legend(loc="upper left")
    fig.suptitle("goSPL strong scaling")
    fig.tight_layout()
    p1 = os.path.join(outdir, "scaling_speedup.png")
    fig.savefig(p1, dpi=130)
    plt.close(fig)

    # --- per-phase mean wall, stacked ----------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    bottom = np.zeros(len(ps))
    x = np.arange(len(ps))
    for n in names:
        vals = np.array([r["phases"].get(n, {}).get("mean", 0.0) for r in records])
        ax.bar(x, vals, bottom=bottom, label=n)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([str(p) for p in ps])
    ax.set_xlabel("MPI ranks")
    ax.set_ylabel("mean wall time per phase (s)")
    ax.set_title("goSPL per-phase wall time vs ranks")
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    p2 = os.path.join(outdir, "scaling_phases.png")
    fig.savefig(p2, dpi=130)
    plt.close(fig)

    print("[analyze] wrote %s and %s" % (p1, p2))


def main(argv=None):
    parser = argparse.ArgumentParser(description="Analyse a goSPL scaling sweep")
    parser.add_argument("indir", help="dir with scaling_p*.json records")
    parser.add_argument("-o", "--outdir", default="results/scaling")
    args = parser.parse_args(argv)

    records = _load(args.indir)
    if not records:
        raise SystemExit("No scaling_p*.json found in %s" % args.indir)

    names = _phase_union(records)
    rows, p0 = analyse(records)
    os.makedirs(args.outdir, exist_ok=True)

    # Console table.
    print("\ngoSPL strong scaling (baseline = %d rank(s))" % p0)
    print("%6s %10s %9s %12s %14s" % ("ranks", "wall(s)", "speedup", "efficiency",
                                      "RSS/rank(MB)"))
    for row in rows:
        print(
            "%6d %10.2f %9.2f %11.0f%% %14.0f"
            % (
                row["nranks"],
                row["run_wall_s"],
                row["speedup"],
                100.0 * row["efficiency"],
                row["rss_max_mb_per_rank"],
            )
        )

    write_csv(rows, names, records, os.path.join(args.outdir, "scaling_summary.csv"))
    write_md(rows, p0, os.path.join(args.outdir, "scaling_summary.md"))
    plot(rows, names, records, args.outdir)
    print("\n[analyze] wrote CSV + Markdown to %s" % args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
