#!/usr/bin/env python
"""
Build the combined HPC strong-scaling figure for the user guide from the
``scaling_<N>km.csv`` summaries in this folder (one per global mesh resolution).

Each CSV is a goSPL scaling sweep (``scripts/scaling/analyze_scaling.py`` output:
columns ``nranks, run_wall_s, speedup, efficiency, …``). Speedup and efficiency
are normalised to the SMALLEST rank count in each sweep (no 1-CPU run -- a global
mesh is memory-infeasible at low core counts), so the ideal line is ``p / p0``.

Drop in ``scaling_8km.csv`` / ``scaling_5km.csv`` and re-run; the figure picks up
every resolution automatically. Output: ``docs/images/scaling_hpc.png``.

    python docs/user_guide/scaling/make_scaling_figure.py
"""

import os
import csv
import glob

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "..", "images", "scaling_hpc.png"))

# Per-resolution metadata (node count for the legend). cell_width (km) -> label.
MESH = {5: "5 km · 23.6 M nodes", 8: "8 km · 9.2 M nodes", 10: "10 km · 5.9 M nodes"}
# Distinct, colour-blind-friendly styling per resolution.
STYLE = {
    5:  dict(color="#1b9e77", marker="o"),
    8:  dict(color="#d95f02", marker="s"),
    10: dict(color="#7570b3", marker="^"),
}

# Publication-quality defaults (no LaTeX dependency -- safe on Read the Docs).
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "legend.fontsize": 10.5,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.6,
    "axes.axisbelow": True,
    "savefig.dpi": 200,
    "figure.dpi": 120,
})


def _load(path):
    rows = list(csv.DictReader(open(path)))
    rows.sort(key=lambda r: int(r["nranks"]))
    p = [int(r["nranks"]) for r in rows]
    sp = [float(r["speedup"]) for r in rows]
    eff = [float(r["efficiency"]) for r in rows]
    return p, sp, eff


def main():
    files = sorted(glob.glob(os.path.join(HERE, "scaling_*km.csv")),
                   key=lambda f: -int("".join(c for c in os.path.basename(f) if c.isdigit())))
    if not files:
        raise SystemExit("no scaling_<N>km.csv files in %s" % HERE)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2))
    all_p = set()
    baselines = set()

    for f in files:
        km = int("".join(c for c in os.path.basename(f) if c.isdigit()))
        p, sp, eff = _load(f)
        all_p.update(p)
        baselines.add(p[0])
        st = STYLE.get(km, dict(color="0.3", marker="D"))
        lbl = MESH.get(km, "%d km" % km)
        ax1.plot(p, sp, "-", **st, label=lbl, lw=1.8, ms=6, mec="white", mew=0.6)
        ax2.plot(p, eff, "-", **st, label=lbl, lw=1.8, ms=6, mec="white", mew=0.6)

    pmin, pmax = min(all_p), max(all_p)

    # --- speedup (log-log) with the ideal p/p0 reference(s) ---
    for p0 in sorted(baselines):
        xs = [pmin, pmax]
        ax1.plot(xs, [x / p0 for x in xs], "--", color="0.55", lw=1.0, zorder=0,
                 label="ideal (p / %d)" % p0)
    ax1.set_xscale("log", base=2)
    ax1.set_yscale("log", base=2)
    ticks = sorted(all_p)
    ax1.set_xticks(ticks)
    ax1.set_xticklabels([str(t) for t in ticks])
    yt = [t / min(baselines) for t in ticks]
    ax1.set_yticks(yt)
    ax1.set_yticklabels(["%g" % v for v in yt])
    ax1.set_xlabel("MPI ranks (cores)")
    ax1.set_ylabel("speedup  (relative to smallest run)")
    ax1.set_title("(a) Strong-scaling speedup")
    ax1.legend(frameon=True, framealpha=0.9)

    # --- parallel efficiency (semilog-x) ---
    ax2.axhline(1.0, ls="--", color="0.55", lw=1.0, zorder=0, label="ideal (100 %)")
    ax2.set_xscale("log", base=2)
    ax2.set_xticks(ticks)
    ax2.set_xticklabels([str(t) for t in ticks])
    ax2.set_ylim(0.0, 1.15)
    ax2.set_xlabel("MPI ranks (cores)")
    ax2.set_ylabel("parallel efficiency")
    ax2.set_title("(b) Parallel efficiency")
    ax2.legend(frameon=True, framealpha=0.9)

    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
