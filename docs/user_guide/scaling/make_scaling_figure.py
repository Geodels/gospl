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
OUT_PHASES = os.path.normpath(
    os.path.join(HERE, "..", "..", "images", "scaling_phases.png"))

# Per-phase scaling figure: the TOP-LEVEL phases that partition the wall-clock
# (the `flow_*` sub-phases nest inside `flow` and would double-count). `flexure`
# is the genuinely flat serial floor; `flow_fill` (pit filling) scales partway
# -- the local priority-flood fill is parallel -- then flattens toward its
# serial rank-0 spillover-graph component, so it is broken out too. (column,
# label).
PHASES = [
    ("phase_sed_mean_s", "sediment routing"),
    ("phase_sea_mean_s", "marine deposition"),
    ("phase_flow_mean_s", "flow accumulation"),
    ("phase_erosion_mean_s", "river incision (SPL)"),
    ("phase_hillslope_mean_s", "hillslope diffusion"),
    ("phase_flow_fill_mean_s", "pit filling"),
    ("phase_flexure_mean_s", "flexure (serial)"),
]

# Per-resolution metadata (node count for the legend). cell_width (km) -> label.
MESH = {5: "5 km · 23.7 M nodes", 8: "8 km · 9.2 M nodes", 10: "10 km · 5.9 M nodes"}
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


def _label_ticks(ticks, ratio=1.4):
    """Thin a dense rank list to a well-spaced label subset (≥`ratio` apart on
    the log axis), always keeping the first and last. Markers are still plotted
    at every rank — this only controls which ticks get a printed label."""
    ticks = sorted(set(ticks))
    out = [ticks[0]]
    for t in ticks[1:]:
        if t >= out[-1] * ratio:
            out.append(t)
    if out[-1] != ticks[-1]:
        if ticks[-1] < out[-1] * ratio and len(out) > 1:
            out[-1] = ticks[-1]   # too close to the max -> replace, don't crowd
        else:
            out.append(ticks[-1])
    return out


def _load_phases(path):
    rows = list(csv.DictReader(open(path)))
    rows.sort(key=lambda r: int(r["nranks"]))
    p = [int(r["nranks"]) for r in rows]
    series = {}
    for col, _ in PHASES:
        if col in rows[0]:
            series[col] = [float(r[col]) for r in rows]
    return p, series


def phase_figure(files):
    """Per-phase strong scaling: one panel per mesh, log-log phase time vs
    ranks, with a 1/p ideal guide. Shows the compute phases tracking the ideal
    and the serial floors (flexure, pit-fill) flattening out."""
    cmap = plt.get_cmap("tab10")
    colours = {col: cmap(i) for i, (col, _) in enumerate(PHASES)}

    fig, axes = plt.subplots(1, len(files), figsize=(5.6 * len(files), 5.0),
                             squeeze=False)
    axes = axes[0]
    for ax, f in zip(axes, files):
        km = int("".join(c for c in os.path.basename(f) if c.isdigit()))
        p, series = _load_phases(f)
        ticks = sorted(set(p))
        # ideal 1/p guide anchored at the largest phase value at the baseline.
        anchor = max(v[0] for v in series.values())
        ax.plot([p[0], p[-1]], [anchor, anchor * p[0] / p[-1]], "--",
                color="0.55", lw=1.0, zorder=0, label="ideal $\\propto 1/p$")
        for col, lbl in PHASES:
            if col in series:
                ax.plot(p, series[col], "-o", color=colours[col], label=lbl,
                        lw=1.6, ms=5, mec="white", mew=0.5)
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        lticks = _label_ticks(ticks)
        ax.set_xticks(lticks)
        ax.set_xticklabels([str(t) for t in lticks], rotation=45, ha="right")
        ax.set_xlabel("MPI ranks (cores)")
        ax.set_title(MESH.get(km, "%d km" % km))
    axes[0].set_ylabel("phase wall-clock per call (s)")
    axes[0].legend(frameon=True, framealpha=0.9, fontsize=9.5, ncol=1,
                   loc="lower left")
    fig.tight_layout()
    fig.savefig(OUT_PHASES, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT_PHASES)


def main():
    files = sorted(glob.glob(os.path.join(HERE, "scaling_*km.csv")),
                   key=lambda f: -int("".join(c for c in os.path.basename(f) if c.isdigit())))
    if not files:
        raise SystemExit("no scaling_<N>km.csv files in %s" % HERE)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2))
    all_p = set()
    baselines = set()
    max_eff = 1.0

    for f in files:
        km = int("".join(c for c in os.path.basename(f) if c.isdigit()))
        p, sp, eff = _load(f)
        all_p.update(p)
        baselines.add(p[0])
        max_eff = max(max_eff, max(eff))
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
    lticks = _label_ticks(ticks)
    ax1.set_xticks(lticks)
    ax1.set_xticklabels([str(t) for t in lticks], rotation=45, ha="right")
    yt = [t / min(baselines) for t in lticks]
    ax1.set_yticks(yt)
    ax1.set_yticklabels(["%g" % v for v in yt])
    ax1.set_xlabel("MPI ranks (cores)")
    ax1.set_ylabel("speedup  (relative to smallest run)")
    ax1.set_title("(a) Strong-scaling speedup")
    ax1.legend(frameon=True, framealpha=0.9)

    # --- parallel efficiency (semilog-x) ---
    ax2.axhline(1.0, ls="--", color="0.55", lw=1.0, zorder=0, label="ideal (100 %)")
    ax2.set_xscale("log", base=2)
    ax2.set_xticks(lticks)
    ax2.set_xticklabels([str(t) for t in lticks], rotation=45, ha="right")
    ax2.set_ylim(0.0, max(1.15, max_eff * 1.08))
    ax2.set_xlabel("MPI ranks (cores)")
    ax2.set_ylabel("parallel efficiency")
    ax2.set_title("(b) Parallel efficiency")
    ax2.legend(frameon=True, framealpha=0.9)

    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT)

    phase_figure(files)


if __name__ == "__main__":
    main()
