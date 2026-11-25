#!/usr/bin/env python
"""
How the cost of goSPL's options scales with the mesh size and the number of
MPI ranks, and the figures of docs/user_guide/cost.rst.

    python scripts/option_scaling.py base.yml --workdir DIR \
        [--spacings 1400,700,350,175] [--ranks 4] \
        [--strong-spacing 350 --strong-ranks 1,2,4,8,16] \
        [--options base,dual,prov,...] [--steps 6]
    python scripts/option_scaling.py --plot DIR --out docs/_static/cost

The base input must be a flat run whose ``domain: npdata`` mesh (keys ``v``,
``c``, ``z``) and ``climate: map`` share one npz file, as the
``mixed_carbonate_clastic`` margin of goSPL-examples. For each spacing the
domain is re-meshed with an equilateral lattice and the elevation and the
other per-node keys of that file are interpolated linearly from the original
mesh (no new topographic detail on finer meshes: the cost, not the landscape,
is measured). Every (spacing, ranks) pair runs the chosen variants of
``scripts/option_cost.py`` and keeps their results in ``DIR/<dx>_np<n>/``;
``--plot`` draws the figures from them.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))


def remesh(src, dst, dx):
    """Equilateral lattice of spacing ``dx`` over the bounding box of the
    mesh ``src`` (npz), every per-node key interpolated linearly."""
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
    from scipy.spatial import Delaunay

    d = np.load(src)
    v = d["v"]
    (x0, y0), (x1, y1) = v[:, :2].min(0), v[:, :2].max(0)
    dy = dx * np.sqrt(3.0) / 2.0
    ys = np.linspace(y0, y1, int(round((y1 - y0) / dy)) + 1)
    pts = []
    for j, y in enumerate(ys):
        xs = np.arange(x0, x1 + 0.5 * dx, dx)
        if j % 2:
            xs = np.concatenate([[x0], xs[:-1] + 0.5 * dx, [x1]])
        xs = np.clip(xs, x0, x1)
        pts.append(np.column_stack([xs, np.full(xs.size, y)]))
    xy = np.unique(np.vstack(pts), axis=0)
    cells = Delaunay(xy).simplices
    out = {"v": np.column_stack([xy, np.zeros(len(xy))]), "c": cells.astype(np.int64)}
    tri = Delaunay(v[:, :2])
    for k in d.files:
        if k in ("v", "c") or d[k].shape[0] != v.shape[0]:
            continue
        f = LinearNDInterpolator(tri, d[k])(xy)
        bad = ~np.isfinite(f).all(axis=-1) if f.ndim > 1 else ~np.isfinite(f)
        if bad.any():
            f[bad] = NearestNDInterpolator(v[:, :2], d[k])(xy[bad])
        out[k] = f
    np.savez(dst, **out)
    return len(xy)


def run(base, work, spacings, ranks_list, options, steps):
    cfg = yaml.safe_load(open(base))
    here = os.path.dirname(os.path.abspath(base))
    mesh = os.path.join(here, cfg["domain"]["npdata"][0] + ".npz")
    for dx in spacings:
        tag = "dx%d" % dx
        mname = "_scaling_%s" % tag
        npts = remesh(mesh, os.path.join(here, mname + ".npz"), dx)
        c = yaml.safe_load(open(base))
        c["domain"]["npdata"][0] = mname
        for cl in c.get("climate", []):
            if "map" in cl:
                cl["map"][0] = mname
        yml = os.path.join(here, "_scaling_%s.yml" % tag)
        yaml.safe_dump(c, open(yml, "w"), sort_keys=False)
        for n in ranks_list:
            out = os.path.join(work, "%s_np%d" % (tag, n))
            print("== dx %g m (%d nodes), %d ranks" % (dx, npts, n), flush=True)
            subprocess.run([sys.executable, os.path.join(HERE, "option_cost.py"), yml,
                            "-n", str(n), "--steps", str(steps), "--only", ",".join(options),
                            "--workdir", out], cwd=here, check=False)
        os.remove(yml)
        os.remove(os.path.join(here, mname + ".npz"))


def collect(work):
    sys.path.insert(0, HERE)
    import option_cost as oc

    rows = []
    for d in sorted(os.listdir(work)):
        if not d.startswith("dx") or "_np" not in d:
            continue
        dx, n = d[2:].split("_np")
        for name, label, _ in oc.VARIANTS:
            f = os.path.join(work, d, "res_%s.json" % name)
            if not os.path.exists(f):
                continue
            r = json.load(open(f))
            r["step_s"], r["phases"] = oc._stats(os.path.join(work, d, "sum_%s.jsonl" % name))
            r.update(name=name, label=label, dx=float(dx), ranks=int(n))
            rows.append(r)
    return rows


LABELS = {"base": "base", "dual": "+ stratigraphy, dual lithology",
          "prov": "+ provenance (3 classes)", "gravity": "+ gravity currents",
          "carbonate": "+ carbonates", "pelagic": "+ pelagic rain",
          "gw_geochem": "+ groundwater and geochemistry", "soil": "+ soil",
          "flexure": "+ flexure", "strata": "+ stratigraphy"}


def plot(work, out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = collect(work)
    json.dump(rows, open(os.path.join(work, "scaling.json"), "w"), indent=1)
    names = [n for n in LABELS if any(r["name"] == n for r in rows)]
    cmap = plt.get_cmap("tab10")
    col = {n: cmap(i % 10) for i, n in enumerate(names)}
    ranks = sorted({r["ranks"] for r in rows})
    size_n = max(ranks, key=lambda n: len({r["dx"] for r in rows if r["ranks"] == n}))
    sz = [r for r in rows if r["ranks"] == size_n]

    # 1. time per step against the mesh size
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    for n in names:
        rr = sorted([r for r in sz if r["name"] == n], key=lambda r: r["mpoints"])
        if not rr:
            continue
        ax[0].loglog([r["mpoints"] for r in rr], [r["step_s"] for r in rr], "o-",
                     color=col[n], label=LABELS[n], lw=1.6 if n == "base" else 1.1)
        base = {r["mpoints"]: r["step_s"] for r in sz if r["name"] == "base"}
        if n != "base":
            ax[1].semilogx([r["mpoints"] for r in rr if r["mpoints"] in base],
                           [100 * r["step_s"] / base[r["mpoints"]] for r in rr
                            if r["mpoints"] in base], "o-", color=col[n], label=LABELS[n])
    m = np.array(sorted({r["mpoints"] for r in sz}), dtype=float)
    b0 = min((r for r in sz if r["name"] == "base"), key=lambda r: r["mpoints"])
    ax[0].loglog(m, b0["step_s"] * m / b0["mpoints"], "k:", lw=0.8, label="linear in nodes")
    ax[0].set_xlabel("mesh nodes")
    ax[0].set_ylabel("wall time per step (s)")
    ax[0].set_title("time per step, %d ranks" % size_n, fontsize=10)
    ax[0].legend(fontsize=7, frameon=False)
    ax[1].axhline(100, color="k", lw=0.8)
    ax[1].set_xlabel("mesh nodes")
    ax[1].set_ylabel("% of the base step")
    ax[1].set_title("cost relative to the base run", fontsize=10)
    for a in ax:
        a.grid(alpha=0.25, lw=0.5, which="both")
    fig.savefig(out + "_size.png", dpi=150)

    # 2. strong scaling at one mesh size
    st_dx = max({r["dx"] for r in rows}, key=lambda dx: len({r["ranks"] for r in rows
                                                              if r["dx"] == dx}))
    st = [r for r in rows if r["dx"] == st_dx]
    if len({r["ranks"] for r in st}) > 1:
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
        for n in names:
            rr = sorted([r for r in st if r["name"] == n], key=lambda r: r["ranks"])
            if len(rr) < 2:
                continue
            p = [r["ranks"] for r in rr]
            t = [r["step_s"] for r in rr]
            ax[0].loglog(p, t, "o-", color=col[n], label=LABELS[n])
            ax[1].semilogx(p, [100 * t[0] * p[0] / (ti * pi) for ti, pi in zip(t, p)], "o-",
                           color=col[n], label=LABELS[n])
        p = np.array(sorted({r["ranks"] for r in st}), dtype=float)
        ax[0].set_xlabel("MPI ranks")
        ax[0].set_ylabel("wall time per step (s)")
        ax[0].set_title("strong scaling, %d nodes" % st[0]["mpoints"], fontsize=10)
        ax[0].legend(fontsize=7, frameon=False)
        ax[1].axhline(100, color="k", lw=0.8)
        ax[1].set_xlabel("MPI ranks")
        ax[1].set_ylabel("parallel efficiency (%)")
        ax[1].set_ylim(0, 130)
        from matplotlib.ticker import NullFormatter

        for a in ax:
            a.grid(alpha=0.25, lw=0.5, which="both")
            a.set_xticks(p)
            a.set_xticklabels([str(int(x)) for x in p])
            a.xaxis.set_minor_formatter(NullFormatter())
        fig.savefig(out + "_ranks.png", dpi=150)

    # 3. memory per rank against the nodes per rank
    fig, ax = plt.subplots(figsize=(5.6, 4.2), constrained_layout=True)
    for n in names[::-1]:
        rr = sorted([r for r in rows if r["name"] == n],
                    key=lambda r: r["mpoints"] / r["ranks"])
        if not rr:
            continue
        ax.loglog([r["mpoints"] / r["ranks"] for r in rr],
                  [max(r["rss0"], r["rss_other"]) / 2**20 for r in rr], "o",
                  color=col[n], label=LABELS[n], ms=7 if n == "base" else 4,
                  mfc="none" if n == "base" else col[n], zorder=3 if n == "base" else 2)
    ax.set_xlabel("mesh nodes per rank")
    ax.set_ylabel("peak memory of the largest rank (MB)")
    ax.set_title("memory per rank", fontsize=10)
    ax.grid(alpha=0.25, lw=0.5, which="both")
    ax.legend(fontsize=7, frameon=False)
    fig.savefig(out + "_memory.png", dpi=150)
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("input", nargs="?", help="base goSPL input YAML")
    ap.add_argument("--workdir", required=False)
    ap.add_argument("--spacings", default="1400,700,350,175")
    ap.add_argument("--ranks", default="4")
    ap.add_argument("--strong-spacing", type=float, default=None)
    ap.add_argument("--strong-ranks", default="1,2,4,8,16")
    ap.add_argument("--options", default="base,dual,prov,gravity,carbonate,pelagic,"
                                          "gw_geochem,soil,flexure")
    ap.add_argument("--steps", type=int, default=6)
    ap.add_argument("--plot", default=None, metavar="DIR")
    ap.add_argument("--out", default="option_scaling")
    a = ap.parse_args(argv)
    if a.plot:
        plot(a.plot, a.out)
        return 0
    opts = a.options.split(",")
    work = os.path.abspath(a.workdir)
    os.makedirs(work, exist_ok=True)
    run(a.input, work, [float(x) for x in a.spacings.split(",")],
        [int(x) for x in a.ranks.split(",")], opts, a.steps)
    if a.strong_spacing:
        run(a.input, work, [a.strong_spacing],
            [int(x) for x in a.strong_ranks.split(",") if int(x) not in
             [int(y) for y in a.ranks.split(",")]], opts, a.steps)
    return 0


if __name__ == "__main__":
    sys.exit(main())
