#!/usr/bin/env python
"""
Plot strong-scaling speedup and parallel efficiency from a goSPL scaling
summary CSV (the ``scaling_summary.csv`` written by ``analyze_scaling.py``,
or any CSV with the same columns).

No single-CPU run is required. Large global models are memory-infeasible at
low rank counts, so speedup/efficiency are normalised to the SMALLEST rank
count present in the sweep (the baseline ``p0``):

    speedup     S(p) = T(p0) / T(p)            (S(p0) = 1)
    efficiency  E(p) = S(p) * p0 / p           (E(p0) = 1, i.e. 100 %)
    ideal line  = p / p0                       (NOT p)

State "baseline = p0 ranks" in any caption -- this is the standard convention
when the 1-rank run cannot be measured.

If the CSV also carries ``speedup_abs_amdahl`` / ``efficiency_abs_amdahl``
(an Amdahl-law fit extrapolated to a virtual 1 rank), those are overlaid as a
dashed "absolute" curve. The Amdahl framing avoids the >100 % efficiency that
the relative-to-p0 curve can show when the baseline run is itself memory-
bandwidth bound -- but it is an EXTRAPOLATION, so label it as such.

Usage
-----
    python plot_scaling.py scaling_summary.csv [out.png]

Writes ``scaling.png`` (two panels: speedup, efficiency) by default.
"""

import csv
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def _load(path):
    rows = list(csv.DictReader(open(path)))
    if not rows:
        raise SystemExit("No rows in %s" % path)
    rows.sort(key=lambda r: int(r["nranks"]))
    return rows


def plot(rows, out_png):
    p = [int(r["nranks"]) for r in rows]
    sp = [float(r["speedup"]) for r in rows]
    eff = [float(r["efficiency"]) for r in rows]
    p0 = p[0]
    ideal = [pi / p0 for pi in p]

    has_abs = "speedup_abs_amdahl" in rows[0] and "efficiency_abs_amdahl" in rows[0]
    if has_abs:
        sp_abs = [float(r["speedup_abs_amdahl"]) for r in rows]
        eff_abs = [float(r["efficiency_abs_amdahl"]) for r in rows]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))

    # speedup (log-log; ideal is the p/p0 diagonal)
    a1.plot(p, ideal, "k--", lw=1, label="ideal (linear)")
    a1.plot(p, sp, "o-", color="tab:blue",
            label="measured (rel. to %d ranks)" % p0)
    if has_abs:
        a1.plot(p, sp_abs, "s--", color="tab:green", alpha=0.7,
                label="Amdahl absolute (vs 1)")
    a1.set_xscale("log", base=2)
    a1.set_yscale("log", base=2)
    a1.set_xticks(p)
    a1.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    a1.set_yticks(ideal)
    a1.get_yaxis().set_major_formatter(plt.ScalarFormatter())
    a1.set_xlabel("MPI ranks")
    a1.set_ylabel("speedup")
    a1.set_title("Strong-scaling speedup")
    a1.grid(True, which="both", alpha=0.3)
    a1.legend()

    # efficiency (1.0 = ideal; >1 = superlinear, baseline is bandwidth-bound)
    a2.axhline(1.0, color="k", ls="--", lw=1, label="ideal (100 %)")
    a2.plot(p, eff, "o-", color="tab:blue",
            label="measured (rel. to %d ranks)" % p0)
    if has_abs:
        a2.plot(p, eff_abs, "s--", color="tab:green", alpha=0.7,
                label="Amdahl absolute")
    a2.set_xscale("log", base=2)
    a2.set_xticks(p)
    a2.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    a2.set_ylim(0, max(1.3, max(eff) * 1.1))
    a2.set_xlabel("MPI ranks")
    a2.set_ylabel("parallel efficiency")
    a2.set_title("Parallel efficiency")
    a2.grid(True, alpha=0.3)
    a2.legend()

    fig.suptitle("goSPL strong scaling (baseline = %d ranks)" % p0)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    print("wrote", out_png)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Plot goSPL strong-scaling speedup + efficiency from a "
        "scaling summary CSV (baseline = smallest rank count)."
    )
    parser.add_argument("csv", help="scaling_summary.csv (or compatible)")
    parser.add_argument("out", nargs="?", default="scaling.png",
                        help="output PNG (default: scaling.png)")
    args = parser.parse_args(argv)
    plot(_load(args.csv), args.out)


if __name__ == "__main__":
    main()
