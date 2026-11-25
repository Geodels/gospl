#!/usr/bin/env python
"""
Cost of goSPL's options: wall time per step and memory, each option added to
a base run, as a percentage of the base.

    python scripts/option_cost.py base.yml [-n 4] [--steps 8] [--only a,b]
                                  [--json out.json] [--workdir DIR]
    python scripts/option_cost.py --report DIR     # table from a finished run
    python scripts/option_cost.py --plot DIR1 [DIR2] --out fig.png   # bar chart
    python scripts/option_cost.py --table DIR1 DIR4 --out table.rst   # docs table

The base input should be a regional run with land, rivers and a sea (the
numbers in docs/user_guide/cost.rst come from the mixed_carbonate_clastic
margin of goSPL-examples, ~47,000 nodes). Each variant edits a copy of the
base (blocks added, keys changed), runs ``--steps`` steps with ``mpirun -n``
and records:

* the median wall time of a step, from the ``--summary`` record of every step
  after the first (the first one sets up solvers and caches; the median keeps
  one slow step, a solver rescue or a busy machine, from dominating);
* the time of the option's own phases (``--summary`` phases: ``flexure``,
  ``groundwater``, ``carbonate``, ...), which does not depend on how the
  option changes the landscape the other processes work on;
* the set-up time (``Model`` construction);
* the peak resident memory, on rank 0 and the largest on the other ranks
  (``ru_maxrss``).

Relative values are what to read: absolute times depend on the machine, the
mesh and the rank count. Paths in the base input are kept relative to its
folder (the runs start there); outputs go to ``--workdir``.
"""

from __future__ import annotations

import argparse
import copy
import json
import os
import subprocess
import sys
import tempfile

import yaml

_CHILD = r'''
import json, resource, sys, time
from mpi4py import MPI
from gospl.model import Model
yml, summ, out = sys.argv[1:4]
comm = MPI.COMM_WORLD
t0 = time.time()
m = Model(yml, verbose=False, showlog=False, summary=summ)
t1 = time.time()
try:
    m.runProcesses()
    npts = m.mpoints
finally:
    m.destroy()
scale = 1 if sys.platform == "darwin" else 1024
rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale
allr = comm.gather(rss, root=0)
if comm.rank == 0:
    json.dump({"init_s": t1 - t0, "rss0": allr[0],
               "rss_other": max(allr[1:]) if len(allr) > 1 else allr[0],
               "mpoints": int(npts), "ranks": comm.size}, open(out, "w"))
'''

# --------------------------------------------------------------------------
# Option blocks. Each variant is (name, label, function(cfg) -> None).
# --------------------------------------------------------------------------
STRATA = {"dual": True, "coarse": {"phi0": 0.49, "z0": 3700.0},
          "fine": {"phi0": 0.63, "z0": 1960.0, "k_factor": 1.5},
          "bedrock_coarse_frac": 0.6, "fine_diff_factor": 2.0}
CARB = {"insolation": 400.0, "factories": [
    {"name": "euphotic", "growth": 500.0, "extinction": 0.2, "saturation": 60.0,
     "framework": 0.5, "sand": 0.7,
     "response": [{"control": "depth", "x": [0.0, 20.0, 30.0], "f": [1.0, 1.0, 0.0]}]},
    {"name": "oligophotic", "growth": 300.0, "extinction": 0.05, "saturation": 60.0,
     "framework": 0.0, "sand": 0.5,
     "response": [{"control": "depth", "x": [0.0, 50.0, 80.0], "f": [1.0, 1.0, 0.0]}]}]}
GW = {"Ksat": 3.65e4, "specific_yield": 0.1, "aquifer_base": 50.0,
      "infiltration": 0.02, "conserve_baseflow": True}
DURI = {"form_rate": 1.0e-4, "max_thickness": 5.0, "fringe_depth": 3.0,
        "fringe_width": 2.0, "armor_max": 0.9}
GEOCHEM = {"conserve": True, "river_load": True, "marine_coupling": True,
           "species": [
               {"name": "carbonate", "weatherability": 1.0, "c_sat": 1.0,
                "precip_rate": 1.0, "solid_volume": 1.0},
               {"name": "silica", "weatherability": 0.5, "c_sat": 2.0,
                "precip_rate": 1.0, "solid_volume": 1.0}]}
SOIL = {"mode": "regolith", "soilK": 8.0e-6, "maxProd": 50.0e-6, "depthProd": 0.5,
        "uniform": 1.0}
WAVES = {"climate": [{"from": 90.0, "weight": 1.0}], "decay": 1.0e6,
         "shoal_depth": 20.0, "shoal_length": 5.0e3}


def _strat(c):
    c["time"]["strat"] = float(c["time"]["dt"])


def _strata(c):
    _strat(c)
    c["strata"] = copy.deepcopy(STRATA)


def _carb(c):
    _strata(c)
    c["carbonate"] = copy.deepcopy(CARB)


def _gw(c, duri=False, geo=False):
    g = copy.deepcopy(GW)
    if duri or geo:
        g["duricrust"] = copy.deepcopy(DURI)
    if geo:
        g["geochem"] = copy.deepcopy(GEOCHEM)
    c["groundwater"] = g


VARIANTS = [
    ("base", "base: linear SPL with deposition, linear hillslope, marine deposition",
     lambda c: None),
    ("nlspl", "non-linear SPL (n = 1.5)", lambda c: c["spl"].update({"n": 1.5})),
    ("spl_G", "SPL with transport-limited deposition (G = 0.5)",
     lambda c: c["spl"].update({"G": 0.5})),
    ("hill_nl", "non-linear hillslope (critical slope 0.6)",
     lambda c: c["diffusion"].update({"hillslopeSc": 0.6})),
    ("soil", "soil / regolith production (soil-aware SPL)",
     lambda c: c.update({"soil": copy.deepcopy(SOIL)})),
    ("flexure", "flexural isostasy (flat FEM)",
     lambda c: c.update({"flexure": {"method": "fem", "thick": 30.0e3, "rhoc": 2300.0,
                                     "rhoa": 3300.0, "young": 65e9, "nu": 0.25,
                                     "bcN": "Mirror", "bcE": "0Slope0Shear",
                                     "bcS": "Mirror", "bcW": "0Slope0Shear"}})),
    ("orography", "orographic rain",
     lambda c: c.update({"orography": {"wind_speed": 15.0, "wind_dir": 270.0,
                                       "conv_time": 200.0, "fall_time": 200.0,
                                       "precip_base": 0.5, "precip_min": 0.01,
                                       "rainfall_frequency": 1.0}})),
    ("ice", "diagnostic glacial model",
     lambda c: c.update({"ice": {"hela": 150.0, "hice": 260.0, "hterm": 100.0}})),
    ("strata", "stratigraphy (one layer per step)", _strat),
    ("dual", "stratigraphy + dual lithology (sand / mud)", _strata),
    ("prov", "stratigraphy + dual + provenance tracers (3 classes)", lambda c: (
        _strata(c), c.update({"provenance": {"classes": 3, "uniform": 0}}))),
    ("gw", "groundwater (water table)", lambda c: _gw(c)),
    ("gw_duri", "groundwater + duricrust", lambda c: _gw(c, duri=True)),
    ("gw_geochem", "groundwater + duricrust + geochemistry (2 species, river load)",
     lambda c: _gw(c, geo=True)),
    ("gravity", "stratigraphy + dual + gravity currents", lambda c: (
        _strata(c), c.update({"gravity": {"fraction": 0.3, "coarse_factor": 0.3,
                                          "runout": 1.0e4, "runout_fine": 2.0e4}}))),
    ("carbonate", "stratigraphy + dual + carbonate production (2 factories)", _carb),
    ("carb_exposure", "carbonate + open-water index (exposure)", lambda c: (
        _carb(c), c["carbonate"].update({"exposure_length": 3.0e3}),
        c["carbonate"]["factories"][0]["response"].append(
            {"control": "exposure", "x": [0.0, 3.0], "f": [0.2, 1.0]}))),
    ("carb_turbidity", "carbonate + turbidity control", lambda c: (
        _carb(c), c["carbonate"]["factories"][0]["response"].append(
            {"control": "turbidity", "length": 1.0e4, "x": [6.0e-3, 1.2e-2],
             "f": [1.0, 0.0]}))),
    ("pelagic", "stratigraphy + dual + pelagic rain", lambda c: (
        _strata(c), c.update({"pelagic": [{"start": float(c["time"]["start"]),
                                           "carbonate": 5.0, "lysocline": 300.0,
                                           "ccd": 450.0, "clay": 1.0,
                                           "clay_decay": 3.0e4}]}))),
    ("waves", "carbonate + waves (wave control)", lambda c: (
        _carb(c), c.update({"waves": copy.deepcopy(WAVES)}),
        c["carbonate"]["factories"][0]["response"].append(
            {"control": "wave", "x": [0.0, 0.5], "f": [0.2, 1.0]}))),
    ("coastal", "carbonate + waves + coastal erosion", lambda c: (
        _carb(c), c.update({"waves": copy.deepcopy(WAVES),
                            "coastal_erosion": {"rate": 1.0e-3, "depth_scale": 5.0,
                                                "cliff_height": 5.0, "max_erosion": 2.0,
                                                "rubble": "carb_sand"}}))),
    ("karst", "carbonate + karst dissolution", lambda c: (
        _carb(c), c.update({"karst": {"rate": 5.0e-5, "focus": 0.3}}))),
]


# Phases (``--summary``) that belong to an option: their time is the option's
# own cost. Options without one (erosion laws, lithologies, gravity currents
# inside the marine step) are read from the change of the whole step.
OWN = {"nlspl": ["erosion"], "spl_G": ["erosion"], "hill_nl": ["hillslope"],
       "soil": ["erosion", "hillslope"], "flexure": ["flexure"], "orography": ["forcing"],
       "ice": ["ice", "till"], "strata": ["strat"], "gw": ["groundwater"],
       "gw_duri": ["groundwater"], "gw_geochem": ["groundwater"],
       "carbonate": ["carbonate"], "carb_exposure": ["carbonate"],
       "carb_turbidity": ["carbonate"], "pelagic": ["pelagic"],
       "waves": ["carbonate"], "coastal": ["carbonate", "coastal erosion"],
       "karst": ["karst"]}


def _stats(summ):
    """Median step time and median time per phase over the steps after the
    first."""
    import statistics

    rows = [json.loads(l) for l in open(summ) if '"wall_s"' in l and '"end"' not in l]
    rows = rows[1:] if len(rows) > 1 else rows
    keys = set().union(*[r.get("phases", {}).keys() for r in rows])
    ph = {k: statistics.median(r.get("phases", {}).get(k, 0.0) for r in rows) for k in keys}
    return statistics.median(r["wall_s"] for r in rows), ph


def report(work, out=None):
    """Table of a finished run (``res_*.json`` and ``sum_*.jsonl`` in
    ``work``), in the order of VARIANTS, relative to the base."""
    res = []
    for name, label, _ in VARIANTS:
        f = os.path.join(work, "res_%s.json" % name)
        if not os.path.exists(f):
            continue
        d = json.load(open(f))
        d["step_s"], d["phases"] = _stats(os.path.join(work, "sum_%s.jsonl" % name))
        d.update(name=name, label=label)
        res.append(d)
    b = res[0]
    print("%-15s %8s %6s %8s %7s %9s %9s" % ("option", "step s", "step%", "own s",
                                             "own%", "mem0%", "memOth%"))
    for d in res:
        own = sum(d["phases"].get(k, 0.0) for k in OWN.get(d["name"], []))
        if d["name"] in ("hill_nl", "soil", "nlspl", "spl_G"):
            own -= sum(b["phases"].get(k, 0.0) for k in OWN[d["name"]])
        d["own_s"] = max(own, 0.0)
        print("%-15s %8.3f %6.0f %8.3f %7.0f %+9.0f %+9.0f" % (
            d["name"], d["step_s"], 100 * d["step_s"] / b["step_s"], d["own_s"],
            100 * d["own_s"] / b["step_s"], 100 * (d["rss0"] / b["rss0"] - 1),
            100 * (d["rss_other"] / b["rss_other"] - 1)))
    if out:
        json.dump(res, open(out, "w"), indent=1)
    return res


SHORT = {"base": "base", "nlspl": "non-linear SPL", "spl_G": "transport-limited SPL",
         "hill_nl": "non-linear hillslope", "soil": "soil", "flexure": "flexure",
         "orography": "orography", "ice": "ice", "strata": "stratigraphy",
         "dual": "dual lithology", "prov": "provenance (3 classes)", "gw": "groundwater",
         "gw_duri": "+ duricrust", "gw_geochem": "+ geochemistry",
         "gravity": "gravity currents", "carbonate": "carbonates",
         "carb_exposure": "+ exposure", "carb_turbidity": "+ turbidity",
         "pelagic": "pelagic rain", "waves": "+ waves", "coastal": "+ coastal erosion",
         "karst": "+ karst"}


# What each option is added to, and a note for the documentation table.
ADDED = {"dual": "stratigraphy", "prov": "dual", "gw_duri": "groundwater",
         "gw_geochem": "duricrust", "gravity": "dual", "pelagic": "dual",
         "carbonate": "dual", "carb_exposure": "carbonates",
         "carb_turbidity": "carbonates", "waves": "carbonates",
         "coastal": "carbonates", "karst": "carbonates"}
NOTES = {"nlspl": "changes the relief: the marine step does less work here",
         "soil": "soil-aware SPL and soil diffusion (SNES, TS)",
         "flexure": "operator factorised once",
         "orography": "lowers the rain here: less flow and marine work",
         "strata": "memory grows with the layers (below)",
         "dual": "a second (mud) routing in the cascades",
         "prov": "one routed column per class",
         "gw": "AMG head solve",
         "gw_geochem": "one direct (MUMPS) solve per species",
         "gravity": "its own network and solve, inside the marine step",
         "pelagic": "clay decay solve, and carb_mud added to the marine routing",
         "carbonate": "three more lithologies in the marine routing",
         "carb_exposure": "one screened-Poisson solve",
         "waves": "one upwind solve per wave direction"}
GROUPS = [("Erosion and hillslope laws", ["nlspl", "spl_G", "hill_nl", "soil"]),
          ("Climate, ice and the lithosphere", ["flexure", "orography", "ice"]),
          ("Stratigraphy", ["strata", "dual", "prov"]),
          ("Groundwater and geochemistry", ["gw", "gw_duri", "gw_geochem"]),
          ("Marine processes", ["gravity", "pelagic", "carbonate", "carb_exposure",
                                "carb_turbidity", "waves", "coastal", "karst"])]


def table(dirs, out):
    """The documentation table (RST list-table) of two finished runs, one
    rank and several ranks: % of the base step and the option's own phases
    (on the second run)."""
    import io
    import contextlib

    runs = []
    for d in dirs:
        with contextlib.redirect_stdout(io.StringIO()):
            runs.append({x["name"]: x for x in report(d)})
    r1, rn = runs[0], runs[-1]
    pct = lambda r, n, k="step_s": "%.0f" % (100 * r[n][k] / r["base"]["step_s"]) \
        if n in r else ""
    L = [".. list-table::", "   :header-rows: 1", "   :widths: 28 13 9 9 9 32", "",
         "   * - Option", "     - added to", "     - %d rank" % r1["base"]["ranks"],
         "     - %d ranks" % rn["base"]["ranks"], "     - own phases", "     - notes",
         "   * - base", "     -", "     - 100 (%.2f s)" % r1["base"]["step_s"],
         "     - 100 (%.2f s)" % rn["base"]["step_s"], "     -",
         "     - mostly the marine deposition and the flow accumulation"]
    for title, names in GROUPS:
        L += ["   * - **%s**" % title] + ["     -"] * 5
        for n in names:
            if n not in rn:
                continue
            L += ["   * - %s" % SHORT[n], "     - %s" % ADDED.get(n, "base"),
                  "     - %s" % pct(r1, n), "     - %s" % pct(rn, n),
                  "     - %s" % pct(rn, n, "own_s"), "     - %s" % NOTES.get(n, "")]
    open(out, "w").write("\n".join(L) + "\n")


def plot(dirs, out):
    """Bar chart of the step time of every option (% of the base), whole
    step and the option's own phases, one bar per finished run directory."""
    import io
    import contextlib

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    runs = []
    for d in dirs:
        with contextlib.redirect_stdout(io.StringIO()):
            r = report(d)
        runs.append(({x["name"]: x for x in r}, "%d rank%s" % (r[0]["ranks"],
                                                             "s" if r[0]["ranks"] > 1 else "")))
    names = [n for n, _, _ in VARIANTS if n in runs[0][0]]
    fig, ax = plt.subplots(figsize=(10, 0.32 * len(names) + 1.6), constrained_layout=True)
    y = np.arange(len(names))[::-1]
    h = 0.8 / len(runs)
    cols = ["#457b9d", "#e76f51", "#2a9d8f", "#8d99ae"]
    vmax = 0.0
    for k, (d, lab) in enumerate(runs):
        b = d["base"]["step_s"]
        tot = [100 * d[n]["step_s"] / b if n in d else 0 for n in names]
        own = [100 * d[n].get("own_s", 0.0) / b if n in d else 0 for n in names]
        off = (0.5 * (len(runs) - 1) - k) * h
        ax.barh(y + off, tot, h, color=cols[k], alpha=0.4, label="whole step, %s" % lab)
        ax.barh(y + off, own, h, color=cols[k], label="option's own phases, %s" % lab)
        vmax = max(vmax, max(tot))
    ax.axvline(100, color="k", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([SHORT.get(n, n) for n in names])
    ax.set_xlabel("%% of the base step (%d nodes)" % runs[0][0]["base"]["mpoints"])
    ax.set_title("Time per step: each option added to the base, (+) to the option above",
                 fontsize=10)
    ax.set_xlim(0, 1.05 * vmax)
    ax.grid(axis="x", alpha=0.25, lw=0.5)
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    fig.savefig(out, dpi=150)


def base_config(path, steps):
    cfg = yaml.safe_load(open(path))
    # keep only the core blocks: every option is switched on by a variant
    keep = ("name", "domain", "time", "spl", "diffusion", "sea", "climate", "output")
    cfg = {k: v for k, v in cfg.items() if k in keep}
    t = cfg["time"]
    t["end"] = float(t["start"]) + steps * float(t["dt"])
    t["tout"] = t["end"] - float(t["start"])
    t.pop("strat", None)
    cfg["spl"].pop("G", None)
    cfg["spl"].pop("n", None)
    return cfg


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("input", nargs="?", help="base goSPL input YAML")
    ap.add_argument("--report", default=None, metavar="DIR",
                    help="print the table of a finished run (its --workdir)")
    ap.add_argument("--plot", nargs="+", default=None, metavar="DIR",
                    help="bar chart of finished runs (one bar per DIR)")
    ap.add_argument("--table", nargs="+", default=None, metavar="DIR",
                    help="RST table of finished runs (one rank, then several)")
    ap.add_argument("--out", default="option_cost.png", help="output of --plot/--table")
    ap.add_argument("-n", type=int, default=4, help="MPI ranks (default 4)")
    ap.add_argument("--steps", type=int, default=8, help="steps per run (default 8)")
    ap.add_argument("--only", default=None, help="comma-separated variant names")
    ap.add_argument("--json", default=None, help="write the results here")
    ap.add_argument("--workdir", default=None, help="outputs (default: a temp dir)")
    a = ap.parse_args(argv)
    if a.report:
        report(a.report, a.json)
        return 0
    if a.plot:
        plot(a.plot, a.out)
        return 0
    if a.table:
        table(a.table, a.out)
        return 0

    base = os.path.abspath(a.input)
    here = os.path.dirname(base)
    work = os.path.abspath(a.workdir or tempfile.mkdtemp(prefix="gospl_cost_"))
    os.makedirs(work, exist_ok=True)
    child = os.path.join(work, "_child.py")
    open(child, "w").write(_CHILD)
    want = None if a.only is None else set(a.only.split(","))
    res = []
    for name, label, edit in VARIANTS:
        if want is not None and name not in want and name != "base":
            continue
        cfg = base_config(base, a.steps)
        edit(cfg)
        cfg["output"] = {"dir": os.path.join(work, "out_" + name), "makedir": False}
        yml = os.path.join(here, "_cost_%s.yml" % name)
        yaml.safe_dump(cfg, open(yml, "w"), sort_keys=False)
        summ = os.path.join(work, "sum_%s.jsonl" % name)
        out = os.path.join(work, "res_%s.json" % name)
        for f in (summ, out):
            if os.path.exists(f):
                os.remove(f)
        cmd = ["mpirun", "-n", str(a.n), sys.executable, child, yml, summ, out]
        r = subprocess.run(cmd, cwd=here, capture_output=True, text=True)
        os.remove(yml)
        if r.returncode != 0 or not os.path.exists(out):
            print("%-15s FAILED: %s" % (name, (r.stderr or r.stdout)[-400:]), flush=True)
            continue
        d = json.load(open(out))
        d["step_s"], d["phases"] = _stats(summ)
        d.update(name=name, label=label)
        res.append(d)
        b = res[0]
        print("%-15s step %7.3f s (%4.0f %%)  init %6.1f s  mem rank0 %6.0f MB (%+4.0f %%)"
              "  others %6.0f MB (%+4.0f %%)" % (
                  name, d["step_s"], 100 * d["step_s"] / b["step_s"], d["init_s"],
                  d["rss0"] / 2**20, 100 * (d["rss0"] / b["rss0"] - 1),
                  d["rss_other"] / 2**20, 100 * (d["rss_other"] / b["rss_other"] - 1)),
              flush=True)
    if a.json:
        json.dump(res, open(a.json, "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
