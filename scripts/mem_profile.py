#!/usr/bin/env python
"""
Per-rank memory (RSS) of a goSPL run, by setup stage and by step phase.

    mpirun -n N python scripts/mem_profile.py input.yml [--json out.json]

For every stage (each mixin ``__init__`` in ``Model.__init__`` order, the mesh
build sub-steps, and the main step phases) it reports the PEAK resident memory
seen while the stage ran, on rank 0 and as the maximum over the other ranks.
Rank 0 is reported apart because it builds the serial mesh structure and holds
the rank-0-only tables, so it is usually the heaviest rank.

How to read it
--------------
* Run the same input at two rank counts. A stage whose peak on the OTHER ranks
  does not shrink as ranks increase holds something sized by the whole mesh on
  every rank (AGENTS.md > No whole-mesh arrays on every rank).
* RSS rarely drops after memory is freed, so the first stage that reaches a
  level "owns" it; later stages show the same value. Use the per-stage growth,
  not the absolute value, to attribute memory.
* The import baseline (PETSc, MPI, scipy, VTK) is ~0.55 GB per rank and does
  not depend on the mesh.

Sampling uses ``psutil`` when installed, else ``/proc/self/statm`` (Linux);
without either only the final peak (``ru_maxrss``) is available.
"""

from __future__ import annotations

import argparse
import functools
import json
import os
import resource
import sys
import threading
import time


def _rss_reader():
    try:
        import psutil

        proc = psutil.Process()
        return lambda: proc.memory_info().rss
    except ImportError:
        pass
    if os.path.exists("/proc/self/statm"):
        page = os.sysconf("SC_PAGE_SIZE")

        def statm():
            with open("/proc/self/statm") as f:
                return int(f.read().split()[1]) * page

        return statm
    scale = 1 if sys.platform == "darwin" else 1024     # bytes on macOS, kB on Linux
    return lambda: resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale


# Methods timed when present (setup sub-steps and the main step phases).
_PHASES = [
    "_meshfrom_cell_list", "_meshStructure", "readStratLayers",
    "applyForces", "getTectonics", "flowAccumulation", "erodepSPL",
    "sedChange", "seaChange", "getHillslope", "applyFlexure", "outputMesh",
]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("input", help="goSPL input YAML")
    ap.add_argument("--json", default=None, help="write the per-rank table here")
    args = ap.parse_args(argv)

    import gospl  # noqa: F401  (PETSc initialisation)
    import gospl.model as gm
    from mpi4py import MPI

    rss = _rss_reader()
    stack, peaks, order = [], {}, []
    lock = threading.Lock()

    def sampler():
        while True:
            r = rss()
            with lock:
                for lab in stack:
                    if r > peaks[lab]:
                        peaks[lab] = r
            time.sleep(0.005)

    def wrap(cls, name, label):
        f = cls.__dict__[name]

        @functools.wraps(f)
        def g(*a, **k):
            r0 = rss()
            with lock:
                if label not in peaks:
                    peaks[label] = r0
                    order.append(label)
                stack.append(label)
            try:
                return f(*a, **k)
            finally:
                r1 = rss()
                with lock:
                    stack.pop()
                    peaks[label] = max(peaks[label], r1)

        setattr(cls, name, g)

    Model = gm.Model
    for base in Model.__mro__[1:]:
        if base is not object and "__init__" in base.__dict__:
            wrap(base, "__init__", "init:" + base.__name__)
    for name in _PHASES:
        for base in Model.__mro__:
            if name in base.__dict__:
                wrap(base, name, name)
                break
    wrap(Model, "runProcesses", "runProcesses")

    threading.Thread(target=sampler, daemon=True).start()
    base_rss = rss()
    model = Model(args.input, verbose=False)
    model.runProcesses()
    sizes = (model.mpoints, model.lpoints)
    model.destroy()

    comm = MPI.COMM_WORLD
    rows = comm.gather({"base": base_rss, "peaks": dict(peaks), "order": list(order),
                        "mpoints": sizes[0], "lpoints": sizes[1]}, root=0)
    if comm.rank != 0:
        return 0

    MB = 2.0 ** 20
    labels = rows[0]["order"]
    print(f"np={comm.size}  mpoints={rows[0]['mpoints']:,}  "
          f"lpoints (rank 0)={rows[0]['lpoints']:,}  import baseline {rows[0]['base'] / MB:.0f} MB")
    print(f"{'stage':28s} {'rank 0 peak':>12s} {'others max':>11s} {'others mean':>12s}  (MB)")
    table = []
    for lab in labels:
        p0 = rows[0]["peaks"].get(lab, 0)
        po = [r["peaks"][lab] for r in rows[1:] if lab in r["peaks"]]
        omax = max(po) if po else float("nan")
        omean = sum(po) / len(po) if po else float("nan")
        table.append({"stage": lab, "rank0": p0, "others_max": omax, "others_mean": omean})
        print(f"{lab:28s} {p0 / MB:12.0f} {omax / MB:11.0f} {omean / MB:12.0f}")
    if args.json:
        with open(args.json, "w") as f:
            json.dump({"np": comm.size, "mpoints": rows[0]["mpoints"], "stages": table}, f, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
