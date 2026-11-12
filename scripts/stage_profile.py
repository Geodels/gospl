#!/usr/bin/env python
"""
Wall-clock time per sub-stage of a goSPL run, per time step, rank 0 and the
maximum over ranks.

    mpirun -n N python scripts/stage_profile.py input.yml

Use a short window (two or three steps) of the run you want to explain. It
times the marine routing pieces (`_matOcean` and inside it the smoother,
the serial `epsfill`, the whole-mesh gather), the marine cascade and coupled
deposition, the marine diffusion, the sediment-load and flow routing solves,
and the SPL eroder. The `--profile` phase profiler gives the coarse split
(flow / erosion / sed / sea); this script goes one level down, which is what
a stage that costs ~1000 s at scale but seconds on a smaller mesh needs.

Indented names run inside the one above them. Times are wall clock
(``time.perf_counter``), not the CPU time the verbose log prints.
"""

from __future__ import annotations

import functools
import sys
import time

# (module, class, method, label); module-level functions use class None.
_TARGETS = [
    ("gospl.sed.seaplex", "SEAMesh", "seaChange", "seaChange (marine stage)"),
    ("gospl.sed.seaplex", "SEAMesh", "_distanceCoasts", "  _distanceCoasts"),
    ("gospl.sed.seaplex", "SEAMesh", "_matOcean", "  _matOcean"),
    ("gospl.sed.hillslope", "hillSLP", "_solveSmooth", "    marine smoother solve"),
    ("gospl.sed.seaplex", None, "epsfill", "    epsfill (rank 0 serial)"),
    ("gospl.mesher.unstructuredmesh", "UnstMesh", "_gatherGlobalOnRoot", "    whole-mesh gather to rank 0"),
    ("gospl.sed.seaplex", "SEAMesh", "_distOcean", "  _distOcean (marine cascade)"),
    ("gospl.sed.seaplex", "SEAMesh", "_depMarineSystem", "  _depMarineSystem (Gmar)"),
    ("gospl.sed.hillslope", "hillSLP", "_diffuseOcean", "  _diffuseOcean"),
    ("gospl.sed.sedplex", "SEDMesh", "sedChange", "sedChange (continental)"),
    ("gospl.sed.sedplex", "SEDMesh", "_getSedFlux", "  _getSedFlux (sediment load)"),
    ("gospl.flow.flowplex", "FAMesh", "flowAccumulation", "flowAccumulation"),
    ("gospl.eroder.SPL", "SPL", "erodepSPL", "erodepSPL"),
    ("gospl.sed.hillslope", "hillSLP", "getHillslope", "getHillslope"),
]


def main(argv=None):
    import importlib

    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        print(__doc__)
        return 2

    import gospl  # noqa: F401  (PETSc initialisation)
    from mpi4py import MPI
    import gospl.model as gm

    T = {}

    def timed(label, f):
        @functools.wraps(f)
        def g(*a, **k):
            t = time.perf_counter()
            try:
                return f(*a, **k)
            finally:
                T[label] = T.get(label, 0.0) + time.perf_counter() - t
        return g

    order = []
    for modname, clsname, name, label in _TARGETS:
        mod = importlib.import_module(modname)
        if clsname is None:
            if hasattr(mod, name):
                setattr(mod, name, timed(label, getattr(mod, name)))
                order.append(label)
            continue
        cls = getattr(mod, clsname, None)
        if cls is not None and name in cls.__dict__:
            setattr(cls, name, timed(label, cls.__dict__[name]))
            order.append(label)

    model = gm.Model(argv[0], verbose=False)
    T.clear()
    t0 = time.perf_counter()
    model.runProcesses()
    total = time.perf_counter() - t0
    steps = max(1, round((model.tEnd - model.tStart) / model.dt))
    allT = MPI.COMM_WORLD.gather(T, root=0)
    size, mpoints = MPI.COMM_WORLD.size, model.mpoints
    model.destroy()

    if MPI.COMM_WORLD.rank == 0:
        print(f"np={size}  mpoints={mpoints:,}  steps={steps}  wall per step {total / steps:.2f} s")
        print(f"{'stage':38s} {'rank 0 (s/step)':>16s} {'max rank (s/step)':>18s}")
        for label in order:
            if any(label in t for t in allT):
                r0 = T.get(label, 0.0) / steps
                mx = max(t.get(label, 0.0) for t in allT) / steps
                print(f"{label:38s} {r0:16.2f} {mx:18.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
