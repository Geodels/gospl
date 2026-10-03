"""
Machine-readable per-step run summary (JSON Lines).

Why
---
Diagnosing a goSPL run usually starts from a verbose log: free-form prints from
a dozen modules, partly rank-0 only, hard to compare between runs and costly to
read. With ``Model(..., summary="run.jsonl")`` (or ``gospl -i in.yml --summary
run.jsonl``) rank 0 appends ONE JSON object per time step to that file, so a
run can be checked, plotted or diffed by a script (``gospl-inspect --summary``)
and pasted from an HPC system as a single blob.

What a record holds (all values global, i.e. identical whatever the rank count
up to solver round-off):

* ``t``, ``dt``, ``step``, ``wall_s`` (rank-0 wall time of the step) and
  ``phases`` (rank-0 wall seconds per profiler phase during the step);
* ``elev`` min / max / mean, ``FA_max``, ``sealevel``;
* ``nonfinite``: count of NaN/inf over owned nodes in elevation, cumulative
  erosion-deposition and discharge -- the first number to check;
* ``volume``: this step's eroded (<= 0) and deposited (>= 0) volume and their
  sum ``net`` (m^3; from the change in ``cumED`` times cell area). On a closed
  domain ``net`` is ~0; on an open one it is minus the boundary outflux;
* ``ksp``: per-class flow-KSP counters (``fatal`` = main discharge solve,
  ``cascade`` = downstream cascade solves, ``other`` = every other caller of
  ``_solve_KSP``): number of solves, total and max iterations, primary
  failures;
* ``events``: notable solver outcomes (fallback failures with the un-drained
  cell count and whether they ponded or aborted, cascade stop reasons and pass
  counts, SNES fallbacks and soil sub-stepping).

A final ``{"end": true, ...}`` record carries the total wall time.

Cost and MPI
------------
Off by default. When off, every hook is a single ``getattr`` returning
``None``. When on, :meth:`RunSummary.step_end` performs a handful of
collectives (Vec min/max/sum, one ``allreduce``) and MUST be called on every
rank -- ``Model.runProcesses`` does so unconditionally. The solver hooks only
append rank-identical values (KSP reasons and iteration counts are global) to
a Python list, so they communicate nothing.
"""

from __future__ import annotations

import json

import numpy as np
from mpi4py import MPI

MPIrank = MPI.COMM_WORLD.Get_rank()


# ---------------------------------------------------------------------------
# Hooks called from the process modules (no-ops unless a summary is active)
# ---------------------------------------------------------------------------


def record(model, kind, **data):
    """Append a solver event to the current step (no-op when summary is off)."""
    ev = getattr(model, "_summaryEvents", None)
    if ev is None:
        return
    item = {"kind": kind}
    for k, v in data.items():
        if isinstance(v, (np.integer,)):
            v = int(v)
        elif isinstance(v, (np.floating,)):
            v = float(v)
        item[k] = v
    ev.append(item)


def count_ksp(model, key, its, reason):
    """Accumulate one KSP solve into the per-step counters."""
    stats = getattr(model, "_summaryKSP", None)
    if stats is None:
        return
    s = stats.setdefault(key, {"solves": 0, "its": 0, "max_its": 0, "fails": 0})
    s["solves"] += 1
    s["its"] += int(its)
    s["max_its"] = max(s["max_its"], int(its))
    if reason < 0:
        s["fails"] += 1


# ---------------------------------------------------------------------------
# Per-step writer
# ---------------------------------------------------------------------------


class RunSummary(object):
    """Write one JSON line per time step on rank 0.

    :arg model: the goSPL ``Model``
    :arg path: output file (truncated on a fresh run, appended on a restart)
    """

    def __init__(self, model, path):
        self.path = str(path)
        self.step = 0
        self._t0 = None
        self._prof0 = {}
        self._cum_prev = None
        self._run_t0 = MPI.Wtime()
        model._summaryEvents = []
        model._summaryKSP = {}
        if MPIrank == 0:
            mode = "a" if getattr(model, "rStep", 0) > 0 else "w"
            with open(self.path, mode) as f:
                from gospl import __version__

                f.write(json.dumps({
                    "start": True,
                    "version": __version__,
                    "input": getattr(model, "finput", None),
                    "nranks": MPI.COMM_WORLD.Get_size(),
                    "mpoints": int(getattr(model, "mpoints", 0)),
                    "flat": bool(getattr(model, "flatModel", False)),
                    # closed = no drainable outlet (a sphere, an all-wall box):
                    # only then must the per-step volume budget close to ~0.
                    "closed": not bool(getattr(model, "_domainHasOutlet", True)),
                    "tStart": float(getattr(model, "tStart", 0.0)),
                    "tEnd": float(getattr(model, "tEnd", 0.0)),
                    "dt": float(getattr(model, "dt", 0.0)),
                }) + "\n")

    # -- helpers ---------------------------------------------------------
    @staticmethod
    def _owned(model):
        return model.inIDs == 1

    def _cum_local(self, model):
        return model.cumEDLocal.getArray().copy()

    # -- per step --------------------------------------------------------
    def step_begin(self, model):
        """Snapshot timers and the cumulative deposit (rank-local, no comm)."""
        self._t0 = MPI.Wtime()
        self._prof0 = dict(getattr(model.profiler, "times", {}))
        if self._cum_prev is None:
            self._cum_prev = self._cum_local(model)
        model._summaryEvents.clear()
        model._summaryKSP.clear()

    def step_end(self, model):
        """Reduce the step's diagnostics and write one record. COLLECTIVE."""
        owned = self._owned(model)
        area = np.asarray(model.larea)

        # Elevation / discharge extremes (collective Vec reductions).
        hmin = float(model.hGlobal.min()[1])
        hmax = float(model.hGlobal.max()[1])
        hmean = float(model.hGlobal.sum()) / max(int(model.mpoints), 1)
        famax = float(model.FAG.max()[1])

        # Non-finite counts + this step's eroded / deposited volume.
        h = model.hLocal.getArray()
        cum = self._cum_local(model)
        fa = model.FAL.getArray()
        d = cum - self._cum_prev
        self._cum_prev = cum
        dv = np.where(owned, d * area, 0.0)
        dv = np.where(np.isfinite(dv), dv, 0.0)
        local = np.array([
            np.count_nonzero(~np.isfinite(h[owned])),
            np.count_nonzero(~np.isfinite(cum[owned])),
            np.count_nonzero(~np.isfinite(fa[owned])),
            dv[dv < 0.0].sum(),
            dv[dv > 0.0].sum(),
        ], dtype=np.float64)
        glob = np.zeros_like(local)
        MPI.COMM_WORLD.Allreduce(local, glob, op=MPI.SUM)

        if MPIrank != 0:
            return

        now = MPI.Wtime()
        prof = getattr(model.profiler, "times", {})
        phases = {k: round(v - self._prof0.get(k, 0.0), 4)
                  for k, v in prof.items() if v - self._prof0.get(k, 0.0) > 0.0}
        rec = {
            "step": self.step,
            "t": float(model.tNow),
            "dt": float(model.dt),
            "wall_s": round(now - self._t0, 4),
            "phases": phases,
            "elev": {"min": hmin, "max": hmax, "mean": hmean},
            "FA_max": famax,
            "sealevel": float(getattr(model, "sealevel", 0.0)),
            "nonfinite": {"elev": int(glob[0]), "cumED": int(glob[1]),
                          "FA": int(glob[2])},
            "volume": {"eroded": float(glob[3]), "deposited": float(glob[4]),
                       "net": float(glob[3] + glob[4])},
            "ksp": {k: dict(v) for k, v in model._summaryKSP.items()},
            "events": list(model._summaryEvents),
        }
        evap = getattr(model, "evapLoss", None)
        if evap is not None:
            rec["evap_loss_total"] = float(evap)
        with open(self.path, "a") as f:
            f.write(json.dumps(rec) + "\n")
        self.step += 1

    def finish(self, model):
        if MPIrank == 0:
            with open(self.path, "a") as f:
                f.write(json.dumps({"end": True, "steps": self.step,
                                    "wall_s": round(MPI.Wtime() - self._run_t0, 3)})
                        + "\n")
