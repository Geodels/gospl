"""
Lightweight wall-clock phase profiler for goSPL.

Motivation
----------
The pre-existing per-phase ``process_time()`` prints in goSPL measure
per-process *CPU* time, gated behind ``verbose``. CPU time is the wrong clock
for a parallel run: under MPI it either hides time a rank spends blocked in a
collective (sleep-wait) or counts busy-wait spinning as if it were useful
compute. Either way it makes load imbalance invisible, which is exactly the
thing you need to see when scaling across CPUs.

This module provides :class:`Profiler`, a small wall-clock timer that:

* uses ``MPI.Wtime()`` (the MPI standard wall clock),
* accumulates elapsed seconds and call counts per named phase,
* does **no** collective communication inside the time loop — the only
  reduction happens once, in :meth:`report` — so per-step overhead is two
  ``Wtime`` calls per phase,
* is a true no-op when disabled (``phase`` returns a shared ``nullcontext``),
* at the end reduces min/mean/max across ranks and prints a load-imbalance
  table on rank 0, also writing a machine-readable ``profile.json`` consumed by
  the scaling harness (``scripts/scaling/``).

Usage
-----
::

    self.profiler = Profiler(enabled=True)
    with self.profiler.phase("flow"):
        self.flowAccumulation()
    ...
    self.profiler.report(outputDir, total_wall=elapsed)

When disabled the ``with`` block costs nothing measurable.
"""

import os
import json
from contextlib import nullcontext

import numpy as np
from mpi4py import MPI


class Profiler(object):
    """Accumulate per-phase wall-clock time and report cross-rank statistics.

    :arg enabled: when ``False`` every method is a no-op (zero overhead).
    :arg comm: MPI communicator (defaults to ``MPI.COMM_WORLD``).
    """

    class _Timer(object):
        """Context manager that charges its elapsed wall time to one phase."""

        __slots__ = ("_prof", "_name", "_t0")

        def __init__(self, prof, name):
            self._prof = prof
            self._name = name

        def __enter__(self):
            self._t0 = MPI.Wtime()
            return self

        def __exit__(self, exc_type, exc_val, exc_tb):
            dt = MPI.Wtime() - self._t0
            self._prof._accumulate(self._name, dt)
            return False

    def __init__(self, enabled=False, comm=None):
        self.enabled = bool(enabled)
        self.comm = comm if comm is not None else MPI.COMM_WORLD
        self.times = {}
        self.counts = {}
        self._order = []
        # Shared zero-cost context returned when profiling is off.
        self._null = nullcontext()

    def _accumulate(self, name, dt):
        if name not in self.times:
            self.times[name] = 0.0
            self.counts[name] = 0
            self._order.append(name)
        self.times[name] += dt
        self.counts[name] += 1

    def phase(self, name):
        """Return a context manager timing ``name`` (no-op when disabled)."""
        if not self.enabled:
            return self._null
        return self._Timer(self, name)

    def start(self, name):
        """Imperative timer start; pair with :meth:`stop`."""
        if not self.enabled:
            return
        self._starts = getattr(self, "_starts", {})
        self._starts[name] = MPI.Wtime()

    def stop(self, name):
        """Imperative timer stop; pair with :meth:`start`."""
        if not self.enabled:
            return
        t0 = getattr(self, "_starts", {}).get(name)
        if t0 is not None:
            self._accumulate(name, MPI.Wtime() - t0)

    def reduce(self, total_wall=None):
        """Reduce the per-phase timings across all ranks.

        Collective: every rank MUST call this. Returns a dict (identical on
        every rank) keyed by phase name with min/mean/max/imbalance/calls, plus
        a ``__meta__`` entry. Returns ``None`` when disabled.
        """
        if not self.enabled:
            return None

        comm = self.comm
        size = comm.Get_size()

        # Union of phase names across ranks (a rank may skip a conditional
        # phase, e.g. flexure-off), preserving first-seen order on rank 0.
        gathered = comm.allgather(self._order)
        names = []
        seen = set()
        for lst in gathered:
            for n in lst:
                if n not in seen:
                    seen.add(n)
                    names.append(n)

        local = np.array([self.times.get(n, 0.0) for n in names], dtype=np.float64)
        tmin = local.copy()
        tmax = local.copy()
        tsum = local.copy()
        comm.Allreduce(MPI.IN_PLACE, tmin, op=MPI.MIN)
        comm.Allreduce(MPI.IN_PLACE, tmax, op=MPI.MAX)
        comm.Allreduce(MPI.IN_PLACE, tsum, op=MPI.SUM)
        tmean = tsum / size

        lcnt = np.array([self.counts.get(n, 0) for n in names], dtype=np.int64)
        cmax = lcnt.copy()
        comm.Allreduce(MPI.IN_PLACE, cmax, op=MPI.MAX)

        if total_wall is not None:
            twall = np.array([total_wall], dtype=np.float64)
            comm.Allreduce(MPI.IN_PLACE, twall, op=MPI.MAX)
            total_wall = float(twall[0])

        out = {"__meta__": {"nranks": size, "total_wall": total_wall}}
        for i, n in enumerate(names):
            mean = float(tmean[i])
            mx = float(tmax[i])
            out[n] = {
                "min": float(tmin[i]),
                "mean": mean,
                "max": mx,
                "sum": float(tsum[i]),
                # max/mean: 1.0 = perfectly balanced, >1 = imbalance.
                "imbalance": (mx / mean) if mean > 0.0 else 0.0,
                "calls": int(cmax[i]),
            }
        return out

    def report(self, outputDir=None, filename="profile.json",
               total_wall=None, to_stdout=True):
        """Reduce, print a table (rank 0), and write ``profile.json``.

        Collective: every rank MUST call this (it calls :meth:`reduce`).
        """
        stats = self.reduce(total_wall=total_wall)
        if stats is None:
            return None

        if self.comm.Get_rank() != 0:
            return stats

        meta = stats["__meta__"]
        phases = [(k, v) for k, v in stats.items() if k != "__meta__"]
        # Sort by mean wall time, descending — biggest cost first.
        phases.sort(key=lambda kv: kv[1]["mean"], reverse=True)

        sum_mean = sum(v["mean"] for _, v in phases)
        denom = meta["total_wall"] if meta["total_wall"] else sum_mean

        if to_stdout:
            nproc = meta["nranks"]
            print("", flush=True)
            print(
                "================ goSPL wall-clock profile (P=%d ranks) ================"
                % nproc,
                flush=True,
            )
            print(
                "%-14s %10s %7s %10s %10s %8s"
                % ("phase", "mean(s)", "%", "max(s)", "imbal", "calls"),
                flush=True,
            )
            print("-" * 64, flush=True)
            for name, v in phases:
                pct = (100.0 * v["mean"] / denom) if denom > 0 else 0.0
                print(
                    "%-14s %10.3f %6.1f%% %10.3f %10.2f %8d"
                    % (name, v["mean"], pct, v["max"], v["imbalance"], v["calls"]),
                    flush=True,
                )
            print("-" * 64, flush=True)
            if meta["total_wall"]:
                unacc = meta["total_wall"] - sum_mean
                print("%-14s %10.3f" % ("TOTAL(wall)", meta["total_wall"]), flush=True)
                print(
                    "%-14s %10.3f %6.1f%%"
                    % ("unaccounted", unacc, 100.0 * unacc / meta["total_wall"]),
                    flush=True,
                )
            print("=" * 64, flush=True)

        if outputDir is not None:
            try:
                path = os.path.join(outputDir, filename)
                with open(path, "w") as fh:
                    json.dump(stats, fh, indent=2)
                if to_stdout:
                    print("Wrote profile to %s" % path, flush=True)
            except (IOError, OSError) as err:
                if to_stdout:
                    print("Could not write %s: %s" % (filename, err), flush=True)

        return stats
