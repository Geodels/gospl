#!/usr/bin/env python3
"""
goSPL parallel-scaling driver.

Runs a goSPL model for a bounded number of timesteps with the wall-clock phase
profiler enabled, then writes one machine-readable JSON record per MPI rank
count. Feed several of those records (from runs at different ``mpirun -np``) to
``analyze_scaling.py`` to get strong-/weak-scaling speedup, efficiency and a
per-phase breakdown.

Why a dedicated driver (rather than just running the model)?
  * It **caps the step count** (``--steps``) so a scaling sweep costs a fixed,
    small amount of compute regardless of the input's ``time:end``.
  * It can **isolate compute from I/O** (``--io off``) by pushing the output
    cadence past the (shortened) end time, so the scaling numbers reflect the
    solver/communication cost, not HDF5 write bandwidth. Run once with
    ``--io on`` to measure the I/O contribution too.
  * It records **peak RSS per rank** alongside the timings, so you can see the
    memory-per-rank trend (the global-array anti-patterns show up here).

Invocation (under mpirun/srun)::

    mpirun -np 8 python run_scaling.py -i input.yml --steps 20 --io off \
        --outdir scaling_results

Produces ``scaling_results/scaling_p0008.json``.

All ranks participate (the profiler reduction is collective); only rank 0
writes the JSON.
"""

import os
import sys
import json
import time
import argparse
import platform
import resource

from mpi4py import MPI

# goSPL writes to the cwd-relative paths in the YAML, so run from the input's
# directory (the PBS/Slurm wrappers cd there).
from gospl.model import Model


def _peak_rss_mb():
    """Peak resident set size of this process, in MiB (platform-aware)."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # ru_maxrss is bytes on macOS, kibibytes on Linux.
    if platform.system() == "Darwin":
        return rss / (1024.0 * 1024.0)
    return rss / 1024.0


def main(argv=None):
    parser = argparse.ArgumentParser(description="goSPL scaling driver")
    parser.add_argument("-i", "--input", required=True, help="goSPL input YAML")
    parser.add_argument(
        "--steps",
        type=int,
        default=20,
        help="number of timesteps to run (caps tEnd); <=0 runs to the YAML end",
    )
    parser.add_argument(
        "--io",
        choices=("on", "off"),
        default="off",
        help="'off' suppresses per-step HDF5 output to isolate compute",
    )
    parser.add_argument(
        "--outdir",
        default="scaling_results",
        help="directory for the scaling_p<N>.json record (rank 0)",
    )
    parser.add_argument(
        "--tag",
        default="",
        help="optional label stored in the record (e.g. 'native' / 'container')",
    )
    args = parser.parse_args(argv)

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # ---- model construction (counts as "init", measured separately) ----------
    t0 = MPI.Wtime()
    model = Model(args.input, verbose=False, profile=True)
    init_wall = MPI.Wtime() - t0

    # ---- bound the run -------------------------------------------------------
    if args.steps and args.steps > 0:
        # Cap the run at `steps` timesteps, but never extend past the input's
        # own end time: the forcing interpolants (sea level, tectonics, rain)
        # are only defined over the YAML's time range, so running beyond it
        # raises a scipy "above the interpolation range" ValueError.
        capped = model.tStart + args.steps * model.dt
        model.tEnd = min(model.tEnd, capped)
        if rank == 0 and capped > model.tEnd:
            print(
                "[scaling] --steps %d (%.0f) exceeds the input end time; "
                "capping at tEnd=%.0f" % (args.steps, capped, model.tEnd),
                flush=True,
            )
    if args.io == "off":
        # Push every output trigger past the (shortened) end time so visModel
        # and the stratal writer never fire. The single step-0 write that
        # `saveTime == tStart` would normally force is also suppressed.
        big = model.tEnd + abs(model.tEnd) + 1.0e12
        model.saveTime = big
        model.saveStrat = big

    # ---- timed run -----------------------------------------------------------
    comm.Barrier()
    t0 = MPI.Wtime()
    model.runProcesses()
    run_wall = MPI.Wtime() - t0

    # Phase breakdown (collective: every rank calls reduce()). Tolerate gospl
    # builds that predate the wall-clock profiler (no `model.profiler`): the
    # sweep still reports wall time / speedup / RSS, just without the per-phase
    # split. All ranks run the same build, so this branch is consistent across
    # the communicator and the collective reduce stays balanced.
    profiler = getattr(model, "profiler", None)
    if profiler is not None:
        stats = profiler.reduce(total_wall=run_wall)
    else:
        if rank == 0:
            print(
                "[scaling] this gospl build has no wall-clock profiler; "
                "recording totals only (no per-phase breakdown).",
                flush=True,
            )
        stats = {}

    # Peak memory across ranks.
    rss = _peak_rss_mb()
    rss_arr = comm.gather(rss, root=0)

    # Reduce the two wall times to their max (slowest rank governs).
    run_wall = comm.allreduce(run_wall, op=MPI.MAX)
    init_wall = comm.allreduce(init_wall, op=MPI.MAX)

    model.destroy()

    if rank == 0:
        record = {
            "nranks": size,
            "input": os.path.abspath(args.input),
            "steps": args.steps,
            "io": args.io,
            "tag": args.tag,
            "host": platform.node(),
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "init_wall_s": init_wall,
            "run_wall_s": run_wall,
            "rss_max_mb_per_rank": max(rss_arr),
            "rss_sum_mb": sum(rss_arr),
            "phases": {k: v for k, v in stats.items() if k != "__meta__"},
        }
        os.makedirs(args.outdir, exist_ok=True)
        path = os.path.join(args.outdir, "scaling_p%04d.json" % size)
        with open(path, "w") as fh:
            json.dump(record, fh, indent=2)
        print(
            "[scaling] P=%d  run=%.2fs  init=%.2fs  rss/rank=%.0fMB  -> %s"
            % (size, run_wall, init_wall, max(rss_arr), path),
            flush=True,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
