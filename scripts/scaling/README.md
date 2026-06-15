# goSPL parallel-scaling harness

Tools to measure how goSPL scales across CPU counts and where the wall-clock
time goes. Built around the wall-clock phase profiler (`gospl/tools/profiler.py`,
enabled with `Model(..., profile=True)` or YAML `output: profile: true`).

| File | Role |
|---|---|
| `run_scaling.py` | mpirun/srun entry point: runs a capped number of steps with profiling, writes `scaling_p<N>.json` (timings + speedup inputs + peak RSS). |
| `gadi_scaling.pbs` | PBS job for **one** rank count on NCI Gadi (native module build *or* Singularity container). |
| `submit_sweep.sh` | qsub's `gadi_scaling.pbs` once per rank count (1,2,4,…,48). |
| `analyze_scaling.py` | Reads all `scaling_p*.json`, computes speedup/efficiency, writes CSV + Markdown + plots to `results/scaling/`. |

## What it measures

- **Strong scaling**: same problem, more ranks → speedup `S(p)=T(base)/T(p)` and
  efficiency `E(p)=S·base/p`.
- **Per-phase wall time** (`flow`, `erosion`, `sed`, `sea`, `hillslope`,
  `flexure`, `ice`, `till`, `strat`, `forcing`, `tectonics`, `output`) with a
  cross-rank **load-imbalance** ratio (`max/mean`; 1.0 = balanced).
- **Peak RSS per rank** — surfaces the global-array (`mpoints`) memory patterns
  that don't shrink as you add ranks.

`--io off` (the default) pushes the output cadence past the shortened end time
so the numbers reflect compute + communication, not HDF5 bandwidth. Do one
sweep with `--io on` to quantify the I/O share.

## Quick local check (workstation)

```bash
cd tests/fixtures           # any dir whose YAML's relative paths resolve
mpirun -n 2 python ../../scripts/scaling/run_scaling.py \
    -i minimal.yml --steps 10 --io off --outdir /tmp/sc
mpirun -n 1 python ../../scripts/scaling/run_scaling.py \
    -i minimal.yml --steps 10 --io off --outdir /tmp/sc
python ../../scripts/scaling/analyze_scaling.py /tmp/sc -o /tmp/sc/out
```

> The bundled fixtures (~3.8k nodes) are **far too small** to scale past a few
> ranks — partitions become dominated by halo/communication. Use them only to
> check the harness runs. For real numbers, point `-i` at a
> production-resolution global-sphere input.

## On NCI Gadi

1. Edit `gadi_scaling.pbs`: set `#PBS -P <project>`, `-l storage`, and (native
   mode) the `GOSPL_VENV` path or (container mode) `CONTAINER` `.sif` path.
2. Submit a sweep (absolute input path; pick a representative production input):

   ```bash
   # native build (mirrors hpc_setup/nciRun.pbs modules)
   PROJECT=q97 GOSPL_VENV=$HOME/envi_gospl/bin/activate \
     ./submit_sweep.sh /scratch/q97/$USER/run/input.yml 1 2 4 8 16 24 48

   # or the Singularity container
   PROJECT=q97 MODE=container CONTAINER=/scratch/q97/$USER/gospl-hpc-v2026.6.13.sif \
     ./submit_sweep.sh /scratch/q97/$USER/run/input.yml 1 2 4 8 16 24 48
   ```

   Each rank count is a separate job sized `mem = ranks × MEM_PER_CPU_GB`
   (default 4 GB/rank). Counts >48 span multiple nodes — Gadi's `normal` queue
   allocates in whole 48-core nodes, so use 48, 96, 192, … there.

   Add `DEPEND=1` to chain the jobs (one at a time) for the cleanest timing on a
   shared system.

3. After they finish:

   ```bash
   python analyze_scaling.py /scratch/q97/$USER/run/scaling_results -o results/scaling
   ```

   → console table, `results/scaling/scaling_summary.{csv,md}`, and
   `scaling_speedup.png` / `scaling_phases.png`.

## Reading the results

- **Efficiency falling off a cliff at some p** → that's where communication /
  serial sections dominate. Cross-reference the per-phase plot: a phase whose
  *mean* time stops dropping (or grows) with more ranks is the culprit.
- **High imbalance (`max/mean` ≫ 1)** on a phase → load imbalance (uneven
  partition work, e.g. all the ocean/ice on a few ranks).
- **`rss_max_mb_per_rank` flat or rising with p** → a global-sized
  (`mpoints`) allocation that doesn't decompose (the known flexure / advection /
  marine `_matOcean` patterns — see the parallel-performance plan).

Use `--steps` large enough that per-step cost dwarfs init, but small enough to
keep the sweep cheap; 20–50 steps is usually plenty.
