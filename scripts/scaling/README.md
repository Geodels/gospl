# goSPL parallel-scaling harness

Tools to measure how goSPL scales across CPU counts and where the wall-clock
time goes. Built around the wall-clock phase profiler (`gospl/tools/profiler.py`,
enabled with `Model(..., profile=True)` or YAML `output: profile: true`).

| File | Role |
|---|---|
| `run_scaling.py` | mpirun/srun entry point: runs a capped number of steps with profiling, writes `scaling_p<N>.json` (timings + speedup inputs + peak RSS). `--steps` caps downward only (never past the input's forcing range); `--verbose` enables per-phase + solver monitors; `--set ATTR=VALUE` overrides model attributes (e.g. `--set soil_solver=qn`) without editing the YAML. |
| `local_sweep.sh` | **Workstation** driver: runs ranks 1–6 (default) sequentially under `mpirun`, writes `scaling_p<N>.json`, then auto-runs `analyze_scaling.py`. The local counterpart of `submit_sweep.sh`. |
| `gadi.pbs` | PBS job for **one** rank count on NCI Gadi (native intel-mpi module build *or* Singularity container). Launches `RANKS` MPI ranks — deliberately **not** `NCPUS` (a reserved Gadi variable auto-set to the node core count, which silently shadowed any override above 48) — and sets `I_MPI_HYDRA_BRANCH_COUNT` so Intel-MPI spans multiple nodes (without it, >48-rank jobs all ran on one node and overwrote `scaling_p0048.json`). |
| `submit_sweep.sh` | qsub's `gadi.pbs` once per rank count, sizing each job's `ncpus`/`mem`, optionally chaining them (`DEPEND=1`) for clean timing. |
| `analyze_scaling.py` | Reads all `scaling_p*.json`, computes speedup/efficiency, writes CSV + Markdown + plots to `results/scaling/`. |
| `plot_scaling.py` | Plots speedup + efficiency from a `scaling_summary.csv` **alone** (when the per-run `scaling_p*.json` are no longer around). Baselines to the smallest rank count; overlays the Amdahl-absolute curve if present. |

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

Easiest — the 1→6 rank sweep + analysis in one go (run from the input's dir):

```bash
STEPS=5 IO=off bash scripts/scaling/local_sweep.sh input.yml
# -> scaling_results/scaling_p000N.json + scaling_results/summary/{csv,md,png}
```

Or drive `run_scaling.py` directly (e.g. to profile / tune a single config):

```bash
# one step, verbose, force a solver choice without editing the YAML:
mpirun -n 1 python scripts/scaling/run_scaling.py \
    -i input.yml --steps 1 --io off --outdir /tmp/sc --verbose --set soil_solver=qn
```

> The bundled fixtures (~3.8k nodes) are **far too small** to scale past a few
> ranks — partitions become dominated by halo/communication. Use them only to
> check the harness runs. For real numbers, point `-i` at a
> production-resolution global-sphere input. On a single workstation, strong
> scaling also plateaus at a few cores (shared memory bandwidth); the true
> verdict needs a multi-node HPC run.

## On NCI Gadi

1. Edit `gadi.pbs`: set `#PBS -P <project>`, `-l storage`, and (native mode) the
   `GOSPL_VENV` path or (container mode) the `CONTAINER` `.sif` path. The native
   path loads the intel-mpi/petsc/hdf5/netcdf module stack (mirrors
   `hpc_setup/nciRun.pbs`).
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
   (default 4 GB/rank, raise it for large meshes). Counts >48 span multiple
   nodes — Gadi's `normal` queue allocates in whole 48-core nodes, so use 48,
   96, 144, 192, … there.

   Add `DEPEND=1` to chain the jobs (one at a time) for the cleanest timing on a
   shared system.

3. After they finish:

   ```bash
   python analyze_scaling.py /scratch/q97/$USER/run/scaling_results -o results/scaling
   ```

   → console table, `results/scaling/scaling_summary.{csv,md}`, and
   `scaling_speedup.png` / `scaling_phases.png`.

## Plotting the results

`analyze_scaling.py` writes the plots directly. If you only kept the summary CSV
(the per-run `scaling_p*.json` are gone), regenerate the speedup/efficiency
figure from the CSV alone:

```bash
python plot_scaling.py scaling_summary.csv scaling.png
```

> **No 1-CPU baseline needed.** A production global model is memory-infeasible at
> low rank counts (it OOMs well before 1 rank — e.g. the 5.9M-node mesh needs
> ~6.5 GB/rank and cannot start under ~24 ranks). Both `analyze_scaling.py` and
> `plot_scaling.py` therefore **baseline to the smallest rank count in the
> sweep** (`p0`): speedup `S(p)=T(p0)/T(p)`, efficiency `E(p)=S·p0/p`, ideal line
> `p/p0`. State "baseline = `p0` ranks" in any caption. Efficiency can read
> **>100 %** at the first few steps up from `p0` when the `p0` run is itself
> memory-bandwidth bound — that is real (superlinear) and expected; the
> Amdahl-absolute columns (`*_abs_amdahl`) give the extrapolated vs-1-rank view
> that avoids the artifact.

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
