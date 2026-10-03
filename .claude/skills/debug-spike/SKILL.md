---
name: debug-spike
description: Find which goSPL process produced an elevation spike, a needle, a km-scale deposit or erosion pit, NaNs, or a mass imbalance. Use when a run shows an unphysical elevation extreme, a sudden jump in zmax/zmin, NaN/inf in outputs, or a solver that suddenly slows down.
---

# Debug an elevation spike / mass imbalance

The failure class: a km-scale needle at one cell, noticed several outputs later,
when the output alone cannot say which stage produced it. Background:
`gospl/tools/AGENTS.md` > Debugging an elevation spike, `docs/dev/FIXED_BUGS.md`.

## 1. Locate it in time and space (no rerun needed)
```bash
gospl-inspect <outdir>              # largest |dz| between outputs + node id + xyz
gospl-inspect <outdir> --json       # same, machine-readable
```
Note the first output step where the jump appears and the input-mesh node id.

## 2. Rerun the window with the instruments on
Shorten the input to a few steps before and after the jump (copy the YAML; set
`time: start/end` around it, `tout = dt` so every step is written, restart from the
nearest output with `time: rstep: <output index>` if the run is long), then:
```bash
gospl -i window.yml -v --summary run.jsonl > run.log 2>&1
gospl-inspect --summary run.jsonl   # nonfinite, solver failures, elev jumps, budget drift
grep -n "zprobe\|mass rescale\|fDep at the\|RAISED THE GLOBAL MAXIMUM" run.log
```
`zprobe` (verbose only) prints the global elevation extremes after each stage and
tags the stage that raised the maximum; `reportVolume` gives each unbounded deposit's
volume AND peak thickness.

## 3. Check the usual causes, in this order
1. **`fDep` cap** (`spl: G` > 0): if `zprobe` reports many cells pinned at the 0.99
   cap, the coupled SPL block is near-singular. Rule: `G ≲ rain/2` (rain in m/yr);
   the criterion is mesh-resolution independent. Also look for a step-by-step rise of
   "Solve SPL accounting for sediment deposition" time.
2. **Unbounded deposition sites** (each conserves mass when routing has nowhere else
   to go): closed-sink deposit `_closedDepo` (always active on a sphere), marine
   force-deposit at terminal sinks and the cascade-exit residual, the `_diffuseOcean`
   global mass rescale. `zprobe` instruments these.
3. **Flow solver**: `run.jsonl` events `flow_ksp_fallback_failed` (ponded / zeroed /
   abort, with `nbad` and `worst_mesh_id`) and `flow_cascade` with outcome `stall`.
   A large un-drained region is usually a topography problem: run
   `python scripts/fill_mesh_pits.py --check <mesh.npz>` before touching solvers.
4. **Soil / eroder SNES**: `soil_substep` (converged via Δt/N, or discarded),
   `*_snes_fallback`, `nlspl_snes_failed`.
5. **Not instrumented by zprobe**: ice/till, soil production and creep,
   groundwater/duricrust, tectonics/advection, nlSPL/soilSPL erosion. Add a
   `probeZ(model, "tag")` call around the suspect stage (collective; gate it ONLY on
   `self.verbose`) and extend the zprobe docstring's scope list.

## 4. Parallel-only?
If np=1 is clean and np>1 spikes, switch to the `parallel-check` skill: the cause is
usually a partition-dependent operator (ghost rows), not the physics.

## 5. Fix and guard
Write a regression test that is red without the fix (tests/README.md), add a
`docs/dev/FIXED_BUGS.md` entry (symptom, cause, fix, guard), and a
`docs/dev/CHANGELOG_DEV.md` row.
