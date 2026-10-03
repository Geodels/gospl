---
name: parallel-check
description: Diagnose or rule out a goSPL MPI problem - a run that hangs at np>1 (deadlock), results that change with the number of ranks, a spike along partition boundaries, or a solver that is slow only at some rank counts. Also run it after adding or moving any collective (localToGlobal, Allreduce, ksp.solve, vec.norm/sum, garbage_cleanup) or a matrix assembly.
---

# Parallel check

The two bug classes and their rules: root `AGENTS.md` > MPI contract. Flow-specific
partition invariance (drainage graph, near-singular cells): `gospl/flow/AGENTS.md`.

## A. Static check first (seconds)
```bash
python scripts/lint_mpi_collectives.py          # collectives under rank-local guards
```
Every finding is a question: either reduce the rank-local flag before it gates the
collective (`flag = MPI.COMM_WORLD.allreduce(local, op=MPI.LOR)`), compute the
collective on every rank and gate only the print/IO, or, if the condition is provably
global (config scalar, per-pit array from `pitParams`, an already-reduced value),
annotate the guard line `# mpi-lint: ok <why it is global>`.
The lint cannot see through helper calls: also read the callers of anything you
changed.

## B. Multi-rank tests (~70 s)
```bash
pytest tests/ -m mpi
```
A hang shows up as a failure after `MPI_TIMEOUT` (240 s; `GOSPL_TEST_MPI_TIMEOUT`).
Serial always passes for a deadlock; that asymmetry is the signature.

## C. Compare decompositions node by node
```bash
python scripts/ab_partition.py <input.yml> -n 4 --steps 2        # np=1 vs np=4
python scripts/ab_partition.py <input.yml> -n 4 --steps 6        # does it GROW?
python scripts/ab_partition.py <input.yml> -n 4 --field iceHL --json
```
A healthy model does not agree node-by-node: KSP round-off flips near-tie routing.
Baseline on `tests/fixtures/minimal.yml`, np=1 vs np=2, 2 steps: elevation max rel
~4e-6, discharge max rel ~0.3 at ~175 nodes. Suspect a real bug when the difference
is orders of magnitude above the baseline for your input, is localised along the
partition boundary, appears as non-finite values on one side only, or grows with
`--steps`.

## D. Matrix or solver?
If a field is partition-dependent, swap the suspect solve's PC for an exact one via
its options prefix (e.g. `PETSC_OPTIONS="-<prefix>ksp_type preonly -<prefix>pc_type lu"`,
MUMPS at np>1) and rerun C. If the difference survives an exact solve, the
**matrix** is partition-dependent: an assembly writing GHOST rows with
`INSERT_VALUES` (AGENTS.md > the #2 class). Restrict the loop to owned rows
(`self.glIDs`) or use one of the three safe patterns. Only chase KSP/PC/tolerances
once assembly is confirmed owned-rows-only.

## E. A hang in a real run
- With np>1 an uncaught exception now aborts every rank (`_install_mpi_abort_excepthook`);
  a true hang is a collective mismatch. Run the window with `-v` and compare the last
  line printed per rank (`mpirun --output-filename logs`), then look for the first
  collective one rank reaches and another does not.
- A serial rank-0 step inside `Model.__init__` (e.g. the global DH grid) also
  presents as an init-time hang at np>1 when it is slow.

## F. Guard it
Add an `mpi`-marked test in `tests/test_parallel.py` (use `_helpers.mpi_child_env()`
and `timeout=MPI_TIMEOUT`) that hangs or fails without the fix.
