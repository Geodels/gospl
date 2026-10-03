---
name: add-feature
description: Add a new opt-in physical process, model option, output field, forcing or solver to goSPL without breaking the init order, MPI safety, the destroy list or byte-identical behaviour when off. Use when implementing a new feature, a new YAML key/block, a new output or a new PETSc solver.
---

# Add a feature to goSPL

Read the root `AGENTS.md` and the subsystem `AGENTS.md` for every directory you will
touch. Step-specific guides: `docs/HOW_TO_ADD_FORCING.md`, `docs/HOW_TO_ADD_OUTPUT.md`.

## 1. Design
- Opt-in, behind a flag parsed from YAML; **byte-identical when off**. Write down the
  conservation law or analytical limit it should satisfy: that is your test.
- Long features get a `docs/DESIGN_<NAME>.md` with phases (see the existing ones).

## 2. Parse the input (`gospl/tools/inputparser.py`, high-risk)
- New block → a `_readX` / `_extraX` method in the chain; the `_extra*` methods are
  mandatory continuations, so call yours from the right parent.
- Defaults via `self._get_param(...)` or `section.get(key, default)`; never
  `try/except KeyError` for a default. Required keys keep their diagnostic + raise.
- A new forcing DataFrame column is read by NAME (`df.at[nb, col]`), never `iloc`.
- Constants go in `gospl/tools/constants.py` AND the AGENTS.md Magic numbers table.

## 3. State and init order (`gospl/model.py`)
- Allocate in the mixin's `__init__` ONLY when the flag is on; set cached solver /
  operator attributes to `None` unconditionally so `destroy_DMPlex` is safe.
- Check the init-order table (AGENTS.md > The Model god-class): you may only read
  attributes set by earlier mixins. Never switch the chain to `super()`.
- Every persistent Vec/Mat/KSP/SNES/TS → the `destroy_DMPlex` list in
  `mesher/unstructuredmesh.py`.

## 4. Numerics
- Hot-path solver → CACHED pattern (create lazily, store `self._x`, own options
  prefix, never destroy per step). Nested fieldsplit → AD-HOC (create/destroy per
  call). Add the row to the AGENTS.md solver table.
- Assembly: owned rows only, or one of the three safe patterns (MPI contract #2).
- Collectives on every rank; a rank-local decision gating one is reduced first.
  Run `python scripts/lint_mpi_collectives.py`.
- Scratch Vecs (`tmp`, `tmpL`, `tmp1`, `h`, `hl`, `dh`, ...): document which you use in
  the docstring and leave them in a defined state.
- Erosion code reads the `i`-suffix drainage arrays (`rcvIDi`, `wghtVali`, `fMati`).
- Keep `Eb`/`EbLocal` in the thickness-rate convention (positive = deposition).

## 5. Outputs and restart
Follow `docs/HOW_TO_ADD_OUTPUT.md`: HDF5 dataset AND XDMF entry, gated on the flag;
restore model-memory state in `outmesh.readData` (robust to older outputs that lack
it). Add a `gospl/tools/runsummary.py` event or field if the feature has a solver
outcome worth surfacing.

## 6. Tests (in the module for the subsystem, see tests/README.md)
- Opt-in parse + allocation test; off-path byte-identity test.
- Conservation / invariant test that is red without the feature's key line.
- `mpi`-marked np=1 vs np=2 test if the feature communicates
  (`tests/test_parallel.py`), and `python scripts/ab_partition.py` on its fixture.
- An analytical benchmark in `benchmarks/` if a closed form exists
  (`benchmarks/AGENTS.md`).
- `pytest tests/ -n 4` and `pytest benchmarks/ -m benchmark` green.

## 7. Docs
User keys in `docs/user_guide/surfproc.rst` (or the relevant page), theory in
`docs/tech_guide/`, every new method in BOTH the `autosummary` and `automethod` lists
of its `docs/api_ref/*_ref.rst`. Update the subsystem `AGENTS.md`, add a
`docs/dev/CHANGELOG_DEV.md` row.
