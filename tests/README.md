# goSPL test suite

Regression tests that guard the invariants in `AGENTS.md`. The analytical
benchmarks live separately in `benchmarks/` (`pytest benchmarks/ -m benchmark`).

## Which tests to run

| You changed… | Run | Time (M-series laptop) |
|---|---|---|
| anything, quick check | `pytest tests/ -n 4 -m "not mpi"` | ~25 s |
| one subsystem | `pytest tests/ -m <marker>` (table below) | 3–15 s |
| a collective, a solver, assembly, halo sync | `pytest tests/ -m mpi` | ~70 s |
| before a commit | `pytest tests/ -n 4` | ~80 s |
| serial reference (what CI runs) | `pytest tests/` | ~2.5 min |

`-n 4` needs `pytest-xdist` (in `environment.yml`). The suite is xdist-safe:
every process gets a private copy of `tests/fixtures/` (see
`_helpers._private_fixtures_dir`), so parallel workers cannot delete each
other's model output.

## Markers

Each `test_<group>.py` module sets `pytestmark`, so every test carries exactly
one subsystem marker. They are registered in `pyproject.toml` (authoritative).

| Marker | Module(s) | Covers |
|---|---|---|
| `parser` | `test_parser.py` | forcing DataFrames, YAML parsing (no mesh) |
| `contracts` | `test_contracts.py` | `rcvIDi` snapshot, scratch Vecs, `Eb` sign |
| `flow` | `test_flow.py`, `test_fill_mesh_pits.py` | flow KSP, pits, cascade, evaporation, `globalngbhs` |
| `sediment` | `test_sediment.py` | mass conservation, marine routing/diffusion, `oFill` |
| `boundaries` | `test_boundaries.py` | `o/f/w/c` edges, closed-sink deposition, cyclic |
| `advection` | `test_advection.py` | advection schemes, tectonic events, stratigraphy moved by horizontal motion (flat, cyclic, sphere; serial) |
| `flexure` | `test_flexure.py` | flat FEM, global DH grid |
| `orography` | `test_orography.py` | orographic rain |
| `ice` | `test_ice.py` | diagnostic glacial model, till, meltwater |
| `soil` | `test_soil.py` | soil/regolith |
| `groundwater` | `test_groundwater.py` | water table, recharge, baseflow, duricrust |
| `geochem` | `test_geochem.py` | Level-B solute transport |
| `strata` | `test_strata.py` | stratigraphy, dual lithology |
| `provenance` | `test_provenance_tracers.py` | in-model provenance tracers |
| `mpi` (+`parallel`) | `test_parallel.py` | `mpirun` np=1 vs np=2: partition invariance, deadlock guards, no whole-mesh copies per rank |
| `analyse` | `test_provenance.py`, `test_gridexport.py`, `test_stratamesh.py`, `test_stratasection.py`, `test_catchment.py` | post-processing tools |
| `tools` | `test_profiler.py`, `test_ela_from_temperature.py` | `gospl.tools` utilities |
| `cli` | `test_cli.py` | the `gospl` command |

The orthogonal `slow` marker (a full `Model` instantiation) predates the
subsystem split and is not applied consistently; prefer the subsystem markers.

## Mesh distribution and horizontal advection guards (2026-11)

| Test | Module | Protects |
|---|---|---|
| `test_mesh_not_replicated_per_rank` | `test_parallel.py` | only rank 0 reads the input mesh; `locIDs`/`glbIDs`/initial elevation delivered by the distribution; no `mCoords` copy and local-sized forcing maps on every rank |
| `test_plate_interp_tables_local` | `test_parallel.py` | plate-interpolation tables are local; np=1 vs np=2 elevation |
| `test_plate_strata_partition_invariant[interp, iioe2]` | `test_parallel.py` | the stratigraphy move fetches source rows from their owners: bitwise partition-invariant |
| `test_plate_advection_moves_stratigraphy` | `test_advection.py` | plate interpolation moves every layer field (independent numpy reference) |
| `test_scheme_advection_moves_stratigraphy[upwind, iioe1, iioe2]` | `test_advection.py` | the schemes move the layers: exact one-cell shift on a regular grid |
| `test_cyclic_scheme_advection_moves_stratigraphy_across_seam[upwind, iioe2]` | `test_advection.py` | layers wrap across the periodic seam |
| `test_cyclic_plate_advection_wraps_across_seam` | `test_advection.py` | plate interpolation moves around the cylinder (surface and layers) |
| `test_sphere_scheme_advection_moves_stratigraphy` | `test_advection.py` | the layer move on a global mesh (rigid rotation) |
| `test_tectonic_event_without_maps_runs` | `test_advection.py` | a tectonic event with only `start`/`end` runs |
| `test_rain_events_prepend_without_futurewarning` | `test_parser.py` | forcing tables built without `pd.concat` (no pandas 2 warning per rank) |

## Writing a test

* Put it in the module for its subsystem. If none fits, add a module and a
  marker (register it in `pyproject.toml` and the table above).
* Use the model fixtures in `conftest.py` (`minimal_model`, `incising_model`,
  …) or `_helpers._gw_model`. Always tear a hand-built `Model` down with
  `model.destroy()` in a `try/finally`.
* Fixture inputs: use `FIXTURES_DIR` from `_helpers`, never the literal
  `tests/fixtures` path (that would bypass the per-process copy). Set
  `GOSPL_TEST_FIXTURES_INPLACE=1` to run in `tests/fixtures/` when you want to
  inspect a test's model output afterwards.
* Multi-rank tests: build the child environment with `_helpers.mpi_child_env()`
  (scrubs the OpenMPI variables that make a nested `mpirun` refuse to launch,
  and sets `OMPI_MCA_rmaps_base_oversubscribe=1`: GitHub's ubuntu runners have
  4 vCPUs but 2 physical cores, and OpenMPI refuses np=3 without it)
  and pass `timeout=MPI_TIMEOUT`. A run past the timeout is almost always a
  deadlock (AGENTS.md > MPI contract). Raise it on a slow machine with
  `GOSPL_TEST_MPI_TIMEOUT=<seconds>`.
* A regression guard should fail without the fix. Check that before you commit
  (revert the fix locally, confirm red, restore).
