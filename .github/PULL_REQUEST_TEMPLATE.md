## Summary

<!-- What does this change and why? Link related issues: Fixes #123 -->

## Type of change

- [ ] Bug fix
- [ ] New feature (opt-in, off by default)
- [ ] Refactor / performance
- [ ] Documentation / examples
- [ ] Build, CI or release

## Checklist

- [ ] `pytest tests/ -n 4` passes
- [ ] Tests added or updated, ideally against a conservation law or analytical benchmark
- [ ] Docs updated (user guide, API, `docs/dev/CHANGELOG_DEV.md`)
- [ ] New features are opt-in and **byte-identical when off**

<details>
<summary>If your change touches MPI, solvers or Fortran</summary>

- [ ] `python scripts/lint_mpi_collectives.py` is clean (any `# mpi-lint: ok` has a reason)
- [ ] Checked with `pytest tests/ -m mpi` and `scripts/ab_partition.py`
- [ ] New Vec/Mat/KSP objects added to the `destroy_DMPlex` list
- [ ] `functions.F90` and `functions.pyf` updated together, `scripts/fcheck.sh` run

</details>

## Testing

<!-- What you ran, on how many ranks, and on which input. -->
