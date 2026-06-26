"""
Command-line entry point for goSPL.

Installed as the ``gospl`` console command (and runnable as ``python -m gospl``)::

    gospl -i input.yml -v
    mpirun -np 8 gospl --input input.yml --verbose

This is exactly the scripted form most users already write::

    from gospl.model import Model
    model = Model(args.input, args.verbose, args.log)   # filename, verbose, showlog
    model.runProcesses()
    model.destroy()

``--fast`` is *not* a CLI flag — fast mode is a property of the run, set with
``domain: fast: true`` in the YAML input.
"""

import argparse


def _parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="gospl",
        description="Run a goSPL landscape-evolution simulation from a YAML "
        "input file. Launch under MPI with e.g. `mpirun -np 8 gospl -i in.yml`.",
    )
    # `--i`/`--v` are explicit aliases (exact match) so they keep working
    # alongside `--version` — an abbreviation `--v` would otherwise be ambiguous
    # between `--verbose` and `--version`. `dest` is pinned so `args.input` /
    # `args.verbose` stay stable regardless of the alias order.
    parser.add_argument("-i", "--i", "--input", dest="input", required=True,
                        metavar="FILE", help="goSPL YAML input file")
    parser.add_argument("-v", "--v", "--verbose", dest="verbose",
                        action="store_true", help="print per-step model progress")
    parser.add_argument("--log", action="store_true",
                        help="write the PETSc solver log summary (showlog)")
    parser.add_argument("--profile", action="store_true",
                        help="record per-phase wall-clock profiling (profile.json)")
    from gospl import __version__
    parser.add_argument("--version", action="version",
                        version="goSPL %s" % __version__)
    return parser.parse_args(argv)


def main(argv=None):
    """Read the input file, run the forward model, then clean up."""
    args = _parse_args(argv)

    # Imported lazily so `gospl --help` / argument errors are fast and any
    # heavy-import failure surfaces with a clear traceback.
    from gospl.model import Model

    model = Model(
        args.input, verbose=args.verbose, showlog=args.log, profile=args.profile
    )
    model.runProcesses()
    model.destroy()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
