"""Enable ``python -m gospl`` — delegates to the :mod:`gospl.cli` entry point."""

from gospl.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
