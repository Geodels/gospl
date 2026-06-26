#!/usr/bin/env python3
"""Thin shim — the tool now lives in the goSPL package and installs as the
``gospl-ela`` command. Kept so existing `python scripts/ela_from_temperature.py`
invocations (and imports of ``derive_ela``) still work."""
from gospl.tools.ela_from_temperature import derive_ela, main  # noqa: F401

if __name__ == "__main__":
    raise SystemExit(main())
