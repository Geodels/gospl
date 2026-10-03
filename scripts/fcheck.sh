#!/usr/bin/env bash
# Run the goSPL suites against a bounds-checked build of the Fortran extension.
#
# An out-of-bounds write in fortran/functions.F90 is silent on macOS and shows
# up on Linux/glibc as a "wandering" SIGABRT at the next malloc, nowhere near
# the culprit (AGENTS.md > fortran/AGENTS.md: the faceVel and globalngbhs
# bugs). gfortran's -fcheck=all turns it into
#   Fortran runtime error: Index 13 of dimension 2 of array 'fvgnid' above upper bound of 12
# with a backtrace into functions.F90, reproducible on any platform.
#
# Usage:
#   scripts/fcheck.sh                     # tests/ (default)
#   scripts/fcheck.sh tests/ -m flow      # any pytest arguments
#   scripts/fcheck.sh --benchmarks        # tests/ + benchmarks/ -m benchmark
#
# Env: BUILD_DIR (default: build/cp<pyver>, the editable-install build dir).
#
# The build is ALWAYS restored to the normal flags on exit (trap), even on
# failure or Ctrl-C, because a bounds-checked extension is several times
# slower. Requires an editable install (pip install -e . --no-build-isolation).
set -euo pipefail

cd "$(dirname "$0")/.."
PYTAG=$(python -c 'import sys; print(f"cp{sys.version_info.major}{sys.version_info.minor}")')
BUILD_DIR=${BUILD_DIR:-build/$PYTAG}
if [[ ! -f "$BUILD_DIR/build.ninja" ]]; then
    echo "fcheck: no meson build dir at $BUILD_DIR (editable install required;" \
         "set BUILD_DIR=...)" >&2
    exit 2
fi

ORIG_ARGS=$(meson introspect "$BUILD_DIR" --buildoptions |
    python -c 'import json,shlex,sys; v=[o["value"] for o in json.load(sys.stdin) if o["name"]=="fortran_args"][0]; print(shlex.join(v) if isinstance(v, list) else v)')

restore() {
    echo "fcheck: restoring fortran_args='${ORIG_ARGS}' and rebuilding" >&2
    meson configure "$BUILD_DIR" "-Dfortran_args=${ORIG_ARGS}" >/dev/null
    ninja -C "$BUILD_DIR" >/dev/null
}
trap restore EXIT

echo "fcheck: rebuilding $BUILD_DIR with -fcheck=all -fbacktrace" >&2
meson configure "$BUILD_DIR" "-Dfortran_args=-fcheck=all -fbacktrace" >/dev/null
ninja -C "$BUILD_DIR" >/dev/null

if [[ "${1:-}" == "--benchmarks" ]]; then
    shift
    python -m pytest tests/ -q -p no:cacheprovider "$@"
    python -m pytest benchmarks/ -m benchmark -q -p no:cacheprovider "$@"
elif [[ $# -gt 0 ]]; then
    python -m pytest -q -p no:cacheprovider "$@"
else
    python -m pytest tests/ -q -p no:cacheprovider
fi
