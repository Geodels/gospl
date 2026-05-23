#!/usr/bin/env bash
# Build script invoked by conda-build for linux-64 and osx-arm64.
# The actual install is driven by meta.yaml's `script:` entry point; this
# file exists so conda-build does not fall back to a default behaviour
# and is a convenient place to hook platform-specific tweaks.

set -euo pipefail

# Make sure the conda-provided gfortran is the one meson detects. On
# osx-arm64 the system /usr/bin/clang shadows the conda compiler shim
# unless we keep the prefix bin on PATH first (conda-build already does
# this, but we re-assert it here to fail loudly if the assumption breaks).
export PATH="${PREFIX}/bin:${PATH}"

# meson-python build via pip (matches the script: in meta.yaml).
"${PYTHON}" -m pip install . \
    --no-build-isolation \
    --no-deps \
    -vv
