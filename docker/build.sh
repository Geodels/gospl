#!/usr/bin/env bash
###############################################################################
# build.sh — build the goSPL HPC container and convert it to a .sif for
# NCI Gadi / Pawsey Setonix. Run on a LOCAL Linux machine with Docker; you
# cannot build on Gadi or Setonix (root is required). See AGENTS.md >
# "Docker / HPC Container".
#
#   ./build.sh v2026.06.11           # build + smoke-test + convert to .sif
#   ./build.sh v2026.06.11 --push    # also push to docker.io/geodels/gospl-hpc
###############################################################################
set -euo pipefail

usage() {
    echo "Usage: $0 <tag> [--push]    e.g. $0 v2026.06.11 --push" >&2
    exit 1
}
[ "$#" -ge 1 ] || usage

TAG="$1"; shift || true
PUSH=0
if [ "$#" -ge 1 ]; then
    [ "$1" = "--push" ] || usage
    PUSH=1
fi

# Strip the leading 'v' for the package version — MUST match meson.build line 4.
GOSPL_VERSION="${TAG#v}"
# Docker Hub repo (override with REGISTRY=... ./build.sh ...). Matches the
# docker-build.yml CI workflow that publishes on tag push.
REGISTRY="${REGISTRY:-docker.io/geodels/gospl-hpc}"
IMAGE="${REGISTRY}:${TAG}"
SIF="gospl-hpc-${TAG}.sif"
# Build context is the repo root (one level up from this script).
CONTEXT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo ">>> Building ${IMAGE}  (GOSPL_VERSION=${GOSPL_VERSION})"
docker build \
    --platform linux/amd64 \
    --build-arg "GOSPL_VERSION=${GOSPL_VERSION}" \
    -f "${CONTEXT_DIR}/docker/Dockerfile" \
    -t "${IMAGE}" \
    "${CONTEXT_DIR}"

echo ">>> Smoke-testing image"
docker run --rm --platform linux/amd64 "${IMAGE}" \
    python -c "import gospl, petsc4py, mpi4py; print('goSPL container build OK')"

echo ">>> Verifying built version matches ${GOSPL_VERSION}"
BUILT_VER="$(docker run --rm --platform linux/amd64 "${IMAGE}" \
    python -c 'import gospl; print(gospl.__version__)')"
# gospl.__version__ is PEP 440-normalized (2026.06.11 -> 2026.6.11); normalize
# the requested version the same way before comparing.
WANT_VER="$(python3 -c "print('.'.join(str(int(p)) if p.isdigit() else p for p in '${GOSPL_VERSION}'.split('.')))")"
if [ "${BUILT_VER}" != "${WANT_VER}" ]; then
    echo "ERROR: container goSPL ${BUILT_VER} != ${WANT_VER} (requested ${GOSPL_VERSION})" >&2
    exit 1
fi

if [ "${PUSH}" -eq 1 ]; then
    echo ">>> Pushing ${IMAGE}"
    docker push "${IMAGE}"
fi

echo ">>> Converting to ${SIF}"
if command -v apptainer >/dev/null 2>&1; then
    BUILDER=apptainer
elif command -v singularity >/dev/null 2>&1; then
    BUILDER=singularity
else
    echo "Neither apptainer nor singularity found; skipping .sif conversion." >&2
    echo "Image is available as ${IMAGE}." >&2
    exit 0
fi

if [ "${PUSH}" -eq 1 ]; then
    # Pull from the registry we just pushed to.
    "${BUILDER}" build "${SIF}" "docker://${IMAGE}"
else
    # Convert straight from the local Docker daemon (no registry round-trip).
    "${BUILDER}" build "${SIF}" "docker-daemon://${IMAGE}"
fi

echo ">>> Done: ${SIF}"
echo ">>> scp this .sif to Gadi/Setonix and point CONTAINER= in the job scripts at it."
