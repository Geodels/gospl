#!/bin/bash
###############################################################################
# Run a goSPL strong-scaling sweep on a local workstation (no PBS/Slurm).
#
# Runs run_scaling.py once per MPI rank count, *sequentially* (each run gets the
# whole machine — concurrent runs would compete for cores and pollute timings),
# writing scaling_p<N>.json into OUTDIR. Once every rank count is done it runs
# analyze_scaling.py on OUTDIR.
#
# This is the workstation counterpart of submit_sweep.sh (which qsubs one PBS
# job per rank count on NCI Gadi).
#
# Usage:
#   ./local_sweep.sh /path/to/input.yml [ranks...]
#
# Examples:
#   ./local_sweep.sh /path/to/input.yml              # 1 2 3 4 5 6 (default)
#   ./local_sweep.sh /path/to/input.yml 1 2 4 6      # custom set
#
# Environment overrides:
#   STEPS         timesteps per run (default 20)
#   IO            on|off (default off — compute-only)
#   OUTDIR        results dir (default $PWD/scaling_results)
#   ANALYZE       1|0 run analyze_scaling.py at the end (default 1)
#   ANALYZE_OUT   analyze output dir (default $OUTDIR/summary)
#   MPIRUN        launcher (default mpirun)
#   PYTHON        python interpreter (default python)
#   OVERSUBSCRIBE 1|0 pass --oversubscribe to Open MPI when ranks > cores
#                 (auto-enabled if a rank count exceeds the detected core count)
#   MPIRUN_ARGS   extra args appended verbatim to every mpirun invocation
###############################################################################
set -euo pipefail

INPUT="${1:?usage: local_sweep.sh <input.yml> [ranks...]}"
shift || true
RANKS=("$@")
if [ "${#RANKS[@]}" -eq 0 ]; then
    RANKS=(1 2 3 4 5 6)
fi

STEPS="${STEPS:-20}"
IO="${IO:-off}"
OUTDIR="${OUTDIR:-$PWD/scaling_results}"
ANALYZE="${ANALYZE:-1}"
ANALYZE_OUT="${ANALYZE_OUT:-$OUTDIR/summary}"
MPIRUN="${MPIRUN:-mpirun}"
PYTHON="${PYTHON:-python}"
MPIRUN_ARGS="${MPIRUN_ARGS:-}"

HERE="$(cd "$(dirname "$0")" && pwd)"
DRIVER="$HERE/run_scaling.py"
ANALYZER="$HERE/analyze_scaling.py"

# Resolve the input to an absolute path: goSPL resolves the YAML's relative
# paths against the cwd, so we cd into the input's directory before each run.
INPUT_DIR="$(cd "$(dirname "$INPUT")" && pwd)"
INPUT_ABS="$INPUT_DIR/$(basename "$INPUT")"
mkdir -p "$OUTDIR"

# Detect logical cores (Linux: nproc; macOS: sysctl) to warn about / handle
# oversubscription. Open MPI 4.x refuses to launch more ranks than slots unless
# --oversubscribe is given.
if command -v nproc >/dev/null 2>&1; then
    CORES="$(nproc)"
elif command -v sysctl >/dev/null 2>&1; then
    CORES="$(sysctl -n hw.logicalcpu 2>/dev/null || echo 0)"
else
    CORES=0
fi

MAX_RANK=0
for n in "${RANKS[@]}"; do [ "$n" -gt "$MAX_RANK" ] && MAX_RANK="$n"; done

OS_FLAG=""
if [ "${OVERSUBSCRIBE:-auto}" = "1" ]; then
    OS_FLAG="--oversubscribe"
elif [ "${OVERSUBSCRIBE:-auto}" = "auto" ] && [ "$CORES" -gt 0 ] && [ "$MAX_RANK" -gt "$CORES" ]; then
    echo "[local_sweep] max rank ($MAX_RANK) > cores ($CORES) — enabling --oversubscribe"
    OS_FLAG="--oversubscribe"
fi

echo "[local_sweep] ranks=[${RANKS[*]}] steps=$STEPS io=$IO cores=$CORES -> $OUTDIR"
echo "[local_sweep] input=$INPUT_ABS"

for n in "${RANKS[@]}"; do
    echo
    echo "[local_sweep] === ranks=$n ==="
    # Run from the input's directory so the YAML's relative paths resolve, and
    # write the JSON record back into the (absolute) OUTDIR.
    ( cd "$INPUT_DIR" && \
      "$MPIRUN" -n "$n" $OS_FLAG $MPIRUN_ARGS \
          "$PYTHON" "$DRIVER" \
              -i "$INPUT_ABS" --steps "$STEPS" --io "$IO" \
              --outdir "$OUTDIR" --tag local )
done

if [ "$ANALYZE" = "1" ]; then
    echo
    echo "[local_sweep] analysing -> $ANALYZE_OUT"
    "$PYTHON" "$ANALYZER" "$OUTDIR" -o "$ANALYZE_OUT"
else
    echo
    echo "[local_sweep] done. Analyse with:"
    echo "  $PYTHON $ANALYZER $OUTDIR -o $ANALYZE_OUT"
fi
