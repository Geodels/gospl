#!/bin/bash
###############################################################################
# Submit a goSPL strong-scaling sweep on NCI Gadi (PBS).
#
# qsub's gadi_scaling.pbs once per rank count, overriding ncpus/mem and passing
# the run parameters via -v. Each job writes scaling_p<N>.json into OUTDIR; once
# all jobs finish, run analyze_scaling.py on OUTDIR.
#
# Usage:
#   ./submit_sweep.sh /path/to/input.yml [ranks...]
#
# Examples:
#   ./submit_sweep.sh /scratch/q97/me/run/input.yml                 # 1 2 4 8 16 24 48
#   ./submit_sweep.sh /scratch/q97/me/run/input.yml 2 4 8 16        # custom set
#
# Environment overrides:
#   STEPS   timesteps per run (default 20)
#   IO      on|off (default off — compute-only)
#   MODE    native|container (default native)
#   OUTDIR  results dir (default $PWD/scaling_results)
#   MEM_PER_CPU_GB  memory request per rank (default 4)
#   WALLTIME  per-job walltime (default 00:30:00)
#   CONTAINER  .sif path (MODE=container)
#   PROJECT  NCI project code for -P (default: from the .pbs file)
#   QUEUE    PBS queue (default normal)
#   DEPEND   if set to 1, chain jobs so only one runs at a time (cleaner timing)
###############################################################################
set -euo pipefail

INPUT="${1:?usage: submit_sweep.sh <input.yml> [ranks...]}"
shift || true
RANKS=("${@:-1 2 4 8 16 24 48}")
# When no ranks were passed, the default above is a single word; re-split it.
if [ "${#RANKS[@]}" -eq 1 ]; then
    read -r -a RANKS <<< "${RANKS[0]}"
fi

STEPS="${STEPS:-20}"
IO="${IO:-off}"
MODE="${MODE:-native}"
OUTDIR="${OUTDIR:-$PWD/scaling_results}"
MEM_PER_CPU_GB="${MEM_PER_CPU_GB:-4}"
WALLTIME="${WALLTIME:-00:30:00}"
QUEUE="${QUEUE:-normal}"
HERE="$(cd "$(dirname "$0")" && pwd)"
PBS="$HERE/gadi_scaling.pbs"

mkdir -p "$OUTDIR"
INPUT_ABS="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"

echo "Sweep: ranks=[${RANKS[*]}] steps=$STEPS io=$IO mode=$MODE -> $OUTDIR"

prev_job=""
for n in "${RANKS[@]}"; do
    mem=$(( n * MEM_PER_CPU_GB ))
    vargs="NCPUS=$n,INPUT=$INPUT_ABS,STEPS=$STEPS,IO=$IO,MODE=$MODE,OUTDIR=$OUTDIR"
    [ -n "${CONTAINER:-}" ] && vargs="$vargs,CONTAINER=$CONTAINER"
    [ -n "${GOSPL_VENV:-}" ] && vargs="$vargs,GOSPL_VENV=$GOSPL_VENV"

    qsub_args=(-q "$QUEUE" -l "ncpus=$n" -l "mem=${mem}GB"
               -l "walltime=$WALLTIME" -v "$vargs" -N "gospl_p${n}")
    [ -n "${PROJECT:-}" ] && qsub_args+=(-P "$PROJECT")
    if [ "${DEPEND:-0}" = "1" ] && [ -n "$prev_job" ]; then
        qsub_args+=(-W "depend=afterany:$prev_job")
    fi

    jobid=$(qsub "${qsub_args[@]}" "$PBS")
    echo "  ranks=$n  mem=${mem}GB  job=$jobid"
    prev_job="$jobid"
done

echo
echo "When all jobs finish:"
echo "  python $HERE/analyze_scaling.py $OUTDIR -o results/scaling"
