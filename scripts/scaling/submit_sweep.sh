#!/bin/bash
###############################################################################
# Submit a goSPL strong-scaling sweep on NCI Gadi (PBS).
#
# qsub's gadi.pbs once per rank count, overriding ncpus/mem and passing
# the run parameters via -v. Each job writes scaling_p<N>.json into OUTDIR; once
# all jobs finish, run analyze_scaling.py on OUTDIR.
#
# Usage:
#   ./submit_sweep.sh /path/to/input.yml [ranks...]
#
# Examples:
#   ./submit_sweep.sh /scratch/q97/me/run/input_30km.yml                 # 1 2 4 8 16 24 48
#   ./submit_sweep.sh /scratch/q97/me/run/input_15km.yml 2 4 8 16 96     # custom set, multi-node ok
#   RES=10km ./submit_sweep.sh /scratch/q97/me/run/input.yml 48 96 192   # explicit resolution tag
#
#
# WALLTIME=02:00:00 \
# WALLTIME_MAP="4=03:00:00,8=02:30:00,96=01:30:00,144=01:00:00,192=00:45:00" \
# PROJECT=do20 GOSPL_VENV=$HOME/gospl-smoke/bin/activate RES=10km \
#     ./submit_sweep.sh input_10km.yml 4 8 16 24 48 96 144 192
#
#
# Environment overrides:
#   RES     resolution tag used to label OUTDIR/job names (default: parsed
#           from the input filename, e.g. input_30km.yml -> 30km; falls back
#           to "run" if nothing matches \d+km)
#   STEPS   timesteps per run (default 20)
#   IO      on|off (default off — compute-only)
#   MODE    native|container (default native)
#   OUTDIR  results dir (default $PWD/scaling_<RES>)
#   MEM_PER_CPU_GB  memory request per rank (default 4)
#   WALLTIME  per-job walltime (default 01:00:00), applied to any rank count
#             not covered by WALLTIME_MAP
#   WALLTIME_MAP  optional per-rank-count walltime overrides, as a comma-
#             separated list of N=hh:mm:ss, e.g.:
#               WALLTIME_MAP="4=02:00:00,8=01:30:00,96=00:20:00,192=00:15:00"
#             Rank counts not listed fall back to WALLTIME.
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

# Resolution tag: explicit RES wins; otherwise sniff it out of the input
# filename (input_30km.yml, mesh-15km.yml, ... -> 30km / 15km); else "run".
if [ -n "${RES:-}" ]; then
    RES_TAG="$RES"
else
    base="$(basename "$INPUT")"
    if [[ "$base" =~ ([0-9]+km) ]]; then
        RES_TAG="${BASH_REMATCH[1]}"
    else
        RES_TAG="run"
    fi
fi

STEPS="${STEPS:-10}"
IO="${IO:-off}"
MODE="${MODE:-native}"
OUTDIR="${OUTDIR:-$PWD/scaling_${RES_TAG}}"
MEM_PER_CPU_GB="${MEM_PER_CPU_GB:-4}"
WALLTIME="${WALLTIME:-01:00:00}"
QUEUE="${QUEUE:-normal}"
HERE="$(cd "$(dirname "$0")" && pwd)"
PBS="${PBS:-$HERE/gadi.pbs}"

# Optional per-rank-count walltime overrides: WALLTIME_MAP="4=02:00:00,96=00:20:00"
declare -A WALLTIME_BY_N=()
if [ -n "${WALLTIME_MAP:-}" ]; then
    IFS=',' read -r -a _wt_pairs <<< "$WALLTIME_MAP"
    for pair in "${_wt_pairs[@]}"; do
        key="${pair%%=*}"
        val="${pair#*=}"
        if [ -z "$key" ] || [ -z "$val" ] || [ "$key" = "$pair" ]; then
            echo "ERROR: malformed WALLTIME_MAP entry '$pair' (expected N=hh:mm:ss)" >&2
            exit 1
        fi
        WALLTIME_BY_N["$key"]="$val"
    done
fi

mkdir -p "$OUTDIR"
INPUT_ABS="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"

echo "Sweep: res=$RES_TAG ranks=[${RANKS[*]}] steps=$STEPS io=$IO mode=$MODE -> $OUTDIR"

CORES_PER_NODE="${CORES_PER_NODE:-48}"

prev_job=""
for n in "${RANKS[@]}"; do
    mem=$(( n * MEM_PER_CPU_GB ))
    # NOTE: passed as RANKS, not NCPUS — NCI/Gadi auto-exports a reserved
    # NCPUS env var into every job equal to the node's physical core count
    # (e.g. 48), which silently shadows any qsub -v NCPUS=<n> override for
    # n > that. gadi.pbs reads RANKS to avoid the collision.
    vargs="RANKS=$n,INPUT=$INPUT_ABS,STEPS=$STEPS,IO=$IO,MODE=$MODE,OUTDIR=$OUTDIR"
    [ -n "${CONTAINER:-}" ] && vargs="$vargs,CONTAINER=$CONTAINER"
    [ -n "${GOSPL_VENV:-}" ] && vargs="$vargs,GOSPL_VENV=$GOSPL_VENV"

    # Gadi's qsub does NOT support PBS Pro's "-l select=..." syntax
    # ("Illegal attribute or resource value" / "doesn't support -l select").
    # Multi-node jobs here just use a larger flat ncpus=/mem= request — Gadi's
    # scheduler spans nodes for you, but only in whole-node multiples above
    # CORES_PER_NODE (it rejects e.g. ncpus=80 outright; 96 = 2 whole nodes
    # is fine). Mem is requested as total (not per-node).
    if [ "$n" -gt "$CORES_PER_NODE" ] && [ $(( n % CORES_PER_NODE )) -ne 0 ]; then
        echo "ERROR: ranks=$n is not a whole multiple of CORES_PER_NODE=$CORES_PER_NODE" >&2
        echo "       Gadi rejects multi-node requests that don't use whole nodes." >&2
        exit 1
    fi
    l_args=(-l "ncpus=$n" -l "mem=${mem}GB")

    job_walltime="${WALLTIME_BY_N[$n]:-$WALLTIME}"

    qsub_args=(-q "$QUEUE" "${l_args[@]}"
               -l "walltime=$job_walltime" -v "$vargs" -N "gospl_${RES_TAG}_p${n}")
    [ -n "${PROJECT:-}" ] && qsub_args+=(-P "$PROJECT")
    if [ "${DEPEND:-0}" = "1" ] && [ -n "$prev_job" ]; then
        qsub_args+=(-W "depend=afterany:$prev_job")
    fi

    jobid=$(qsub "${qsub_args[@]}" "$PBS")
    echo "  ranks=$n  mem=${mem}GB  walltime=$job_walltime  job=$jobid"
    prev_job="$jobid"
done

echo
echo "When all jobs finish:"
echo "  python $HERE/analyze_scaling.py $OUTDIR -o results/scaling_${RES_TAG}"