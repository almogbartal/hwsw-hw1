#!/bin/bash
#
# measure.sh -- perf wrappers for the Kuramoto firefly variants.
#
# Modes:
#   ./measure.sh <prog> [args...]              perf stat -r 10  (default)
#   ./measure.sh --once <prog> [args...]       perf stat (single run)
#   ./measure.sh --record <prog> [args...]     perf record (writes perf.data)
#   ./measure.sh --all <N>                     run every binary at fireflies=N;
#                                              writes opt_<id>_N<N>.stats and
#                                              opt_<id>_N<N>.data into cwd
#   ./measure.sh --check                       sanity-check: run every binary
#                                              at N=200, show convergence line
#
# Event list is shared across modes so numbers are comparable.

set -e

EVENTS="cycles,instructions,task-clock,\
branch-instructions,branch-misses,\
L1-dcache-load-misses,LLC-load-misses,dTLB-load-misses,\
stalled-cycles-frontend,stalled-cycles-backend,\
page-faults,context-switches"

# Binary -> label mapping for --all mode. Labels become the output filenames.
BINS=(unoptimized opt_1 opt_1_2 opt_3 opt_4 opt_5 opt_6 opt_7 opt_1_2_3)
LABELS=(opt_0     opt_1 opt_1_2 opt_3 opt_4 opt_5 opt_6 opt_7 opt_1_2_3)

usage() {
    echo "Usage: $0 [--once|--record|--all|--check] <program|N> [args...]"
    exit 1
}

run_check() {
    local N=200
    for k in "${!BINS[@]}"; do
        local bin="./${BINS[$k]}"
        local lab="${LABELS[$k]}"
        if [ ! -x "$bin" ]; then
            echo "[measure.sh] skipping $bin (not built)" >&2
            continue
        fi
        echo "=== $lab ($bin $N) ==="
        "$bin" "$N" | tail -2
    done
}

run_stat() {
    local repeat="$1"; shift
    perf stat $repeat -e "$EVENTS" "$@"
}

run_record() {
    perf record -F 999 -g -- "$@"
}

run_all() {
    local N="$1"
    [ -z "$N" ] && usage

    for k in "${!BINS[@]}"; do
        local bin="./${BINS[$k]}"
        local lab="${LABELS[$k]}_N${N}"
        if [ ! -x "$bin" ]; then
            echo "[measure.sh] skipping $bin (not built)" >&2
            continue
        fi
        echo "[measure.sh] === $lab ($bin $N) ==="
        # stats: redirect perf's output (stderr) to the stats file
        perf stat -r 10 -e "$EVENTS" "$bin" "$N" 2> "${lab}.stats" > /dev/null
        # record: explicit -o so we control the filename
        perf record -F 999 -g -o "${lab}.data" -- "$bin" "$N" > /dev/null 2>&1
        echo "[measure.sh]   -> ${lab}.stats, ${lab}.data"
    done
    echo "[measure.sh] done. inspect with: perf report -i opt_<label>_N${N}.data"
}

[ $# -lt 1 ] && usage

case "$1" in
    --once)
        shift
        [ $# -lt 1 ] && usage
        run_stat "" "$@"
        ;;
    --record)
        shift
        [ $# -lt 1 ] && usage
        run_record "$@"
        ;;
    --all)
        shift
        run_all "$@"
        ;;
    --check)
        run_check
        ;;
    *)
        run_stat "-r 10" "$@"
        ;;
esac
