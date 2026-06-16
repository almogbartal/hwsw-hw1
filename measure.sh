#!/bin/bash
#
# measure.sh -- perf wrappers for the Kuramoto firefly variants.
#
# Modes:
#   ./measure.sh <prog> [args...]              perf stat -r 10  (single program)
#   ./measure.sh --once <prog> [args...]       perf stat (single run)
#   ./measure.sh --stats <N>                   perf stat on every variant
#                                              at fireflies=N (single run);
#                                              writes opt_<id>_N<N>.stats into cwd
#   ./measure.sh --record <N>                  perf record -F 999 -g on every
#                                              variant at fireflies=N; writes
#                                              opt_<id>_N<N>.data and the
#                                              corresponding .html flame graph
#                                              (via FlameGraph/) into cwd
#   ./measure.sh --check                       sanity-check: run every binary
#                                              at N=200, show convergence line
#
# Program stdout is always visible in --stats / --record / --check. perf's
# own report (stat counters) goes into the .stats file; perf record writes
# its profile to the .data file.

set -e

EVENTS="cycles,instructions,task-clock,\
branch-instructions,branch-misses,\
cache-references,cache-misses,\
L1-dcache-loads,L1-dcache-load-misses,\
dTLB-loads,dTLB-load-misses,\
page-faults,context-switches"

# Binary -> label mapping for batch modes. Labels become the output filenames.
BINS=(unoptimized opt_1 opt_1_2 opt_3 opt_4 opt_5 opt_6 opt_7 opt_1_2_3)
LABELS=(opt_0     opt_1 opt_1_2 opt_3 opt_4 opt_5 opt_6 opt_7 opt_1_2_3)

usage() {
    echo "Usage: $0 [--once <prog> args...|--stats <N>|--record <N>|--check] [<prog> args...]"
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

run_stats_all() {
    local N="$1"
    [ -z "$N" ] && usage
    for k in "${!BINS[@]}"; do
        local bin="./${BINS[$k]}"
        local lab="${LABELS[$k]}_N${N}"
        if [ ! -x "$bin" ]; then
            echo "[measure.sh] skipping $bin (not built)" >&2
            continue
        fi
        echo "[measure.sh] === $lab (perf stat $bin $N) ==="
        # perf's counter report goes to stderr -> file. Program stdout
        # passes through to the terminal naturally.
        perf stat -e "$EVENTS" "$bin" "$N" 2> "${lab}.stats"
        echo "[measure.sh]   -> ${lab}.stats"
    done
}

run_record_all() {
    local N="$1"
    [ -z "$N" ] && usage

    local fg_ok=1
    if [ ! -x FlameGraph/stackcollapse-perf.pl ] || [ ! -x FlameGraph/flamegraph.pl ]; then
        echo "[measure.sh] FlameGraph/ scripts not found -- skipping .html generation" >&2
        echo "[measure.sh]   populate it with: git clone https://github.com/brendangregg/FlameGraph" >&2
        fg_ok=0
    fi

    for k in "${!BINS[@]}"; do
        local bin="./${BINS[$k]}"
        local lab="${LABELS[$k]}_N${N}"
        if [ ! -x "$bin" ]; then
            echo "[measure.sh] skipping $bin (not built)" >&2
            continue
        fi
        echo "[measure.sh] === $lab (perf record $bin $N) ==="
        # perf record writes the profile to ${lab}.data via -o. Program output
        # goes to the terminal naturally.
        perf record -F 999 -g -o "${lab}.data" -- "$bin" "$N"
        echo "[measure.sh]   -> ${lab}.data"

        if [ "$fg_ok" = "1" ]; then
            # Resolve symbols, fold stacks, render the flame graph.
            # .html extension because browsers render SVG content from .html.
            perf script -i "${lab}.data" \
                | FlameGraph/stackcollapse-perf.pl \
                | FlameGraph/flamegraph.pl --title "$lab" --subtitle "fireflies=$N" \
                > "${lab}.html"
            echo "[measure.sh]   -> ${lab}.html"
        fi
    done
}

run_stat_single() {
    local repeat="$1"; shift
    perf stat $repeat -e "$EVENTS" "$@"
}

[ $# -lt 1 ] && usage

case "$1" in
    --once)
        shift
        [ $# -lt 1 ] && usage
        run_stat_single "" "$@"
        ;;
    --stats)
        shift
        run_stats_all "$@"
        ;;
    --record)
        shift
        run_record_all "$@"
        ;;
    --check)
        run_check
        ;;
    *)
        run_stat_single "-r 10" "$@"
        ;;
esac
