#!/bin/bash
#
# gen_table_data.sh -- extract the metrics that appear in tables.md / tables.typst
# directly from the per-variant .stats files. Run after `./measure.sh --stats <N>`.
#
# Output: one block per N (one block per `argv[1]`, or all detected N if none),
# columns = variants, rows = (Time, Speedup, IPC, Instructions(B), Br%, L1%, Cache%).
#
# Portable to macOS bash 3.2 (no associative arrays, no mapfile).

set -e

VARIANTS=(opt_0 opt_1 opt_1_2 opt_3 opt_4 opt_5 opt_6 opt_7 opt_1_2_3)
NVARS=${#VARIANTS[@]}

# Pull one numeric value from a .stats file.  $1=file, $2=grep pattern, $3=awk field.
extract() {
    local f="$1" pat="$2" col="$3"
    grep -E "$pat" "$f" 2>/dev/null | head -1 | awk -v c="$col" '{print $c}'
}

# Discover N values from existing .stats files, or take from argv.
if [ $# -gt 0 ]; then
    N_LIST=("$@")
else
    N_LIST=()
    for n in $(ls *_N*.stats 2>/dev/null \
        | sed -E 's/.*_N([0-9]+)\.stats/\1/' | sort -un); do
        N_LIST+=("$n")
    done
fi

for N in "${N_LIST[@]}"; do
    echo
    echo "=== N=$N ==="

    # Parallel indexed arrays, one slot per VARIANTS[$k].
    TIME=();  IPC=();  INSNS=();  BR=();  L1=();  CACHE=()
    base_time=""
    for k in $(seq 0 $((NVARS - 1))); do
        v="${VARIANTS[$k]}"
        f="${v}_N${N}.stats"
        if [ ! -f "$f" ]; then
            TIME[$k]="-"; IPC[$k]="-"; INSNS[$k]="-"
            BR[$k]="-";   L1[$k]="-";  CACHE[$k]="-"
            continue
        fi
        TIME[$k]=$(extract  "$f" "seconds time elapsed"  1)
        IPC[$k]=$(extract   "$f" "insn per cycle"        4)
        INSNS[$k]=$(extract "$f" "^[[:space:]]*[0-9]+[[:space:]]+instructions[[:space:]]" 1)
        BR[$k]=$(extract    "$f" "% of all branches"     4)
        L1[$k]=$(extract    "$f" "% of all L1-dcache"    4)
        CACHE[$k]=$(extract "$f" "% of all cache refs"   4)
        [ "$v" = "opt_0" ] && base_time="${TIME[$k]}"
    done

    # Header row.
    printf "%-18s" "Metric"
    for k in $(seq 0 $((NVARS - 1))); do printf "%12s" "${VARIANTS[$k]}"; done
    echo

    # Time
    printf "%-18s" "Time (s)"
    for k in $(seq 0 $((NVARS - 1))); do
        t="${TIME[$k]}"
        if [ "$t" = "-" ]; then
            printf "%12s" "-"
        else
            printf "%12s" "$(awk -v x="$t" 'BEGIN{printf "%.2f", x}')"
        fi
    done
    echo

    # Speedup vs opt_0
    printf "%-18s" "Speedup"
    for k in $(seq 0 $((NVARS - 1))); do
        t="${TIME[$k]}"
        if [ -n "$base_time" ] && [ "$t" != "-" ]; then
            printf "%12s" "$(awk -v a="$base_time" -v b="$t" 'BEGIN{printf "%.2fx", a/b}')"
        else
            printf "%12s" "-"
        fi
    done
    echo

    # IPC
    printf "%-18s" "IPC"
    for k in $(seq 0 $((NVARS - 1))); do printf "%12s" "${IPC[$k]}"; done
    echo

    # Instructions in billions
    printf "%-18s" "Instructions (B)"
    for k in $(seq 0 $((NVARS - 1))); do
        x="${INSNS[$k]}"
        if [ "$x" = "-" ]; then
            printf "%12s" "-"
        else
            printf "%12s" "$(awk -v y="$x" 'BEGIN{printf "%.1f", y/1e9}')"
        fi
    done
    echo

    # Branch miss %
    printf "%-18s" "Branch miss %"
    for k in $(seq 0 $((NVARS - 1))); do printf "%12s" "${BR[$k]}"; done
    echo

    # L1 miss %
    printf "%-18s" "L1 miss %"
    for k in $(seq 0 $((NVARS - 1))); do printf "%12s" "${L1[$k]}"; done
    echo

    # Cache miss %
    printf "%-18s" "Cache miss %"
    for k in $(seq 0 $((NVARS - 1))); do printf "%12s" "${CACHE[$k]}"; done
    echo
done
