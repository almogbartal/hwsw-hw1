#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 <program> [args...]"
    exit 1
fi

# "$@" captures all arguments passed to the script
perf stat -r 10 -e branch-instructions,branch-misses,cache-references,cache-misses,cpu-cycles,cpu-clock,page-faults "$@"
