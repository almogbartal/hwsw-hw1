#!/bin/bash

repeat="-r 10"
if [ "$1" = "--once" ]; then
    repeat=""
    shift
fi

if [ $# -lt 1 ]; then
    echo "Usage: $0 [--once] <program> [args...]"
    exit 1
fi

# "$@" captures all arguments passed to the script
perf stat $repeat -e branch-instructions,branch-misses,cache-references,cache-misses,cpu-cycles,cpu-clock,page-faults "$@"
