#!/bin/bash

if [ $# -lt 1 ]; then
    exit 1
fi

if [ ! -d "FlameGraph" ]; then
    git clone https://github.com/brendangregg/FlameGraph.git > /dev/null 2>&1
fi

perf record -F 99 -e cpu-clock -g --call-graph dwarf -o perf.data -- "$@"

EXEC_NAME=$(basename "$1")

perf script -i perf.data | \
    ./FlameGraph/stackcollapse-perf.pl | \
    ./FlameGraph/flamegraph.pl --title="Flame Graph - $EXEC_NAME" > "${EXEC_NAME}_flamegraph.html"
