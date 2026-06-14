#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 <program> [args...]"
    exit 1
fi

perf record -F 99 -g "$@"

perf script | npx flamegraph -o flamegraph.svg
