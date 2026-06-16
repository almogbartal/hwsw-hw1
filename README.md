# Prompt

generate me a kuramoto model firefly synchronization simulation in c that takes as an input parameter the number of fireflies to synchronize. every firefly should be specified by 2 floating point coordinates and they live on a [0,1]x[0,1] square. the fireflies should be placed randomly and uniformly inside the square, and the distance should be euclidean distance. different parameters, like coupling strength etc. can be defined arbitrarily or by looking at other projects. on the screen you should print every 1000 iterations mean phase, median phase and std. all fireflies should be in the same frequency.

write it in a single .c file. 

what can be a criterion to see that the simulation converges to synchronization? what can be printed on the cli to see this without visually rendering all the fireflies? implement it and print this parameter on the screen beside all the other parameters

also the firefly phases should also be initialized as uniformly random on [0:2pi]
make it stop whenever the simulation has converged. initialize the random number generator with value 42.

for every firefly make a struct with id, x,y and phase.

do it single threaded.

to change fireflies phase you should compare each firefly against all the rest


## Reproduce

```
make clean && make all
./measure.sh --stats 1000     # writes opt_<id>_N1000.stats (single run per variant)
./measure.sh --record 1000    # writes opt_<id>_N1000.data
```

Produces `opt_0_N1000.{stats,data}` ... `opt_7_N1000.{stats,data}` plus
`opt_1_2_3_N1000.{stats,data}` in the current directory. The `.stats` files
hold `perf stat -r 10` output, the `.data` files are `perf record -F 999 -g`
captures suitable for `perf report` / `perf annotate`. Re-running at a
different N (e.g. `./measure.sh --all 2000`) produces a parallel set of
`_N2000.{stats,data}` files instead of overwriting.


`measure.sh` flags:
- `./measure.sh <prog> [args]` -- `perf stat -r 10` on a single program
- `./measure.sh --once <prog> [args]` -- single `perf stat` run on one program
- `./measure.sh --stats <N>` -- `perf stat` on every variant at fireflies=N (single run each)
- `./measure.sh --record <N>` -- `perf record -F 999 -g` on every variant at fireflies=N, plus an `.html` flame graph rendered from each `.data` file (requires `FlameGraph/` populated -- `git clone https://github.com/brendangregg/FlameGraph`)
- `./measure.sh --check` -- quick correctness pass: run each variant at N=200, show the convergence line

Event list (shared by all modes):
```
cycles, instructions, task-clock,
branch-instructions, branch-misses,
cache-references, cache-misses,
L1-dcache-loads, L1-dcache-load-misses,
dTLB-loads, dTLB-load-misses,
page-faults, context-switches
```
