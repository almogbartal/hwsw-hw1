# Prompt

generate me a kuramoto model firefly synchronization simulation in c that takes as an input parameter the number of fireflies to synchronize. every firefly should be specified by 2 floating point coordinates and they live on a [0,1]x[0,1] square. the fireflies should be placed randomly and uniformly inside the square, and the distance should be euclidean distance. different parameters, like coupling strength etc. can be defined arbitrarily or by looking at other projects. on the screen you should print every 1000 iterations mean phase, median phase and std. all fireflies should be in the same frequency.

write it in a single .c file. 

what can be a criterion to see that the simulation converges to synchronization? what can be printed on the cli to see this without visually rendering all the fireflies? implement it and print this parameter on the screen beside all the other parameters

also the firefly phases should also be initialized as uniformly random on [0:2pi]
make it stop whenever the simulation has converged. initialize the random number generator with value 42.

for every firefly make a struct with id, x,y and phase.

do it single threaded.

to change fireflies phase you should compare each firefly against all the rest

# Optimizations

`firefly_sync_optimized.c` is a drop-in replacement for the baseline with the
same parameters, RNG seed, and convergence criterion. Three optimizations are
applied; together they eliminate all transcendental calls from the O(N^2) hot
loop, halve the weight-matrix memory, and remove a redundant trig pass per
iteration.

## 1. Trig-free inner loop via sin-of-difference identity

The baseline computes `sin(theta_j - theta_i)` inside the inner `j` loop, so
each iteration does N^2 `sin` calls.

Using the identity

    sin(t_j - t_i) = sin(t_j) * cos(t_i) - cos(t_j) * sin(t_i)

the coupling sum factors as

    sum_j W_ij * sin(t_j - t_i)
        = cos(t_i) * sum_j W_ij * sin(t_j)  -  sin(t_i) * sum_j W_ij * cos(t_j)
        = c_i * S_i  -  s_i * C_i

where `s_j = sin(t_j)`, `c_j = cos(t_j)` are computed once per iteration (N
trig calls total) and `S = W*s`, `C = W*c` are two matrix-vector products that
contain only multiplies and adds. The hot O(N^2) inner loop becomes
trig-free, which both removes the cost of the libm calls and exposes the
loop to straightforward SIMD/FMA codegen.

## 2. Reuse `sin`/`cos` for the order parameter

The Kuramoto order parameter needs `sx = sum_j cos(t_j)` and
`sy = sum_j sin(t_j)`. The baseline does this in its own loop, separate from
the dynamics update. Since optimization #1 already builds the same `s[j]`
and `c[j]` tables, we just accumulate `sx` and `sy` during that pass and drop
the separate trig loop entirely. Saves N `sin` and N `cos` calls per
iteration with no extra memory.

## 3. Symmetric W: store only the strict upper triangle

The weight matrix is symmetric (`W_ij = W_ji`, distance is symmetric) and the
diagonal is zero (no self-coupling). The baseline stores all N^2 entries; we
store only the `N*(N-1)/2` strict upper-triangle entries in a packed array.
This halves the memory footprint of `W` (the dominant allocation) and halves
the number of W reads per iteration.

The matvec is rewritten as a pair walk: for each unordered pair `(i,j)` with
`i < j`, read `w = W_ij` once and contribute to both rows

    S[i] += w * s[j];   S[j] += w * s[i];
    C[i] += w * c[j];   C[j] += w * c[i];

Row sums (`rowsum[i]`) are computed during precomputation by the same pair
walk, so the symmetric storage extends cleanly to initialization.

# Ablation

Each ablation variant changes **one thing** vs. the baseline so each
optimization can be measured in isolation.

| Label   | File           | Optimization isolated                                | Key perf event to watch                |
|---------|----------------|------------------------------------------------------|----------------------------------------|
| opt_0   | `unoptimized.c`| (none -- control)                                    | --                                     |
| opt_1   | `opt_1.c`      | sin-of-difference identity (no trig in inner loop)   | `instructions`, `cycles` (IPC)         |
| opt_1_2 | `opt_1_2.c`    | opt_1 + reuse `s[]`/`c[]` for the order parameter    | `instructions` vs opt_1 (drops 2N)     |
| opt_3   | `opt_3.c`      | packed upper-triangle W + pair-walk matvec           | `LLC-load-misses` at large N           |
| opt_4   | `opt_4.c`      | drop `if (i == j) continue` from the inner loop      | `instructions`, IPC (vectorization)    |
| opt_5   | `opt_5.c`      | store W as `float` (half the W footprint)            | `L1-dcache-load-misses`, `LLC` at mid N |
| opt_6   | `opt_6.c`      | Structure-of-Arrays firefly layout                   | `L1-dcache-load-misses`, IPC           |
| opt_7   | `opt_7.c`      | Z-order (Morton) reindexing of fireflies             | `L1-dcache-load-misses`, `LLC`         |
| opt_1_2_3 | `opt_1_2_3.c` | #1 + #2 + #3 combined                              | every event vs opt_0                   |

Notes:
- opt_1_2 isn't "opt_2 alone" because optimization #2 (reusing `s[]`/`c[]`
  for the order parameter) only makes sense on top of #1, which builds those
  tables in the first place.
- opt_3, by switching to an `i<j` pair walk, naturally has no `i==j` branch.
  That's intrinsic to the storage layout, not an extra optimization.
- opt_5's dynamics differ by a tiny amount of float rounding per W entry;
  the convergence trajectory stays well within the 1e-4 threshold.
- opt_6 (SoA) and opt_7 (Z-order) leave the math identical to baseline --
  they are pure data-layout changes targeting L1.

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

To see flamegraph-style hotspots for a specific variant:
```
perf report -i opt_1_2_3_N1000.data
perf annotate -i opt_1_2_3_N1000.data
```

`measure.sh` flags:
- `./measure.sh <prog> [args]` -- `perf stat -r 10` on a single program
- `./measure.sh --once <prog> [args]` -- single `perf stat` run on one program
- `./measure.sh --stats <N>` -- `perf stat` on every variant at fireflies=N (single run each)
- `./measure.sh --record <N>` -- `perf record -F 999 -g` on every variant at fireflies=N
- `./measure.sh --check` -- quick correctness pass: run each variant at N=200, show the convergence line

Event list (shared by all modes):
```
cycles, instructions, task-clock,
branch-instructions, branch-misses,
L1-dcache-load-misses, LLC-load-misses, dTLB-load-misses,
stalled-cycles-frontend, stalled-cycles-backend,
page-faults, context-switches
```
