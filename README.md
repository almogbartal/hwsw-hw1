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
