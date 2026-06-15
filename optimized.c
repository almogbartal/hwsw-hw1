/*
 * firefly_sync_optimized.c -- Kuramoto firefly sync, optimized single-threaded.
 *
 * Same model, same parameters, same RNG seed as firefly_sync_unoptimized.c.
 * Optimizations vs. the baseline:
 *
 *   #1  Trig-free inner loop via the sin-of-difference identity
 *           sin(t_j - t_i) = sin(t_j)*cos(t_i) - cos(t_j)*sin(t_i)
 *       Per iteration build s[j]=sin(t_j), c[j]=cos(t_j) once (N trig calls),
 *       then form S = W*s and C = W*c. The per-firefly coupling sum becomes
 *           sum_j W_ij sin(t_j - t_i) = c_i * S_i - s_i * C_i
 *       so the O(N^2) inner loop contains no transcendental calls at all.
 *
 *   #2  Reuse those s[] and c[] tables to compute the Kuramoto order
 *       parameter r and mean phase psi (sx = sum c, sy = sum s). The baseline
 *       did a separate sin/cos pass for r; we get it for free.
 *
 *   #3  W is symmetric (distance is symmetric) and W_ii = 0, so we keep only
 *       the strict upper triangle in a packed array of length N*(N-1)/2.
 *       Half the memory, half the W reads. The matvec walks each unordered
 *       pair once and contributes to both rows: S[i]+=w*s[j], S[j]+=w*s[i].
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846

typedef struct {
    int    id;
    double x;
    double y;
    double phase;
} Firefly;

static const double OMEGA       = 1.0;
static const double K           = 1.0;
static const double SIGMA       = 0.5;
static const double DT          = 0.005;
static const double R_CONVERGE  = 0.9999;
static const long   MAX_ITERS   = 50000000L;
static const int    PRINT_EVERY = 1000;

static double urand(void) {
    return rand() / (RAND_MAX + 1.0);
}

static double wrap_2pi(double a) {
    a = fmod(a, 2.0 * PI);
    if (a < 0.0) a += 2.0 * PI;
    return a;
}

static double wrap_pi(double a) {
    a = fmod(a + PI, 2.0 * PI);
    if (a <= 0.0) a += 2.0 * PI;
    return a - PI;
}

static int cmp_double(const void *pa, const void *pb) {
    double a = *(const double *)pa, b = *(const double *)pb;
    return (a > b) - (a < b);
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: %s <num_fireflies>\n", argv[0]);
        return 1;
    }
    int N = atoi(argv[1]);
    if (N < 2) {
        fprintf(stderr, "need at least 2 fireflies\n");
        return 1;
    }

    srand(42);

    size_t Wsize = (size_t)N * (size_t)(N - 1) / 2;

    Firefly *f         = malloc((size_t)N * sizeof *f);
    double  *new_phase = malloc((size_t)N * sizeof *new_phase);
    double  *W         = malloc(Wsize     * sizeof *W);
    double  *rowsum    = malloc((size_t)N * sizeof *rowsum);
    double  *dev       = malloc((size_t)N * sizeof *dev);
    double  *s         = malloc((size_t)N * sizeof *s);
    double  *c         = malloc((size_t)N * sizeof *c);
    double  *S         = malloc((size_t)N * sizeof *S);
    double  *C         = malloc((size_t)N * sizeof *C);
    if (!f || !new_phase || !W || !rowsum || !dev || !s || !c || !S || !C) {
        fprintf(stderr, "out of memory (N too large for the O(N^2) weight cache?)\n");
        return 1;
    }

    for (int i = 0; i < N; ++i) {
        f[i].id    = i;
        f[i].x     = urand();
        f[i].y     = urand();
        f[i].phase = 2.0 * PI * urand();
    }

    /* precompute weights for unordered pairs (i,j), i<j, into the packed
       upper triangle. Accumulate rowsum for both i and j simultaneously. */
    for (int i = 0; i < N; ++i) rowsum[i] = 0.0;
    size_t idx = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dx = f[i].x - f[j].x;
            double dy = f[i].y - f[j].y;
            double d  = sqrt(dx * dx + dy * dy);
            double w  = exp(-d / SIGMA);
            W[idx++]  = w;
            rowsum[i] += w;
            rowsum[j] += w;
        }
    }
    for (int i = 0; i < N; ++i)
        if (rowsum[i] <= 0.0) rowsum[i] = 1.0;

    printf("# Kuramoto firefly synchronization (single threaded, optimized)\n");
    printf("# N=%d  OMEGA=%.3f  K=%.3f  SIGMA=%.3f  DT=%.3f  seed=42\n",
           N, OMEGA, K, SIGMA, DT);
    printf("# convergence criterion: order parameter r >= %.6f\n", R_CONVERGE);
    printf("# %-10s %-10s %-12s %-12s %-12s\n",
           "iter", "r", "mean_phase", "median_phase", "circ_std");

    long iter = 0;
    int  converged = 0;

    for (iter = 0; iter <= MAX_ITERS; ++iter) {
        /* opt #2: one sin/cos pass feeds both the order parameter and the
           dynamics update below. */
        double sx = 0.0, sy = 0.0;
        for (int i = 0; i < N; ++i) {
            double ti = f[i].phase;
            double si = sin(ti);
            double ci = cos(ti);
            s[i] = si;
            c[i] = ci;
            sx  += ci;
            sy  += si;
        }
        double r   = sqrt(sx * sx + sy * sy) / (double)N;
        double psi = atan2(sy, sx);
        if (psi < 0.0) psi += 2.0 * PI;

        if (iter % PRINT_EVERY == 0 || r >= R_CONVERGE) {
            double cstd = (r > 0.0) ? sqrt(fmax(0.0, -2.0 * log(r))) : INFINITY;

            for (int i = 0; i < N; ++i)
                dev[i] = wrap_pi(f[i].phase - psi);
            qsort(dev, (size_t)N, sizeof *dev, cmp_double);
            double med_dev = (N % 2) ? dev[N / 2]
                                     : 0.5 * (dev[N / 2 - 1] + dev[N / 2]);
            double median_phase = wrap_2pi(psi + med_dev);

            printf("  %-10ld %-10.6f %-12.6f %-12.6f %-12.6f\n",
                   iter, r, psi, median_phase, cstd);
        }

        if (r >= R_CONVERGE) { converged = 1; break; }

        /* opt #1 + #3: trig-free matvecs S = W*s, C = W*c with the symmetric
           weight stored once per pair. Each pair contributes to both rows. */
        for (int i = 0; i < N; ++i) { S[i] = 0.0; C[i] = 0.0; }
        idx = 0;
        for (int i = 0; i < N; ++i) {
            double si = s[i], ci = c[i];
            double sum_s = 0.0, sum_c = 0.0;
            for (int j = i + 1; j < N; ++j) {
                double w = W[idx++];
                sum_s += w * s[j];
                sum_c += w * c[j];
                S[j]  += w * si;
                C[j]  += w * ci;
            }
            S[i] += sum_s;
            C[i] += sum_c;
        }

        /* opt #1: sum_j W_ij sin(t_j - t_i) = c_i*S_i - s_i*C_i. */
        for (int i = 0; i < N; ++i) {
            double acc    = c[i] * S[i] - s[i] * C[i];
            double dtheta = OMEGA + K * acc / rowsum[i];
            new_phase[i]  = wrap_2pi(f[i].phase + DT * dtheta);
        }
        for (int i = 0; i < N; ++i) f[i].phase = new_phase[i];
    }

    if (converged)
        printf("# CONVERGED after %ld iterations (r >= %.6f)\n", iter, R_CONVERGE);
    else
        printf("# stopped at MAX_ITERS=%ld without reaching r >= %.6f\n",
               MAX_ITERS, R_CONVERGE);

    free(f); free(new_phase); free(W); free(rowsum); free(dev);
    free(s); free(c); free(S); free(C);
    return 0;
}
