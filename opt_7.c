/*
 * opt_7.c -- Kuramoto firefly sync, opt #7 only (Z-order reindexing).
 *
 * Diff vs unoptimized: after generating uniformly random positions in [0,1]^2,
 * sort fireflies along a Z-order (Morton) curve. Spatially close fireflies
 * end up at close indices, so the heavy contributions of the matvec (largest
 * weights = nearest neighbors) hit cache lines that are near each other in
 * memory. Pure data layout -- the dynamics, the W storage, and the inner
 * loop structure are unchanged.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define PI 3.14159265358979323846

typedef struct {
    int    id;
    double x;
    double y;
    double phase;
} Firefly;

typedef struct {
    uint32_t morton;
    Firefly  f;
} FireflyM;

static const double OMEGA       = 1.0;
static const double K           = 1.0;
static const double SIGMA       = 0.5;
static const double DT          = 0.005;
static const double R_CONVERGE  = 0.9999;
static const long   MAX_ITERS   = 50000000L;
static const int    PRINT_EVERY = 1000;

static double urand(void) { return rand() / (RAND_MAX + 1.0); }
static double wrap_2pi(double a) { a = fmod(a, 2.0*PI); if (a < 0.0) a += 2.0*PI; return a; }
static double wrap_pi(double a)  { a = fmod(a + PI, 2.0*PI); if (a <= 0.0) a += 2.0*PI; return a - PI; }
static int cmp_double(const void *pa, const void *pb) {
    double a = *(const double *)pa, b = *(const double *)pb;
    return (a > b) - (a < b);
}

/* Spread the low 16 bits of v into the even bit positions of a 32-bit word. */
static uint32_t spread_bits(uint32_t v) {
    v &= 0x0000FFFFu;
    v = (v | (v << 8)) & 0x00FF00FFu;
    v = (v | (v << 4)) & 0x0F0F0F0Fu;
    v = (v | (v << 2)) & 0x33333333u;
    v = (v | (v << 1)) & 0x55555555u;
    return v;
}
static uint32_t morton2d(uint16_t x, uint16_t y) {
    return spread_bits(x) | (spread_bits(y) << 1);
}
static int cmp_morton(const void *a, const void *b) {
    uint32_t ma = ((const FireflyM *)a)->morton;
    uint32_t mb = ((const FireflyM *)b)->morton;
    return (ma > mb) - (ma < mb);
}

int main(int argc, char **argv) {
    if (argc != 2) { fprintf(stderr, "usage: %s <num_fireflies>\n", argv[0]); return 1; }
    int N = atoi(argv[1]);
    if (N < 2) { fprintf(stderr, "need at least 2 fireflies\n"); return 1; }
    srand(42);

    Firefly *f         = malloc((size_t)N * sizeof *f);
    double  *new_phase = malloc((size_t)N * sizeof *new_phase);
    double  *W         = malloc((size_t)N * (size_t)N * sizeof *W);
    double  *rowsum    = malloc((size_t)N * sizeof *rowsum);
    double  *dev       = malloc((size_t)N * sizeof *dev);
    FireflyM *fm       = malloc((size_t)N * sizeof *fm);
    if (!f || !new_phase || !W || !rowsum || !dev || !fm) {
        fprintf(stderr, "out of memory\n"); return 1;
    }

    /* Initialize, compute Morton code, sort by Morton, copy back into f[]. */
    for (int i = 0; i < N; ++i) {
        fm[i].f.id    = i;
        fm[i].f.x     = urand();
        fm[i].f.y     = urand();
        fm[i].f.phase = 2.0 * PI * urand();
        uint16_t xi = (uint16_t)(fm[i].f.x * 65535.0);
        uint16_t yi = (uint16_t)(fm[i].f.y * 65535.0);
        fm[i].morton  = morton2d(xi, yi);
    }
    qsort(fm, (size_t)N, sizeof *fm, cmp_morton);
    for (int i = 0; i < N; ++i) f[i] = fm[i].f;
    free(fm);

    /* Build W after reordering -- close-in-index ~ close-in-space. */
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i == j) { W[(size_t)i * N + j] = 0.0; continue; }
            double dx = f[i].x - f[j].x;
            double dy = f[i].y - f[j].y;
            double d  = sqrt(dx*dx + dy*dy);
            double w  = exp(-d / SIGMA);
            W[(size_t)i * N + j] = w;
            sum += w;
        }
        rowsum[i] = (sum > 0.0) ? sum : 1.0;
    }

    printf("# Kuramoto firefly synchronization (opt_7: Z-order reindexing)\n");
    printf("# N=%d  OMEGA=%.3f  K=%.3f  SIGMA=%.3f  DT=%.3f  seed=42\n", N, OMEGA, K, SIGMA, DT);
    printf("# convergence criterion: order parameter r >= %.6f\n", R_CONVERGE);
    printf("# %-10s %-10s %-12s %-12s %-12s\n", "iter", "r", "mean_phase", "median_phase", "circ_std");

    long iter = 0;
    int  converged = 0;
    for (iter = 0; iter <= MAX_ITERS; ++iter) {
        double sx = 0.0, sy = 0.0;
        for (int i = 0; i < N; ++i) {
            sx += cos(f[i].phase);
            sy += sin(f[i].phase);
        }
        double r   = sqrt(sx*sx + sy*sy) / (double)N;
        double psi = atan2(sy, sx);
        if (psi < 0.0) psi += 2.0 * PI;

        if (iter % PRINT_EVERY == 0 || r >= R_CONVERGE) {
            double cstd = (r > 0.0) ? sqrt(fmax(0.0, -2.0 * log(r))) : INFINITY;
            for (int i = 0; i < N; ++i) dev[i] = wrap_pi(f[i].phase - psi);
            qsort(dev, (size_t)N, sizeof *dev, cmp_double);
            double med_dev = (N % 2) ? dev[N/2] : 0.5 * (dev[N/2 - 1] + dev[N/2]);
            double median_phase = wrap_2pi(psi + med_dev);
            printf("  %-10ld %-10.6f %-12.6f %-12.6f %-12.6f\n", iter, r, psi, median_phase, cstd);
        }
        if (r >= R_CONVERGE) { converged = 1; break; }

        for (int i = 0; i < N; ++i) {
            double ti  = f[i].phase;
            double acc = 0.0;
            const double *Wi = &W[(size_t)i * N];
            for (int j = 0; j < N; ++j) {
                if (i == j) continue;
                acc += Wi[j] * sin(f[j].phase - ti);
            }
            double dtheta = OMEGA + K * acc / rowsum[i];
            new_phase[i]  = wrap_2pi(ti + DT * dtheta);
        }
        for (int i = 0; i < N; ++i) f[i].phase = new_phase[i];
    }

    if (converged) printf("# CONVERGED after %ld iterations (r >= %.6f)\n", iter, R_CONVERGE);
    else           printf("# stopped at MAX_ITERS=%ld without reaching r >= %.6f\n", MAX_ITERS, R_CONVERGE);

    free(f); free(new_phase); free(W); free(rowsum); free(dev);
    return 0;
}
