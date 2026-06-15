/*
 * fireflies.c  --  Kuramoto-model firefly synchronization (single threaded).
 *
 * Build:  gcc -O2 -o fireflies fireflies.c -lm
 * Run:    ./fireflies <num_fireflies>
 *
 * Model
 * -----
 * N fireflies live on the unit square [0,1]x[0,1] at random uniform positions.
 * Each carries a phase theta in [0,2*pi) (random uniform at t=0) and shares a
 * single common natural frequency OMEGA. They are coupled by the Kuramoto rule
 * with a distance-dependent weight (Euclidean distance), so nearby fireflies
 * pull on each other more strongly than far ones:
 *
 *     w_ij = exp(-d_ij / SIGMA),        d_ij = || (x_i,y_i) - (x_j,y_j) ||_2
 *
 *     dtheta_i/dt = OMEGA + K * ( sum_{j!=i} w_ij * sin(theta_j - theta_i) )
 *                                / ( sum_{j!=i} w_ij )
 *
 * Updated with forward Euler. Because every firefly compares itself against all
 * the rest, the inner loop is O(N^2) per iteration. The weights depend only on
 * the (static) positions, so they are computed once and cached.
 *
 * Convergence criterion (printed to the CLI, no rendering needed)
 * --------------------------------------------------------------
 * The Kuramoto ORDER PARAMETER:
 *
 *     r * e^{i*psi} = (1/N) * sum_j e^{i*theta_j}
 *
 * r lies in [0,1]: r ~ 0 means phases are scattered (incoherent), r -> 1 means
 * all phases coincide (fully synchronized), and psi is the common (mean) phase.
 * r is invariant to the shared OMEGA drift, so it cleanly measures alignment.
 * The run stops once r >= R_CONVERGE. Equivalently the circular standard
 * deviation sqrt(-2 ln r) -> 0 at synchronization; it is printed too.
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

/* ---- parameters (chosen as typical values for this kind of simulation) ---- */
static const double OMEGA       = 1.0;     /* shared natural frequency (rad/t) */
static const double K           = 1.0;     /* coupling strength                */
static const double SIGMA       = 0.5;     /* spatial coupling length scale    */
static const double DT          = 0.005;   /* Euler integration time step      */
static const double R_CONVERGE  = 0.9999;  /* order-parameter sync threshold   */
static const long   MAX_ITERS   = 50000000L;
static const int    PRINT_EVERY = 1000;

/* uniform double in [0,1) */
static double urand(void) {
    return rand() / (RAND_MAX + 1.0);
}

/* wrap angle to [0, 2*PI) */
static double wrap_2pi(double a) {
    a = fmod(a, 2.0 * PI);
    if (a < 0.0) a += 2.0 * PI;
    return a;
}

/* wrap angle to (-PI, PI] */
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

    /* deterministic run: seed RNG with 42 */
    srand(42);

    Firefly *f         = malloc((size_t)N * sizeof *f);
    double  *new_phase = malloc((size_t)N * sizeof *new_phase);
    double  *W         = malloc((size_t)N * (size_t)N * sizeof *W); /* static weights */
    double  *rowsum    = malloc((size_t)N * sizeof *rowsum);
    double  *dev       = malloc((size_t)N * sizeof *dev);
    if (!f || !new_phase || !W || !rowsum || !dev) {
        fprintf(stderr, "out of memory (N too large for the O(N^2) weight cache?)\n");
        return 1;
    }

    /* initialize: uniform random position in the unit square, uniform random
       phase in [0, 2*pi) */
    for (int i = 0; i < N; ++i) {
        f[i].id    = i;
        f[i].x     = urand();
        f[i].y     = urand();
        f[i].phase = 2.0 * PI * urand();
    }

    /* precompute distance-based coupling weights once (positions never change).
       w_ij = exp(-euclidean_distance / SIGMA), self-coupling excluded. */
    for (int i = 0; i < N; ++i) {
        double s = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i == j) { W[(size_t)i * N + j] = 0.0; continue; }
            double dx = f[i].x - f[j].x;
            double dy = f[i].y - f[j].y;
            double d  = sqrt(dx * dx + dy * dy);   /* Euclidean distance */
            double w  = exp(-d / SIGMA);
            W[(size_t)i * N + j] = w;
            s += w;
        }
        rowsum[i] = (s > 0.0) ? s : 1.0;
    }

    printf("# Kuramoto firefly synchronization (single threaded)\n");
    printf("# N=%d  OMEGA=%.3f  K=%.3f  SIGMA=%.3f  DT=%.3f  seed=42\n",
           N, OMEGA, K, SIGMA, DT);
    printf("# convergence criterion: order parameter r >= %.6f\n", R_CONVERGE);
    printf("# %-10s %-10s %-12s %-12s %-12s\n",
           "iter", "r", "mean_phase", "median_phase", "circ_std");

    long iter = 0;
    int  converged = 0;

    for (iter = 0; iter <= MAX_ITERS; ++iter) {
        /* order parameter r and mean phase psi from the current phases */
        double sx = 0.0, sy = 0.0;
        for (int i = 0; i < N; ++i) {
            sx += cos(f[i].phase);
            sy += sin(f[i].phase);
        }
        double r   = sqrt(sx * sx + sy * sy) / (double)N;
        double psi = atan2(sy, sx);
        if (psi < 0.0) psi += 2.0 * PI;

        if (iter % PRINT_EVERY == 0 || r >= R_CONVERGE) {
            /* circular standard deviation: sqrt(-2 ln r) -> 0 as r -> 1 */
            double cstd = (r > 0.0) ? sqrt(fmax(0.0, -2.0 * log(r))) : INFINITY;

            /* circular median: median of deviations from the mean phase,
               folded back onto the circle */
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

        /* synchronous Kuramoto update: each firefly vs. all the rest (O(N^2)) */
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

    if (converged)
        printf("# CONVERGED after %ld iterations (r >= %.6f)\n", iter, R_CONVERGE);
    else
        printf("# stopped at MAX_ITERS=%ld without reaching r >= %.6f\n",
               MAX_ITERS, R_CONVERGE);

    free(f); free(new_phase); free(W); free(rowsum); free(dev);
    return 0;
}
