#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define TWO_PI 6.283185307179586

typedef struct {
    double x;
    double y;
    double phase;
    double omega;
} Firefly;

static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

static double uniform01(void) {
    return (double)rand() / (double)RAND_MAX;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <num_fireflies>\n", argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);
    if (N <= 0) {
        fprintf(stderr, "Number of fireflies must be a positive integer.\n");
        return 1;
    }

    srand((unsigned)time(NULL));

    double K = 2.0;
    double dt = 0.01;
    double sigma = 0.3;
    double mean_omega = 1.0;
    double omega_spread = 0.1;
    int max_iter = 100000;
    int print_every = 1000;
    double sync_threshold = 0.99;
    int sync_stable_required = 3;
    int sync_stable_count = 0;

    Firefly *fireflies = (Firefly *)malloc((size_t)N * sizeof(Firefly));
    double *new_phases = (double *)malloc((size_t)N * sizeof(double));
    double *phase_buf = (double *)malloc((size_t)N * sizeof(double));
    if (!fireflies || !new_phases || !phase_buf) {
        fprintf(stderr, "Allocation failed.\n");
        return 1;
    }

    for (int i = 0; i < N; i++) {
        fireflies[i].x = uniform01();
        fireflies[i].y = uniform01();
        fireflies[i].phase = uniform01() * TWO_PI;
        fireflies[i].omega = mean_omega + (uniform01() - 0.5) * 2.0 * omega_spread;
    }

    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < N; i++) {
            double coupling_sum = 0.0;
            for (int j = 0; j < N; j++) {
                if (i == j) continue;
                double dx = fireflies[i].x - fireflies[j].x;
                double dy = fireflies[i].y - fireflies[j].y;
                double dist = sqrt(dx * dx + dy * dy);
                double weight = exp(-dist / sigma);
                coupling_sum += weight * sin(fireflies[j].phase - fireflies[i].phase);
            }
            double dtheta = fireflies[i].omega + (K / (double)N) * coupling_sum;
            double next = fireflies[i].phase + dt * dtheta;
            next = fmod(next, TWO_PI);
            if (next < 0.0) next += TWO_PI;
            new_phases[i] = next;
        }

        for (int i = 0; i < N; i++) {
            fireflies[i].phase = new_phases[i];
        }

        if ((iter + 1) % print_every == 0) {
            double sum_cos = 0.0;
            double sum_sin = 0.0;
            for (int i = 0; i < N; i++) {
                sum_cos += cos(fireflies[i].phase);
                sum_sin += sin(fireflies[i].phase);
            }
            sum_cos /= (double)N;
            sum_sin /= (double)N;

            double r = sqrt(sum_cos * sum_cos + sum_sin * sum_sin);
            double psi = atan2(sum_sin, sum_cos);
            if (psi < 0.0) psi += TWO_PI;

            double circ_std = (r > 0.0) ? sqrt(-2.0 * log(r)) : INFINITY;

            for (int i = 0; i < N; i++) {
                double d = fireflies[i].phase - psi;
                while (d < -M_PI) d += TWO_PI;
                while (d >= M_PI) d -= TWO_PI;
                phase_buf[i] = d;
            }
            qsort(phase_buf, (size_t)N, sizeof(double), cmp_double);
            double med_centered;
            if (N % 2 == 0) {
                med_centered = 0.5 * (phase_buf[N / 2 - 1] + phase_buf[N / 2]);
            } else {
                med_centered = phase_buf[N / 2];
            }
            double median = med_centered + psi;
            median = fmod(median, TWO_PI);
            if (median < 0.0) median += TWO_PI;

            if (r > sync_threshold) {
                sync_stable_count++;
            } else {
                sync_stable_count = 0;
            }
            const char *flag = (sync_stable_count > 0) ? "  [SYNCHRONIZED]" : "";
            printf("iter %6d  r=%.4f  mean=%.4f  median=%.4f  std=%.4f%s\n",
                   iter + 1, r, psi, median, circ_std, flag);

            if (sync_stable_count >= sync_stable_required) {
                printf("Converged: r > %.2f for %d consecutive checks (%d iterations). Stopping.\n",
                       sync_threshold, sync_stable_required,
                       sync_stable_required * print_every);
                break;
            }
        }
    }

    free(fireflies);
    free(new_phases);
    free(phase_buf);
    return 0;
}
