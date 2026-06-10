#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NUM_FIREFLIES 50000
#define DT 0.01            // Time step
#define K 0.1              // Coupling strength
#define TWO_PI 6.28318530717958647692

typedef struct {
    int id;
    double x;
    double y;
    double phase;
    double frequency;
} Firefly;

// Helper to generate random doubles between 0 and 1
double rand_double() {
    return (double)rand() / (double)RAND_MAX;
}

int main() {
    // Seed random number generator
    srand((unsigned int)time(NULL));

    // Allocate memory on the heap since 50k structs will overflow the stack
    Firefly *fireflies = (Firefly *)malloc(NUM_FIREFLIES * sizeof(Firefly));
    double *next_phases = (double *)malloc(NUM_FIREFLIES * sizeof(double));

    if (fireflies == NULL || next_phases == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return 1;
    }

    printf("Initializing %d fireflies...\n", NUM_FIREFLIES);
    for (int i = 0; i < NUM_FIREFLIES; i++) {
        fireflies[i].id = i;
        fireflies[i].x = rand_double();
        fireflies[i].y = rand_double();
        fireflies[i].phase = rand_double() * TWO_PI; // Phase between 0 and 2*pi
        // Natural frequency clustered around 1.0 Hz with some variation
        fireflies[i].frequency = 1.0 + (rand_double() - 0.5) * 0.1; 
    }

    printf("Starting simulation loop (Running 1 step as a test due to O(N^2) complexity)...\n");
    
    clock_t start_time = clock();

    // --- Single Simulation Step ---
    // Outer loop: Calculate the influence on firefly i
    for (int i = 0; i < NUM_FIREFLIES; i++) {
        double sum_influence = 0.0;
        double xi = fireflies[i].x;
        double yi = fireflies[i].y;
        double phase_i = fireflies[i].phase;

        // Inner loop: Check against all other fireflies j
        for (int j = 0; j < NUM_FIREFLIES; j++) {
            if (i == j) continue;

            // Calculate Euclidean distance
            double dx = xi - fireflies[j].x;
            double dy = yi - fireflies[j].y;
            double dist = sqrt(dx * dx + dy * dy);

            // Avoid division by zero if they overlap perfectly
            if (dist < 1e-6) dist = 1e-6; 

            // Kuramoto coupling scaled inversely by distance 
            // (closer fireflies have a stronger coupling effect)
            sum_influence += sin(fireflies[j].phase - phase_i) / dist;
        }

        // Kuramoto phase transition equation: dTheta/dt = omega + (K/N) * sum(sin(theta_j - theta_i))
        // Standardizing coupling by dividing by NUM_FIREFLIES
        double dphase_dt = fireflies[i].frequency + (K / NUM_FIREFLIES) * sum_influence;
        
        // Update phase using Euler integration
        next_phases[i] = phase_i + dphase_dt * DT;

        // Keep phase within [0, 2*pi]
        if (next_phases[i] > TWO_PI) next_phases[i] -= TWO_PI;
        if (next_phases[i] < 0.0) next_phases[i] += TWO_PI;
    }

    // Apply the synchronized phase updates
    for (int i = 0; i < NUM_FIREFLIES; i++) {
        fireflies[i].phase = next_phases[i];
    }

    clock_t end_time = clock();
    double time_taken = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("Simulation step completed successfully.\n");
    printf("Time elapsed for 1 step: %.2f seconds\n", time_taken);

    // Clean up
    free(fireflies);
    free(next_phases);

    return 0;
}
