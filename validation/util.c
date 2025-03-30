#include "validation.h"

// Global constants
double pi = 3.14159265358979323846;
long double E = 2.7182818284590452353602874;
double errMax = 0;
int Dim = 0;
int bcsNb = 0;
int btsNb = 0;

// Simple random number generator (0 to 1)
double rand_double() {
    return (double)rand() / (double)RAND_MAX;
}

// Generate a random number between a and b
// This is a simplified version that only uses standard C random numbers
double alea(double a, double b, int randCase) {
    // Ignore randCase parameter in this implementation for simplicity
    // We're using only the standard C rand() function
    return a + (b - a) * rand_double();
}

// Normal distribution using Box-Muller transform
double alea_normal(double mean, double std_dev, int option) {
    double x1, x2, w, y1;
    
    do {
        x1 = 2.0 * rand_double() - 1.0;
        x2 = 2.0 * rand_double() - 1.0;
        w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);
    
    w = sqrt(-2.0 * log(w) / w);
    y1 = x1 * w;
    
    if (rand_double() < 0.5) y1 = -y1;
    y1 = y1 * std_dev + mean;
    return y1;
}

// Sign function
double sign(double a) {
    if (a > 0) return 1.0;
    if (a < 0) return -1.0;
    return 0.0;
}

// L-distance between two positions
double distanceL(struct position x, struct position y, int L) {
    int k;
    double d, s = 0;
    int n = x.size < y.size ? x.size : y.size; // Min dimension

    if (L == 0) {
        // L == 0 means L-infinity: maximum difference
        d = 0;
        for (k = 0; k < n; k++) {
            s = fabs(x.x[k] - y.x[k]);
            if (s > d) d = s;
        }
        return d;
    }

    for (k = 0; k < n; k++) {
        s = s + pow(fabs(x.x[k] - y.x[k]), L);
    }
    
    return pow(s, 1.0 / L);
} 