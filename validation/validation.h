#ifndef VALIDATION_H
#define VALIDATION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

// Maximum dimensions
#define DMax 100

// Maximum number of functions
#define funcMax 30

// Common constants
#define zero 1.e-20
#define infinity 1.e20

# define randNbMax 500 // // When using a list of previously generated random numbers

// Global constants
extern double pi;
extern long double E;
extern double errMax;
extern double errMin;
extern int Dim;
extern int bcsNb;
extern int btsNb;

int nCycleMax; //=infinity;		// Length of the cycle, when reading the random numbers on a file
// (BW[2]=3 or 4) Default: infinity
int nCycle;
double randNumber[randNbMax];

// Structure for positions
struct position {
    double x[DMax];
    double f;
    int size;
};

// Structure for search space
struct SS {
    int D;                    // Dimension of the search space
    double min[DMax];            // Lower bounds
    double max[DMax];            // Upper bounds
    struct q {
        double q[DMax];
        int size;
    } q;                    // Quantisation
    int quantisation;                // If >0, quantised search space
    double normalise;                // If >0, normalised search space
};

// Structure for fitness
struct fitness {
    double f[DMax];
    int size;
};

// Structure for problem definition
struct problem {
    int function;        // Type of function
    double epsilon;        // Precision
    double objective;        // Objective value
    int evalMax;        // Maximum number of evaluations
    struct SS SS;        // Search space
};

// Additional structures needed for perf.c
struct param {
    int S;              // Swarm size (number of particles)
    double K;           // Mean connectivity
    double p;           // Individual probability of being informed by the best previous position
    double w;           // First cognitive/confidence coefficient
    double c;           // Second cognitive/confidence coefficient

    int confin;         // Confinement method
    /*
     0 => hyperbolic (start over from the other side), also called "toroidal"
     1 => random with uniform distribution
     2 => random around the violated boundary
     3 => random with destination direction
     4 => bounce back
     5 => bounce back with random direction
     6 => stick to the boundary
     */

    int distrib;        // Initial distribution
    /*
     -1 => uniform in the ball (hypersphere)
     0 => uniform in the hypercube, along the radiuses
     (uniform along the variables, but not uniform in the cube)
     1=> Gaussian (Box-Muller implementation)
     2=> Gaussian (Central Method)
     3=> Stable distribution
     4=> Slash distribution
     */

    // Mean and standard deviation, if applicable
    double mean;
    double sigma;

    // Topology
    int topology;

    // Trace parameter. 0-> no verbose, 1-> verbose, 2-> medium verbose
    int trace;

    // BellsAndWhistle. An array of options
    int BW[10];
};

// Function prototypes
struct problem problemDef(int functionCode);
double perf(struct position x, int function, struct SS SS, double objective);
struct fitness constraint(struct position x, int functCode);
double lennard_jones(struct position x);
double alea(double a, double b, int randCase);
double alea_normal(double mean, double std_dev, int option);
double distanceL(struct position x, struct position y, int L);
double sign(double a);

#endif // VALIDATION_H 