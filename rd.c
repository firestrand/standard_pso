#include "main.h"
//
// Created by Silvers, Travis on 1/7/20.
//
/**
 * Newton-Raphson-Method to calculate g = phi_d.
 * @param D
 * @return
 */
double gammaRd(int D) {
    double x = 1.0;
    for (int i = 0; i < 20; i++) {
        x -= (pow(x, D + 1) - x - 1.0) / ((D + 1) * pow(x, D) - 1);
    }
    return x;
}

/**
 * Generates the first 'nRand' points in the 'R_d' sequence in D dimensional space.
 */
struct vectorList quasiRd(int D, int nRand, double seed) {
    int i, j;
    double out;
    double g = gammaRd(D);
    double alpha[D];
    struct vectorList qRand;

    for(i = 0; i < D; i++){
        alpha[i] = modf(pow(1.0 / g, i + 1), &out);
    }
    for(i = 0; i < D; i++){
        qRand.V[0].v[i] = modf(seed + alpha[i], &out);
    }
    for(j = 1; j < nRand; j++){
        for(i = 0; i < D; i++){
            qRand.V[j].v[i] = modf(qRand.V[j-1].v[i] + alpha[i], &out);
        }
    }
    return qRand;
}

