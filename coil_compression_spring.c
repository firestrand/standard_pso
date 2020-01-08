#include "main.h"

struct fitness constrain(struct position x) {
    static double Fmax = 1000.0;
    static double Fp = 300;
    double Cf;
    double K;
    double sp;
    double lf;

    static double S = 189000.0;
    static double lmax = 14.0;
    static double spm = 6.0;
    static double sw = 1.25;
    static double G = 11500000;
    struct fitness ff = {0};
    ff.size = 1; // Default value

    Cf = 1 + 0.75 * x.x[2] / (x.x[1] - x.x[2]) + 0.615 * x.x[2] / x.x[1];
    K = 0.125 * G * pow(x.x[2], 4) / (x.x[0] * x.x[1] * x.x[1] * x.x[1]);
    sp = Fp / K;
    lf = Fmax / K + 1.05 * (x.x[0] + 2) * x.x[2];

    ff.f[1] = 8 * Cf * Fmax * x.x[1] / (pi * x.x[2] * x.x[2] * x.x[2]) - S;
    ff.f[2] = lf - lmax;
    ff.f[3] = sp - spm;
    ff.f[4] = sw - (Fmax - Fp) / K;

    return ff;
}

double coil_compression_spring(struct position xs) {
    // Coil compression spring  (penalty method)
    // Ref New Optim. Tech. in Eng. p 644
    double c, f, x1, x2, x3;
    struct fitness ff;

    x1 = xs.x[0]; // {1,2, ... 70}
    x2 = xs.x[1];//[0.6, 3]
    x3 = xs.x[2];// relaxed form [0.207,0.5]  dx=0.001
    // In the original problem, it is a list of
    // acceptable values
    // {0.207,0.225,0.244,0.263,0.283,0.307,0.331,0.362,0.394,0.4375,0.5}

    f = pi * pi * x2 * x3 * x3 * (x1 + 2) * 0.25;
    //	f=x2*x3*x3*(x1+2);
    // Constraints
    ff = constrain(xs);

    if (ff.f[1] > 0) {
        c = 1 + ff.f[1];
        f = f * c * c * c;
    }
    if (ff.f[2] > 0) {
        c = 1 + ff.f[2];
        f = f * c * c * c;
    }
    if (ff.f[3] > 0) {
        c = 1 + ff.f[3];
        f = f * c * c * c;
    }
    if (ff.f[4] > 0) {
        c = 1 + pow(10, 10) * ff.f[4];
        f = f * c * c * c;
    }
    if (ff.f[5] > 0) {
        c = 1 + pow(10, 10) * ff.f[5];
        f = f * c * c * c;
    }
    return f;
}
