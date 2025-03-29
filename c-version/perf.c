double perf(struct position x, int function, struct SS SS,
            double objective) {                // Evaluate the fitness value for the particle of rank s
    double beta;
    double c;
    int d;
    double DD;
    double dx1, dx2;
    int grid;
    int i, j;
    int k;
    double min, max;
    int n;
    struct fitness ff = {0};
    double f, p, xd, x1, x2, x3, x4, x5, x6;
    double s11, s12, s21, s22;
    double sum1, sum2;
    double t0, tt, t1;
    double theta;
    double u;
    struct position xs;
    double y1, y2;
    double z;

    // Shifted Parabola/Sphere (CEC 2005 benchmark)
    static double offset_0[30] =
            {
                    -3.9311900e+001, 5.8899900e+001, -4.6322400e+001, -7.4651500e+001, -1.6799700e+001,
                    -8.0544100e+001, -1.0593500e+001, 2.4969400e+001, 8.9838400e+001, 9.1119000e+000,
                    -1.0744300e+001, -2.7855800e+001, -1.2580600e+001, 7.5930000e+000, 7.4812700e+001,
                    6.8495900e+001, -5.3429300e+001, 7.8854400e+001, -6.8595700e+001, 6.3743200e+001,
                    3.1347000e+001, -3.7501600e+001, 3.3892900e+001, -8.8804500e+001, -7.8771900e+001,
                    -6.6494400e+001, 4.4197200e+001, 1.8383600e+001, 2.6521200e+001, 8.4472300e+001
            };
    // Shifted Rosenbrock (CEC 2005 benchmark)
    static double offset_2[30] =
            {
                    8.1023200e+001, -4.8395000e+001, 1.9231600e+001, -2.5231000e+000, 7.0433800e+001,
                    4.7177400e+001, -7.8358000e+000, -8.6669300e+001, 5.7853200e+001, -9.9533000e+000,
                    2.0777800e+001, 5.2548600e+001, 7.5926300e+001, 4.2877300e+001, -5.8272000e+001,
                    -1.6972800e+001, 7.8384500e+001, 7.5042700e+001, -1.6151300e+001, 7.0856900e+001,
                    -7.9579500e+001, -2.6483700e+001, 5.6369900e+001, -8.8224900e+001, -6.4999600e+001,
                    -5.3502200e+001, -5.4230000e+001, 1.8682600e+001, -4.1006100e+001, -5.4213400e+001
            };
    // Shifted Rastrigin (CEC 2005)
    static double offset_3[30] =
            {
                    1.9005000e+000, -1.5644000e+000, -9.7880000e-001, -2.2536000e+000, 2.4990000e+000,
                    -3.2853000e+000, 9.7590000e-001, -3.6661000e+000, 9.8500000e-002, -3.2465000e+000,
                    3.8060000e+000, -2.6834000e+000, -1.3701000e+000, 4.1821000e+000, 2.4856000e+000,
                    -4.2237000e+000, 3.3653000e+000, 2.1532000e+000, -3.0929000e+000, 4.3105000e+000,
                    -2.9861000e+000, 3.4936000e+000, -2.7289000e+000, -4.1266000e+000, -2.5900000e+000,
                    1.3124000e+000, -1.7990000e+000, -1.1890000e+000, -1.0530000e-001, -3.1074000e+000
            };

    // Shifted Schwefel (F2 CEC 2005. Also for F4)
    static double offset_4[30] =
            {
                    3.5626700e+001, -8.2912300e+001, -1.0642300e+001, -8.3581500e+001, 8.3155200e+001,
                    4.7048000e+001, -8.9435900e+001, -2.7421900e+001, 7.6144800e+001, -3.9059500e+001,
                    4.8885700e+001, -3.9828000e+000, -7.1924300e+001, 6.4194700e+001, -4.7733800e+001,
                    -5.9896000e+000, -2.6282800e+001, -5.9181100e+001, 1.4602800e+001, -8.5478000e+001,
                    -5.0490100e+001, 9.2400000e-001, 3.2397800e+001, 3.0238800e+001, -8.5094900e+001,
                    6.0119700e+001, -3.6218300e+001, -8.5883000e+000, -5.1971000e+000, 8.1553100e+001
            };

    // Shifted Griewank (CEC 2005)
    static double offset_5[30] =
            {
                    -2.7626840e+002, -1.1911000e+001, -5.7878840e+002, -2.8764860e+002, -8.4385800e+001,
                    -2.2867530e+002, -4.5815160e+002, -2.0221450e+002, -1.0586420e+002, -9.6489800e+001,
                    -3.9574680e+002, -5.7294980e+002, -2.7036410e+002, -5.6685430e+002, -1.5242040e+002,
                    -5.8838190e+002, -2.8288920e+002, -4.8888650e+002, -3.4698170e+002, -4.5304470e+002,
                    -5.0658570e+002, -4.7599870e+002, -3.6204920e+002, -2.3323670e+002, -4.9198640e+002,
                    -5.4408980e+002, -7.3445600e+001, -5.2690110e+002, -5.0225610e+002, -5.3723530e+002
            };

    // Shifted Ackley (CEC 2005)
    static double offset_6[30] =
            {
                    -1.6823000e+001, 1.4976900e+001, 6.1690000e+000, 9.5566000e+000, 1.9541700e+001,
                    -1.7190000e+001, -1.8824800e+001, 8.5110000e-001, -1.5116200e+001, 1.0793400e+001,
                    7.4091000e+000, 8.6171000e+000, -1.6564100e+001, -6.6800000e+000, 1.4543300e+001,
                    7.0454000e+000, -1.8621500e+001, 1.4556100e+001, -1.1594200e+001, -1.9153100e+001,
                    -4.7372000e+000, 9.2590000e-001, 1.3241200e+001, -5.2947000e+000, 1.8416000e+000,
                    4.5618000e+000, -1.8890500e+001, 9.8008000e+000, -1.5426500e+001, 1.2722000e+000
            };



    //--------- Test
    struct problem pb;
    struct param param;

    // See function 11 (Network)
/*
	 static float bts[5][2]=
	 {
		 {6.8, 9.0},
		 {8.3, 7.9},
		 {6.6, 5.6},
		 {10, 5.4},
		 {8, 3} 
	 };
*/
    static float bts[19][2] =
            {
                    {6,  9},
                    {8,  7},
                    {6,  5},
                    {10, 5},
                    {8,  3},
                    {12, 2},
                    {4,  7},
                    {7,  3},
                    {1,  6},
                    {8,  2},
                    {13, 12},
                    {15, 7},
                    {15, 11},
                    {16, 6},
                    {16, 8},
                    {18, 9},
                    {3,  7},
                    {18, 2},
                    {20, 17}
            };

    float btsPenalty = 100;

    double z1, z2;

    if (SS.normalise > 0) {
        // Back to the real search space
        xs.size = x.size;
        for (d = 0; d < xs.size; d++)
            xs.x[d] = SS.min[d] + (SS.max[d] - SS.min[d]) * x.x[d] / SS.normalise;
    } else xs = x;

    switch (function) {
        case 100: // Parabola (Sphere) CEC 2005
            f = -450;
            for (d = 0; d < xs.size; d++) {
                xd = xs.x[d] - offset_0[d];
                f = f + xd * xd;
            }
            break;

        case 102:  // Rosenbrock
            // Solution point on offset_2 => fitness value = 390
            for (d = 0; d < xs.size; d++) {
                xs.x[d] = xs.x[d] - offset_2[d] + 1;
            }

            f = 390;

            for (d = 1; d < xs.size; d++) {
                tt = xs.x[d - 1] - 1;
                f = f + tt * tt;

                tt = xs.x[d - 1] * xs.x[d - 1] - xs.x[d];
                f = f + 100 * tt * tt;
            }

            //f=log(1+f);
            break;

        case 103: // Rastrigin
            for (d = 0; d < xs.size; d++) {
                xs.x[d] = xs.x[d] - offset_3[d];
            }
            f = -330;
            k = 10;

            for (d = 0; d < xs.size; d++) {
                xd = xs.x[d];
                f = f + xd * xd - k * cos(2 * pi * xd);
            }
            f = f + xs.size * k;
            break;

        case 104: // Schwefel (F2)
        case 107: // Schwefel + noise (F4)
            for (d = 0; d < xs.size; d++) {
                xs.x[d] = xs.x[d] - offset_4[d];
            }

            f = 0;
            for (d = 0; d < xs.size; d++) {
                sum2 = 0.0;
                for (k = 0; k <= d; k++) {
                    sum2 += xs.x[k];
                }
                f += sum2 * sum2;
            }
            if (function == 107) f = f * (1 + 0.4 * fabs(alea_normal(0, 1, 10)));
            f = f - 450;
            break;

        case 105: // Griewank. WARNING: in the CEC 2005 benchmark it is rotated
            sum1 = 0.0;
            sum2 = 1.0;
            f = -180;
            for (d = 0; d < xs.size; d++) {
                xd = xs.x[d] - offset_5[d];
                sum1 += xd * xd;
                sum2 *= cos(xd / sqrt(1.0 + d));
            }
            f = f + 1.0 + sum1 / 4000.0 - sum2;
            break;

        case 106: // Ackley
            f = -140;

            sum1 = 0.0;
            sum2 = 0.0;
            for (d = 0; d < xs.size; d++) {
                xd = xs.x[d] - offset_6[d];
                sum1 += xd * xd;
                sum2 += cos(2.0 * pi * xd);
            }
            sum1 = -0.2 * sqrt(sum1 / xs.size);
            sum2 /= xs.size;
            f = f + 20.0 + E - 20.0 * exp(sum1) - exp(sum2);
            break;

        case -1:        // Constant. For test of biases
            f = 0;
            break;

        case 0:        // Parabola (Sphere)
            f = 0;

            for (d = 0; d < xs.size; d++) {
                xd = xs.x[d] - d;
                f = f + xd * xd;
            }
            break;

        case 1:        // Griewank
            f = 0;
            p = 1;

            for (d = 0; d < xs.size; d++) {
                xd = xs.x[d];
                f = f + xd * xd;
                p = p * cos(xd / sqrt((double) (d + 1)));
            }
            f = f / 4000 - p + 1;
            break;

        case 2:        // Rosenbrock
            f = 0;
            t0 = xs.x[0] + 1;    // Solution on (0,...0) when
            // offset=0
            for (d = 1; d < xs.size; d++) {

                t1 = xs.x[d] + 1;
                tt = 1 - t0;
                f += tt * tt;
                tt = t1 - t0 * t0;
                f += 100 * tt * tt;
                t0 = t1;
            }
            break;

        case 3:        // Rastrigin
            k = 10;
            f = 0;

            for (d = 0; d < xs.size; d++) {
                xd = xs.x[d];
                f = f + xd * xd - k * cos(2 * pi * xd);
            }
            f = f + xs.size * k;
            break;

        case 4:        // 2D Tripod function
            // Note that there is a big discontinuity right on the solution
            // point.
            x1 = xs.x[0];
            x2 = xs.x[1];
            s11 = (1.0 - sign(x1)) / 2;
            s12 = (1.0 + sign(x1)) / 2;
            s21 = (1.0 - sign(x2)) / 2;
            s22 = (1.0 + sign(x2)) / 2;

            //f = s21 * (fabs (x1) - x2); // Solution on (0,0)
            f = s21 * (fabs(x1) + fabs(x2 + 50)); // Solution on (0,-50)
            f = f + s22 * (s11 * (1 + fabs(x1 + 50) +
                                  fabs(x2 - 50)) + s12 * (2 +
                                                          fabs(x1 - 50) +
                                                          fabs(x2 - 50)));

            //f=log(1+f);
            break;

        case 5:  // Ackley
            sum1 = 0;
            sum2 = 0;
            DD = x.size;
            pi = acos(-1);
            for (d = 0; d < x.size; d++) {
                xd = xs.x[d];
                sum1 = sum1 + xd * xd;
                sum2 = sum2 + cos(2 * pi * xd);
            }
            f = -20 * exp(-0.2 * sqrt(sum1 / DD)) - exp(sum2 / DD) + 20 + exp(1);

            break;

        case 6: // Schwefel
            f = 0;
            for (d = 0; d < x.size; d++) {
                xd = xs.x[d];
                f = f - xd * sin(sqrt(fabs(xd)));
            }
            break;

        case 7: // Schwefel 1.2
            f = 0;
            for (d = 0; d < x.size; d++) {
                xd = xs.x[d];
                sum1 = 0;
                for (k = 0; k <= d; k++) sum1 = sum1 + xd;
                f = f + sum1 * sum1;
            }
            break;

        case 8: // Schwefel 2.22
            sum1 = 0;
            sum2 = 1;
            for (d = 0; d < x.size; d++) {
                xd = fabs(xs.x[d]);
                sum1 = sum1 + xd;
                sum2 = sum2 * xd;
            }
            f = sum1 + sum2;
            break;

        case 9: // Neumaier 3
            sum1 = 0;
            sum2 = 1;
            for (d = 0; d < x.size; d++) {
                xd = xs.x[d] - 1;
                sum1 = sum1 + xd * xd;
            }
            for (d = 1; d < x.size; d++) {
                sum2 = sum2 + xs.x[d] * xs.x[d - 1];
            }

            f = sum1 + sum2;
            break;

        case 10: // G3 (constrained)
            // min =0 on (1/sqrt(D), ...)
            f = 1;
            sum1 = 0;
            for (d = 0; d < x.size; d++) {
                xd = xs.x[d];
                f = f * xd;
                sum1 = sum1 + xd * xd;
            }
            f = fabs(1 - pow(x.size, x.size / 2) * f) + x.size * fabs(sum1 - 1);
            break;

        case 11: // Network  btsNb BTS, bcdNb BSC

            f = 0;
            // Constraint: each BTS has one link to one BSC
            for (d = 0; d < btsNb; d++) {
                sum1 = 0;
                for (k = 0; k < bcsNb; k++) sum1 = sum1 + xs.x[d + k * btsNb];
                if (sum1 < 1 - zero || sum1 > 1 + zero) f = f + btsPenalty;

            }
            // Distances
            for (d = 0; d < bcsNb; d++) //For each BCS d
            {
                for (k = 0; k < btsNb; k++) // For each BTS k
                {
                    if (xs.x[k + d * btsNb] < 1) continue;
                    // There is a link between BTS k and BCS d
                    n = bcsNb * btsNb + 2 * d;
                    z1 = bts[k][0] - xs.x[n];
                    z2 = bts[k][1] - xs.x[n + 1];
                    f = f + sqrt(z1 * z1 + z2 * z2);
                }
            }
            break;

        case 12: // Schwefel
            f = 0;
            for (d = 0; d < x.size; d++) {
                xd = xs.x[d];
                f = f - xd * sin(sqrt(fabs(xd)));
            }
            break;

        case 13: // 2D Goldstein-Price function
            x1 = xs.x[0];
            x2 = xs.x[1];

            f = (1 + pow(x1 + x2 + 1, 2) * (19 - 14 * x1 + 3 * x1 * x1 - 14 * x2 + 6 * x1 * x2 + 3 * x2 * x2))
                * (30 + pow(2 * x1 - 3 * x2, 2) *
                        (18 - 32 * x1 + 12 * x1 * x1 + 48 * x2 - 36 * x1 * x2 + 27 * x2 * x2));
            break;

        case 14:  //Schaffer F6
            x1 = xs.x[0];
            x2 = xs.x[1];
            f = 0.5 + (pow(sin(sqrt(x1 * x1 + x2 * x2)), 2) - 0.5) / pow(1.0 + 0.001 * (x1 * x1 + x2 * x2), 2);

            break;

        case 15: // Step
            f = 0;
            for (d = 0; d < x.size; d++) {
                xd = (int) (xs.x[d] + 0.5);
                f = f + xd * xd;
            }
            break;

        case 16: // Schwefel 2.21
            f = 0;
            for (d = 0; d < x.size; d++) {
                xd = fabs(xs.x[d]);
                if (xd > f) f = xd;
            }
            break;

        case 17: // Lennard-Jones
            f = lennard_jones(xs);
            break;

        case 18: // Gear train
            f = pow(1. / 6.931 - x.x[0] * x.x[1] / (x.x[2] * x.x[3]), 2);
            //	f=pow(fabs(1./6.0 -x.x[0]*x.x[1]/(x.x[2]*x.x[3])),2);

            break;

        case 19: // Sine-sine function
            f = 0;
            for (d = 0; d < x.size; d++) {
                xd = xs.x[d];
                f = f - sin(xd) * pow(sin((d + 1) * xd * xd / pi), 20);

            }
            break;

        case 20: // Perm function
            beta = 10;
            f = 0;
            for (k = 0; k < x.size; k++) {
                sum1 = 0;
                for (d = 0; d < x.size; d++) {
                    xd = xs.x[d];
                    sum1 = sum1 + (pow(d + 1, k) + beta) * (pow(xd / (d + 1), k) - 1);
                }
                sum1 = sum1 * sum1;
                f = f + sum1;
            }

            break;

        case 21: // Coil compression spring  (penalty method)
            // Ref New Optim. Tech. in Eng. p 644

            x1 = xs.x[0]; // {1,2, ... 70}
            x2 = xs.x[1];//[0.6, 3]
            x3 = xs.x[2];// relaxed form [0.207,0.5]  dx=0.001
            // In the original problem, it is a list of
            // acceptable values
            // {0.207,0.225,0.244,0.263,0.283,0.307,0.331,0.362,0.394,0.4375,0.5}

            f = pi * pi * x2 * x3 * x3 * (x1 + 2) * 0.25;
            //	f=x2*x3*x3*(x1+2);
            // Constraints
            ff = constraint(xs, function);

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
            break;

        case 22: //Cellular phone
            // Grid 100 x 100. You may modify it for more or less precision
            // (which implies, of course, more or less computation time)
            grid = 100;
            min = 0;
            max = 1; // Warning. Hard coded. Must be the same as in problemDef
            //f=0;

            dx1 = (max - min) / grid;
            dx2 = (max - min) / grid;

            // For each point of the grid, compute the maximum field generated
            // by the stations. The aim is to maximize the smallest maximum
            f = infinity;

            for (i = 0; i <= grid; i++) {
                y1 = min + i * dx1;

                for (j = 0; j <= grid; j++) {
                    y2 = min + j * dx2;
                    z = 0; // Max known field
                    for (d = 0; d < xs.size - 1; d = d + 2) // Loop on station positions
                    {
                        x1 = xs.x[d]; // First coordinate
                        x2 = xs.x[d + 1]; // Second coordinate

                        //	z2=1./((x1-i)*(x1-y1) +(x2-j)*(x2-y2)+1);
                        // Field generated by the station (d, d+1)
                        z2 = 1. / ((x1 - y1) * (x1 - y1) + (x2 - y2) * (x2 - y2) + 0.0001 * dx1);
                        // If higher than already known, keep it
                        if (z2 > z) z = z2;
                    }

                    //	f=f+z;
                    // At this point, the maximum field generated is z
                    // If it is smaller than in previous checked points, keep it
                    if (z < f) f = z;
                }
            }
            // We want it as high as possible
            f = 1. / f; // In order to have something to minimise
            break;

        case 23: // Penalized
            f = pow(sin(pi * xs.x[0]), 2);
            for (d = 1; d < x.size - 1; d++) {
                f = f + pow(xs.x[d], 2) * (1 + pow(sin(3 * pi * xs.x[d + 1]), 2));
            }
            f = 0.1 * (f + pow(xs.x[x.size - 2], 2) * (1 + pow(sin(2 * pi * xs.x[x.size - 1]), 2)));

            for (d = 0; d < x.size; d++) {
                xd = xs.x[d];
                if (xd > 5) {
                    u = 100 * pow(xd - 5, 4);
                    f = f + u;
                }
                if (xd < -5) {
                    u = 100 * pow(-xd - 5, 4);
                    f = f + u;
                }
            }

            break;

        case 24: // Repulsion
            min = 0;
            max = 1; // WARNING. Hard coded. Must be the same as in problemDef

            n = (int) xs.size / Dim; // Number of points

            u = 0.5 * (n - 1); //1./100; // For the repulsive action of the boundaries
            // The smaller u, the nearer to the boundaries can be some points

            f = 0; // Total repulsive force to minimise
            for (i = 0; i < n; i++) // For each charged point ...
            {
                // ... add the repulsive force between the point and the boundaries
                for (d = 0; d < Dim; d++) {
                    f = f + u * (1. / pow(xs.x[d + Dim * i] - min, 2) + 1. / pow(max - xs.x[d + Dim * i], 2));
                }

                // ... add the repulsive force between the charged points
                for (j = 0; j < n; j++) {
                    if (j == i) continue;

                    // Distance^2 between point i and point j
                    z = 0;
                    for (d = 0; d < Dim; d++)
                        z = z + pow(xs.x[d + Dim * i] - xs.x[d + Dim * j], 2);

                    // Repulsion
                    f = f + 1. / z;
                }
            }
            break;

        case 25: // Pressure vessel (penalty method)
        case 1025: // confinement method

            // Ref New Optim. Tech. in Eng. p 638
            // D = 4

            x1 = xs.x[0]; // [1.1,12.5] granularity 0.0625
            x2 = xs.x[1];// [0.6,12.5] granularity 0.0625
            x3 = xs.x[2]; // [0,240]
            x4 = xs.x[3];// [0,240]

            f = 0.6224 * x1 * x3 * x4 + 1.7781 * x2 * x3 * x3 + 3.1611 * x1 * x1 * x4 + 19.84 * x1 * x1 * x3;
            //		f=0.6224*x1*x3*x4 + 1.7781*x2*x3*x3 + 3.1161*x1*x1*x4 + 19.84*x1*x1*x3;
            //  	f=0.6224*x1*x3*x4 + 1.7781*x2*x3*x3 + 3.1661*x1*x1*x4 + 19.84*x1*x1*x3;

            ff = constraint(xs, function);

            if (ff.f[1] > 0) {
                c = 1 + pow(10, 10) * ff.f[1];
                f = f * c * c;
            }
            if (ff.f[2] > 0) {
                c = 1 + ff.f[2];
                f = f * c * c;
            }
            if (ff.f[3] > 0) {
                c = 1 + ff.f[3];
                f = f * c * c;
            }

            break;

        case 26: // Ellipsoidal
            f = 0;

            for (d = 0; d < xs.size; d++) {
                xd = xs.x[d] - d - 1;
                f = f + xd * xd;
            }
            break;

        case 27:        // Quadric
            f = 0;
            for (d = 0; d < xs.size; d++) {
                xd = xs.x[0];
                if (d > 0)
                    for (j = 1; j < d; j++) xd = xd + xs.x[j];

                f = f + xd * xd;
            }
            break;

        case 28:// Frequency modulation sound parameter identification
            theta = 2 * pi / 100;
            f = 0;
            x1 = xs.x[0];
            x2 = xs.x[1];
            x3 = xs.x[2];
            x4 = xs.x[3];
            x5 = xs.x[4];
            x6 = xs.x[5];

            for (d = 1; d <= 100; d++) {
                z = x1 * sin(x2 * d * theta + x3 * sin(x4 * d * theta + x5 * sin(x6 * d * theta)))
                    - sin(5 * d * theta + 1.5 * sin(4.8 * d * theta + 2 * sin(4.9 * d * theta)));
                f = f + z * z;
            }
            break;

        case 999: // for tests
            f = 0;
            for (d = 0; d < xs.size; d++) {
                xd = xs.x[d];
                f = f + pow(xd, 4) - 16 * xd * xd + 5 * xd;
            }
            f = f / xs.size;
            break;

            // Rana
            x1 = xs.x[0];
            x2 = xs.x[1];
            f = x1 * sin(sqrt(fabs(x2 + 1 - x1))) * cos(sqrt(fabs(x2 + 1 + x1)))
                + (x2 + 1) * cos(sqrt(fabs(x2 + 1 - x1))) * sin(sqrt(fabs(x2 + 1 + x1)));
            break;

            // Goldstein & Price function
            x1 = xs.x[0];
            x2 = xs.x[1];
            f = (1 + pow(x1 + x2 + 1, 2) * (19 - 14 * x1 + 13 * x1 * x1 - 14 * x2 + 6 * x1 * x2 + 3 * x2 * x2));
            break;

            // Six-Hump Camel Back
            x1 = xs.x[0];
            x2 = xs.x[1];
            f = 4 * x1 * x1 - 2.1 * pow(x1, 4) + pow(x1, 6) / 3 + x1 * x2 - 4 * x2 * x2 + 4 * pow(x2, 4);

            break;
            // Periodic
            x1 = xs.x[0];
            x2 = xs.x[1];
            f = 1 + pow(sin(x1), 2) + pow(sin(x2), 2) - 0.1 * exp(-x1 * x1 - x2 * x2);

            break;

        case 1004: // For Tripod
            pb.function = 4;
            pb.SS.D = 2;    // Dimension

            // Boundaries
            for (d = 0; d < pb.SS.D; d++) {
                pb.SS.min[d] = -100;
                pb.SS.max[d] = 100;
                pb.SS.q.q[d] = 0;
            }

            pb.evalMax = 100;     // 10000
            pb.epsilon = 0.0001;
            pb.objective = 0;
            //--
            param.BW[0] = 0;
            param.BW[1] = 0;
            param.BW[2] = 4;
            param.BW[3] = 0;
            param.confin = 0;
            param.distrib = 0;
            param.mean = 0.5;
            param.sigma = 1. / 12; // Default: 1./12 (standard deviation of U(0,1))

            param.S = 40;
            param.K = 3;
            param.p = 1 - pow(1 - 1. / param.S, param.K);

            param.w = 1. / (2 * log((double) 2)); // 0.721
            param.c = 0.5 + log((double) 2); // 1.193
            param.topology = 0;
            param.trace = 0;

            nCycle = 0;
            nCycleMax = 3;

            for (d = 0; d < 3; d++) {
                randNumber[d] = xs.x[d];
                printf("%f ", randNumber[d]);
            }
            f = 1;
            for (i = 0; i < 2; i++) // 100
            {
                if (PSO(param, pb).error < pb.epsilon) f = f + 1;
            }
            f = 1. / f; // Minimise the inverse of the success rate

            break;

    }
//------------------
    f = fabs(f - objective);
    if (f < errMin) errMin = f; // For information
    if (f > errMax) { if (f < infinity) errMax = f; else errMax = infinity; } // For information

    return f;
}

//==========================================================
struct fitness constraint(struct position x, int functCode) {
    // ff[0] is defined in perf()
    // Variables specific to Coil compressing spring
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

    switch (functCode) {
        case 21: // Compression Spring
            Cf = 1 + 0.75 * x.x[2] / (x.x[1] - x.x[2]) + 0.615 * x.x[2] / x.x[1];
            K = 0.125 * G * pow(x.x[2], 4) / (x.x[0] * x.x[1] * x.x[1] * x.x[1]);
            sp = Fp / K;
            lf = Fmax / K + 1.05 * (x.x[0] + 2) * x.x[2];

            ff.f[1] = 8 * Cf * Fmax * x.x[1] / (pi * x.x[2] * x.x[2] * x.x[2]) - S;
            ff.f[2] = lf - lmax;
            ff.f[3] = sp - spm;
            ff.f[4] = sw - (Fmax - Fp) / K;
            break;

        case 25: // Pressure vessel
            ff.f[1] = 0.0193 * x.x[2] - x.x[0];
            ff.f[2] = 0.00954 * x.x[2] - x.x[1];
            ff.f[3] = 750 * 1728 - pi * x.x[2] * x.x[2] * (x.x[3] + (4.0 / 3) * x.x[2]);
            break;

    }
    return ff;
}
