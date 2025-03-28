/*
Standard PSO 2011 (from the Particle Swarm Central http://particleswarm.info)
+ some options. In particular
List Based Optimiser 
Mersenne RNG
Quasi-random numbers

 Contact for remarks, suggestions etc.:
Maurice.Clerc@WriteMe.com

	For more details, see ReadMe.txt		
*/

#include "main.h"
#include "wyhash.h"
#include <sys/stat.h>
#include <errno.h>
#include <unistd.h>

// Global variables defined in main.h
long double E;
double errMax = 0, errMin = 0;
double nbRand = 0;
int nBit = 0;
int nCycleMax = 0;
int nCycle = 0;
double pi;
double rMax = 0;
double randNumber[randNbMax];
int randRank = 0;
double randChaos = 0;

// For Network problem
int bcsNb = 0;
int btsNb = 0;

// For Repulsion problem
int Dim = 0;

// Files
FILE *f_run = NULL;
FILE *f_synth = NULL;
FILE *f_rand = NULL;
FILE *f_rand_bin = NULL;
FILE *f_trace = NULL;

// Function to create a directory if it doesn't exist
int create_directory(const char *path) {
    struct stat st = {0};
    if (stat(path, &st) == -1) {
#ifdef _WIN32
        return mkdir(path);
#else
        return mkdir(path, 0755);
#endif
    }
    return 0;
}

// Safe fopen wrapper
FILE* safe_fopen(const char *path, const char *mode) {
    FILE *file = fopen(path, mode);
    if (!file) {
        fprintf(stderr, "Error: Unable to open file %s: %s\n", path, strerror(errno));
        exit(EXIT_FAILURE);
    }
    return file;
}

// Function to get integer argument
int getIntArg(char *argv[], int *i, int argc, const char *flag) {
    if (*i + 1 >= argc) {
        fprintf(stderr, "Missing integer argument for %s\n", flag);
        exit(EXIT_FAILURE);
    }
    (*i)++;
    return atoi(argv[*i]);
}

// Function to get double argument
double getDoubleArg(char *argv[], int *i, int argc, const char *flag) {
    if (*i + 1 >= argc) {
        fprintf(stderr, "Missing double argument for %s\n", flag);
        exit(EXIT_FAILURE);
    }
    (*i)++;
    return atof(argv[*i]);
}

// =================================================
int main(int argc, char **argv) {
    struct position bestBest; // Best position over all runs
    int d;            // Current dimension
    double D;
    double error;            // Current error
    double errorMean;        // Average error
    double errorMin;        // Best result over all runs
    double errorMeanBest[R_max];
    double evalMean;        // Mean number of evaluations
    int functionCode;
    int func[funcMax]; // List of functions (codes) to optimise
    int indFunc;

    int nFailure;        // Number of unsuccessful runs
    double logProgressMean;
    struct param param;
    struct problem pb;
    int randCase;
    int run;
    struct result result;

    time_t seconds;
    char results_dir[256];
    char func_dir[256];
    char file_path[512];
    char base_dir[256];
    FILE *f_summary = NULL;

    int scanNb;
    double success[funcMax];
    double successRate;
    int t;
    double variance;
    float z;
    double zz;

    // Initialize default values
    int runMax = 1;
    double Smean = 40;
    int nbFunc = 1;
    int seed = 1294404794; // Default seed

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-R") == 0) {
            runMax = getIntArg(argv, &i, argc, "-R");
        } else if (strcmp(argv[i], "-S") == 0) {
            Smean = getDoubleArg(argv, &i, argc, "-S");
        } else if (strcmp(argv[i], "-f") == 0) {
            nbFunc = getIntArg(argv, &i, argc, "-f");
            if (nbFunc > funcMax) {
                fprintf(stderr, "Exceeded maximum number of functions.\n");
                exit(EXIT_FAILURE);
            }
        } else if (strcmp(argv[i], "-sd") == 0) {
            seed = getIntArg(argv, &i, argc, "-sd");
        } else {
            fprintf(stderr, "Unknown argument: %s\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }

    // Initialize global constants
    E = exp(1.0);
    pi = acos(-1.0);
    errMax = 0;
    nbRand = 0;

    // Get current working directory
    char cwd[256];
    if (getcwd(cwd, sizeof(cwd)) == NULL) {
        printf("\nWARNING: Could not get current working directory");
        strcpy(cwd, ".");
    }

    // Set base directory for results
    sprintf(base_dir, "%s/results", cwd);

    // Create main results directory
    printf("\nCreating results directory: '%s'", base_dir);
    if (create_directory(base_dir) != 0) {
        printf("\nWARNING: Could not create directory '%s': %s", base_dir, strerror(errno));
    }

    //------------------------------------------------ PARAMETERS
    // Bells and Whistles
    // Not really part of the standard
    // May improve the performance (not always)
    // Sometimes not well mathematically founded (rules of thumbs)
    // * => suggested value

    param.BW[0] = 0;    //	0 => same swarm size for each run
    //	1 => random swarm size around the given mean

    param.BW[1] = 0;    //*	0 => when P=G, use the "standard" method
    // 	1 => when P=G, a specific probabilistic method
    //	2 => when P=G, a more conservative method
    //  3 => when P=G, just look around
    // 4 =>  different weights for X, P, G  (TEST)

    param.BW[2] = 6;    // Randomness options
    // -2nn => Truncated KISS (simulated).
    //			nn is the number of bits you use to define
    //			each random number/ Example: -207 for 7 bits
    // -1 => "native" rand() of the C language
    //* 0 => pseudo-random number generator KISS
    //* 10 => pseudo-random number generator Mersenne 64 bits
    // 1 => quasi-random numbers for initialisation, Sobol sequences
    //      KISS after that
    // 2 => quasi-random numbers for initialisation, Halton sequences
    //      KISS after that
    // 3nn => Read on a list of bits (f_rand_bin.txt).
    //      Normally coming from a "true" (physical) random number generator
    //      (quantic system, atmospheric noise ...)
    //			nn is the number of bits you use to define
    //			each random number/ Example: 307 for 7 bits
    //      Warning: nn must be >=2
    // 4 => Read on a list (f_rand_quasi.txt)
    // 5 => quasi-random Rd and wyRand random numbers
    // 5 => wyRand random numbers
    f_rand_bin = fopen("f_rand_bin.txt", "r"); // A sequence of bits, if BW[2]=3
    f_rand = fopen("Ltest.txt", "r"); //A list of real numbers in ]0,1], if BW[2]=4

    param.BW[3] = 0;    // 1 => random numbering of the particles before each iteration
    // 0 => always the same loop "particle 0 to S-1"

    //--------

    param.confin = 0;    // 0 => keep inside the search space (supposed to be a D-rectangle)
    // 1 => no confinement
    //   WARNING: may be very slow (and bad) for discrete problems

    param.distrib = 0; // -1 => uniform in the hypersphere
    //* 0 => in the hypersphere, uniform along the radius
    // 			(and therefore NOT uniform in the sphere)
    // 1 => Gaussian (Box-Muller method). Warning: infinite loop possible
    // 2 =>	Gaussian (CMS method)
    // 3 => Other stable (CMS, experimental parameters)
    // 4 => Slash distribution (Gaussian BM/Gaussian BM)
    // Useful only if param.distrib>0;
    param.mean = 0.5; //Default: 0.5. For some functions 0 is better, though
    //	Example: shifted Rosenbrock (code 102)
    param.sigma = 1. / 12; // Default: 1./12 (standard deviation of U(0,1))
    // WARNING: the CMS method may not work with randomness option >=2
    if (param.BW[2] >= 2) param.distrib = 0;

    if (argc < 2) {
        Smean = 40; //Swarm size or Mean swarm size (if BW[0]=1). Suggested: 40
    }

    param.K = 3;    // Parameter to compute the probability p for a particle to be an
    // external informant. You may also directly define p (see below),
    // but K is about the mean number of the informants of a particle.
    // Default: 3

    // Confidence coefficients. Default:
    param.w = 1. / (2. * log((double) 2.)); // 0.721
    param.c = 0.5 + log((double) 2.); // 1.193
    param.topology = 0; // 0 => information links as in SPSO 2007 (quasi-random)
    // 1 => variable random ring (EXPERIMENTAL)

    //-------------------------------------------------- False randomnesses
    switch (param.BW[2]) //
    {
        default:
            break;
        case 4: // Prepare a list of false random number, read on a file
            t = 0;
        readRand:
            scanNb = fscanf(f_rand, "%f", &z);
            if (scanNb != EOF) {
                randNumber[t] = z;
                t = t + 1;
                goto readRand;
            }
            nCycleMax = t;
            printf("\n%i false random numbers read on a file", nCycleMax);

            break;
    }
    // ----------------------------------------------- PROBLEM
    param.trace = 1; // If >0 more information is displayed/saved (f_trace.txt)

    // Functions to optimise
    func[0] = 4; // 4
    func[1] = 11;  // 11
    func[2] = 15;  // 15
    func[3] = 17;  // 17
    func[4] = 18; // 18
    func[5] = 20;
    func[6] = 21;
    func[7] = 100;
    func[8] = 102;
    func[9] = 103;
    func[10] = 104;
    func[11] = 105;
    func[12] = 106;

    /* (see problemDef( ) for precise definitions)
        -1  Constant. For test of biases
    0 Parabola (Sphere)
    1 Griewank
    2 Rosenbrock (Banana)
    3 Rastrigin
    4 Tripod (dimension 2)
    5 Ackley
    6 Schwefel
    7 Schwefel 1.2
    8 Schwefel 2.2
    9 Neumaier 3
    10 G3
    11 Network optimisation (Warning: see problemDef() and also perf() for
            problem elements (number of BTS and BSC)
    12 Schwefel
    13 2D Goldstein-Price
    14 Schaffer f6
    15 Step
    16 Schwefel 2.21
    17 Lennard-Jones
    18 Gear train
    19 Sine_sine function
    20 Perm function
    21 Compression Spring
    22 Cellular phone (2D)
    23 Penalized
    24 Repulsion
    25 Pressure Vessel (penalty method)
    26 Ellipsoidal
    27 Quadric
    28 Frequency modulation sound parameter identification

    CEC 2005 benchmark  (no more than 30D. See cec2005data.c)
    100 F1 (shifted Parabola/Sphere)
    102 F6 (shifted Rosenbrock)
    103 F9 (shifted Rastrigin)
    104 F2 Schwefel
    105 F7 Griewank  (NOT rotated)
    106 F8 Ackley  (NOT rotated)
    107 F4 Schwefel + noise

    999 for tests
 

*/

    // Variables for overall summary
    double overallSuccess[funcMax] = {0};
    double overallAvgError[funcMax] = {0};
    double overallAvgIterations[funcMax] = {0};
    double overallAvgEvaluations[funcMax] = {0};

    if (runMax > R_max) {
        runMax = R_max;
        printf("\nWARNING. I can perform only %i runs. See R_max in main.h", R_max);
    }

    for (indFunc = 0; indFunc < nbFunc; indFunc++) // Loop on problems
    {
        functionCode = func[indFunc];

        // Create function-specific directory
        sprintf(func_dir, "%s/f%d", base_dir, functionCode);
        printf("\nCreating function directory: '%s'", func_dir);
        if (create_directory(func_dir) != 0) {
            printf("\nWARNING: Could not create directory %s: %s", func_dir, strerror(errno));
        }

        // Close any previously open files
        if (f_run) {
            fclose(f_run);
            f_run = NULL;
        }
        if (f_synth) {
            fclose(f_synth);
            f_synth = NULL;
        }
        if (f_trace) {
            fclose(f_trace);
            f_trace = NULL;
        }
        if (f_summary) {
            fclose(f_summary);
            f_summary = NULL;
        }

        // Create function-specific output files
        sprintf(file_path, "%s/f_run_seed%d.txt", func_dir, seed);
        printf("\nOpening run file: '%s'", file_path);
        f_run = safe_fopen(file_path, "w");
        if (!f_run) {
            printf("\nWARNING: Could not create file %s: %s", file_path, strerror(errno));
        }

        sprintf(file_path, "%s/f_synth_seed%d.txt", func_dir, seed);
        printf("\nOpening synth file: '%s'", file_path);
        f_synth = safe_fopen(file_path, "w");
        if (!f_synth) {
            printf("\nWARNING: Could not create file %s: %s", file_path, strerror(errno));
        }

        sprintf(file_path, "%s/f_trace_seed%d.txt", func_dir, seed);
        printf("\nOpening trace file: '%s'", file_path);
        f_trace = safe_fopen(file_path, "w");
        if (!f_trace) {
            printf("\nWARNING: Could not create file %s: %s", file_path, strerror(errno));
        }

        // Create summary file for all console output
        sprintf(file_path, "%s/summary_seed%d.txt", func_dir, seed);
        printf("\nOpening summary file: '%s'", file_path);
        f_summary = safe_fopen(file_path, "w");
        if (!f_summary) {
            printf("\nWARNING: Could not create file %s: %s", file_path, strerror(errno));
        }

        // Some information
        printf("\n Function %i ", functionCode);
        if (f_summary) fprintf(f_summary, "\n Function %i ", functionCode);

        // Define the problem
        pb = problemDef(functionCode);
        if (pb.SS.D > DMax) ERROR ("Can't solve it. You should increase DMax");

        // ----------------------------------------------- RUNS
        errorMean = 0;
        evalMean = 0;
        nFailure = 0;
        logProgressMean = 0; // Initialize logProgressMean
        D = pb.SS.D;
        randCase = param.BW[2];
        if (randCase > 300) {
            nBit = randCase - 300;
            randCase = 3;
        } // Decode the number of bits
        if (randCase < -200) {
            nBit = -randCase - 200;
            randCase = -2;
        }

        switch (randCase) {
            default:
                break;

            case 0:
                seed_rand_kiss(seed); // Initialise the RNG KISS for reproducible results
                break;

            case 10: // Mersenne 64 bits
                init_genrand64(seed);
                break;

            case -2: // Truncated KISS (simulated)
                rMax = pow(2, nBit) - 1;
                break;

            case 3: // The file is a string of bits
                //nBit=3;
                rMax = pow(2, nBit) - 1; // For conversion of a nBit string into a number in [0,1]
                break;

            case 4: // The file directly contains the numbers
                nCycle = 0;
                break;

            case 5:// wyRand w/ Rd quasi-random initialization
            case 6:// wyRand
                wysrand(seed);
                break;
        }


        randRank = 0;
        randChaos = 0.02;

        /*
            seconds=time(NULL); // Initialise the RNG KISS more randomly
             printf("\n time %ld",seconds);
            seed_rand_kiss(time(NULL));
    */
        switch (randCase) // "Warm up" the RNG for pseudo-random numbers
        {
            default:
                for (t = 0; t < 10000; t++) zz = alea(0, 1, randCase);
                break;

            case 3:
            case 4:
            case 5:
            case 6:
                break;
        }
//nCycle=4; 
        for (run = 0; run < runMax; run++) {
            if (param.BW[0] == 0) param.S = Smean; // Constant swarm size
            else // Random swarm size "around" the mean value
                param.S = (int) (0.5 *
                                 (0.5 + alea(Smean - D / 2, Smean + D / 2, 0) + alea(Smean - D / 2, Smean + D / 2, 0)));

            param.p = 1. - pow(1. - 1. / param.S, param.K);
//printf("\n p %f",param.p);      
            // (for a "global best" PSO, directly set param.p=1)

            printf("\n Swarm size %i", param.S);
            if (f_summary) fprintf(f_summary, "\n Swarm size %i", param.S);

            result = PSO(param, pb);
            error = result.error;

            if (error > pb.epsilon) // Failure
                nFailure = nFailure + 1;

            if (pb.SS.normalise > 0) {
                for (d = 0; d < pb.SS.D; d++)
                    result.SW.P[result.SW.best].x[d] =
                            pb.SS.min[d] + (pb.SS.max[d] - pb.SS.min[d]) * result.SW.P[result.SW.best].x[d]
                                           / pb.SS.normalise;
            }

            // Initialize bestBest on first run
            if (run == 0) {
                bestBest = result.SW.P[result.SW.best];
            }
                // Memorize the best (useful if more than one run)
            else if (error < bestBest.f) {
                bestBest = result.SW.P[result.SW.best];
            }

            // Result display
            errorMean = errorMean + error;
            printf("\nRun %i. S %i,  Eval %f. Error %e ", run + 1, param.S, result.nEval, error);
            printf(" Mean %e", errorMean / (run + 1));
            zz = 100 * (1 - (double) nFailure / (run + 1));
            printf("  Success  %.2f%%", zz);

            // Write the same info to the summary file
            if (f_summary) {
                fprintf(f_summary, "\nRun %i. S %i,  Eval %f. Error %e ", run + 1, param.S, result.nEval, error);
                fprintf(f_summary, " Mean %e", errorMean / (run + 1));
                fprintf(f_summary, "  Success  %.2f%%", zz);
            }

            // Best position display
            //for (d=0;d<pb.SS.D;d++) printf(" %f",result.SW.P[result.SW.best].x[d]);

            // Save result to f_run
            if (f_run) {
                fprintf(f_run, "\n%i %.1f %.0f %e %e ", run + 1, zz, result.nEval, error, errorMean / (run + 1));
                // Save best position
                for (d = 0; d < pb.SS.D; d++) fprintf(f_run, " %f", result.SW.P[result.SW.best].x[d]);
            } else {
                printf("\nWARNING: Could not write to run file (NULL pointer)");
            }

            // Compute/save some statistical information
            if (run == 0)
                errorMin = error;
            else if (error < errorMin)
                errorMin = error;

            evalMean = evalMean + result.nEval;
            errorMeanBest[run] = error;

            // Safely handle log calculation for error values
            if (error > 0) {
                logProgressMean = logProgressMean - log(error);
            }
        }        // End loop on "run"

        // ---------------------END

        // Display some statistical information
        evalMean = evalMean / (double) runMax;
        errorMean = errorMean / (double) runMax;
        logProgressMean = logProgressMean / (double) runMax;

        printf("\n Eval. (mean)= %f", evalMean);
        printf("\n Error (mean) = %e", errorMean);
        if (f_summary) {
            fprintf(f_summary, "\n Eval. (mean)= %f", evalMean);
            fprintf(f_summary, "\n Error (mean) = %e", errorMean);
        }

        // Variance
        variance = 0;

        for (run = 0; run < runMax; run++)
            variance = variance + pow(errorMeanBest[run] - errorMean, 2);

        variance = sqrt(variance / runMax);
        printf("\n Std. dev. %e", variance);
        printf("\n Log_progress (mean) = %f", logProgressMean);
        if (f_summary) {
            fprintf(f_summary, "\n Std. dev. %e", variance);
            fprintf(f_summary, "\n Log_progress (mean) = %f", logProgressMean);
        }

        // Success rate and minimum value
        printf("\n Failure(s) %i  ", nFailure);
        if (f_summary) fprintf(f_summary, "\n Failure(s) %i  ", nFailure);

        successRate = 100 * (1 - nFailure / (double) runMax);
        printf("Success rate = %.2f%%", successRate);
        if (f_summary) fprintf(f_summary, "Success rate = %.2f%%", successRate);

        success[indFunc] = successRate;

        printf("\n Best min value = %1.20e", errorMin);
        printf("\nPosition of the optimum: ");
        for (d = 0; d < pb.SS.D; d++) printf(" %.20f", bestBest.x[d]);

        if (f_summary) {
            fprintf(f_summary, "\n Best min value = %1.20e", errorMin);
            fprintf(f_summary, "\nPosition of the optimum: ");
            for (d = 0; d < pb.SS.D; d++) fprintf(f_summary, " %.20f", bestBest.x[d]);
        }

        // Save to f_synth file
        if (f_synth) {
            fprintf(f_synth, "%i %i %.0f %e %e %f %f", functionCode, pb.SS.D, successRate,
                    errorMean, variance, evalMean, bestBest.f);
            for (d = 0; d < pb.SS.D; d++) fprintf(f_synth, " %1.20e", bestBest.x[d]);
            fprintf(f_synth, "\n");

            // Specific save for Repulsive problem
            if (pb.function == 24) {
                for (d = 0; d < pb.SS.D - 1; d = d + 2) {
                    fprintf(f_synth, " %1.20e %1.20e", bestBest.x[d], bestBest.x[d + 1]);
                    fprintf(f_synth, "\n");
                }
            }
        } else {
            printf("\nWARNING: Could not write to synth file (NULL pointer)");
        }

        // Repeat informations in summary file
        if (f_summary) {
            fprintf(f_summary, "\n---------");
            fprintf(f_summary, "\n Function: %i", functionCode);
            fprintf(f_summary, "\n Confinement: %s", param.confin == 0 ? "YES" : "NO");
            fprintf(f_summary, "\n Distribution: ");
            switch (param.distrib) {
                case 0:
                    fprintf(f_summary, " uniform");
                    break;
                case 1:
                    fprintf(f_summary, " Gaussian (%f,%f), Box-Muller", param.mean, param.sigma);
                    break;
                case 2:
                    fprintf(f_summary, " Gaussian (%f,%f), CMS", param.mean, param.sigma);
                    break;
                case 3:
                    fprintf(f_summary, " Stable (%f,%f)", param.mean, param.sigma);
                    break;
                case 4:
                    fprintf(f_summary, " Slash (%f,%f)", param.mean, param.sigma);
                    break;
            }

            fprintf(f_summary, "\n BW = (%i, %i, %i, %i)", param.BW[0], param.BW[1],
                    param.BW[2], param.BW[3]);
            fprintf(f_summary, "\n Swarm size: ");
            if (param.BW[0] == 0)
                fprintf(f_summary, "%i", (int) Smean);
            else
                fprintf(f_summary, " mean %i", (int) Smean);
            fprintf(f_summary, "\n K = %i", param.K);
            fprintf(f_summary, "\n w = %f", param.w);
            fprintf(f_summary, "\n c = %f", param.c);
            fprintf(f_summary, "\n %e random numbers have been used", nbRand);
        } else {
            printf("\nWARNING: Could not write to summary file (NULL pointer)");
        }

        // Accumulate data for overall summary
        overallSuccess[indFunc] = successRate;
        overallAvgError[indFunc] = errorMean;
        overallAvgEvaluations[indFunc] = evalMean;

        // Close all files for this function
        printf("\nClosing output files for function %d", functionCode);
        if (f_run) {
            fclose(f_run);
            f_run = NULL;
        }
        if (f_synth) {
            fclose(f_synth);
            f_synth = NULL;
        }
        if (f_trace) {
            fclose(f_trace);
            f_trace = NULL;
        }
        if (f_summary) {
            fclose(f_summary);
            f_summary = NULL;
        }

    } // End "for ind[..."

    FILE *f_overall_summary = NULL;
    char overall_summary_path[512];

    sprintf(overall_summary_path, "%s/overall_summary_seed%d.txt", base_dir, seed);
    f_overall_summary = safe_fopen(overall_summary_path, "w");

    fprintf(f_overall_summary, "Overall Summary:\n");
    fprintf(f_overall_summary, "%-10s %-15s %-15s %-15s\n",
            "Function", "Success (%)", "Avg. Error", "Avg. Evaluations");

    for (indFunc = 0; indFunc < nbFunc; indFunc++) {
        fprintf(f_overall_summary, "f%-9d %-15.2f %-15e %-15.2f\n",
                func[indFunc],
                overallSuccess[indFunc],
                overallAvgError[indFunc],
                overallAvgEvaluations[indFunc]);
    }

    fclose(f_overall_summary);

    printf("\n errMax : %f", errMax);

    // Repeat informations
    printf("\n---------");
    printf("\n Function(s):");
    for (indFunc = 0; indFunc < nbFunc; indFunc++) // Loop on problems
    {
        functionCode = func[indFunc];
        printf(" %i", functionCode);
    }
    printf("\n Confinement: ");
    if (param.confin == 0) printf("YES"); else printf("NO");
    printf("\n Distribution: ");
    switch (param.distrib) {
        case 0:
            printf(" uniform");
            break;
        case 1:
            printf(" Gaussian (%f,%f), Box-Muller", param.mean, param.sigma);
            break;
        case 2:
            printf(" Gaussian (%f,%f), CMS", param.mean, param.sigma);
            break;
        case 3:
            printf(" Stable (%f,%f)", param.mean, param.sigma);
            break;
        case 4:
            printf(" Slash (%f,%f)", param.mean, param.sigma);
            break;
    }

    printf("\n BW = (%i, %i, %i, %i)", param.BW[0], param.BW[1],
           param.BW[2], param.BW[3]);
    printf("\n Swarm size: ");
    if (param.BW[0] == 0) printf("%i", (int) Smean); else printf(" mean %i", (int) Smean);
    printf("\n K = %i", param.K);
    printf("\n w = %f", param.w);
    printf("\n c = %f", param.c);
    printf("\n %e random numbers have been used", nbRand);

    // Close common files
    if (f_rand_bin) fclose(f_rand_bin);
    if (f_rand) fclose(f_rand);

    return 0; // End of main program
}
// ===============================================================
#include "alea.c"
#include "perf.c"
#include "problemDef.c"
#include "PSO.c"
#include "tools.c"


