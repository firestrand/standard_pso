#include "validation.h"

// Define constants for zero and infinity
#define zero 1.e-20
#define infinity 1.e20

// Define global variables declared as extern in validation.h
double errMin = 0;

// Function to generate a random position within the search space
struct position generate_random_position(struct SS SS) {
    struct position pos;
    pos.size = SS.D;
    
    for (int d = 0; d < SS.D; d++) {
        pos.x[d] = alea(SS.min[d], SS.max[d], 0);
    }
    
    return pos;
}

// Function to generate a position at the lower bounds
struct position generate_min_position(struct SS SS) {
    struct position pos;
    pos.size = SS.D;
    
    for (int d = 0; d < SS.D; d++) {
        pos.x[d] = SS.min[d];
    }
    
    return pos;
}

// Function to generate a position at the upper bounds
struct position generate_max_position(struct SS SS) {
    struct position pos;
    pos.size = SS.D;
    
    for (int d = 0; d < SS.D; d++) {
        pos.x[d] = SS.max[d];
    }
    
    return pos;
}

// Function to generate a position at the true optimum for each function
struct position generate_optimal_position(int functionCode, struct SS SS) {
    struct position pos;
    pos.size = SS.D;
    
    // Initialize with zeros as default
    for (int d = 0; d < SS.D; d++) {
        pos.x[d] = 0.0;
    }
    
    // Offset vectors from perf.c for CEC functions
    // Shifted Parabola/Sphere (CEC 2005 benchmark)
    static double offset_0[30] = {
        -3.9311900e+001, 5.8899900e+001, -4.6322400e+001, -7.4651500e+001, -1.6799700e+001,
        -8.0544100e+001, -1.0593500e+001, 2.4969400e+001, 8.9838400e+001, 9.1119000e+000,
        -1.0744300e+001, -2.7855800e+001, -1.2580600e+001, 7.5930000e+000, 7.4812700e+001,
        6.8495900e+001, -5.3429300e+001, 7.8854400e+001, -6.8595700e+001, 6.3743200e+001,
        3.1347000e+001, -3.7501600e+001, 3.3892900e+001, -8.8804500e+001, -7.8771900e+001,
        -6.6494400e+001, 4.4197200e+001, 1.8383600e+001, 2.6521200e+001, 8.4472300e+001
    };
    
    // Shifted Rosenbrock (CEC 2005 benchmark)
    static double offset_2[30] = {
        8.1023200e+001, -4.8395000e+001, 1.9231600e+001, -2.5231000e+000, 7.0433800e+001,
        4.7177400e+001, -7.8358000e+000, -8.6669300e+001, 5.7853200e+001, -9.9533000e+000,
        2.0777800e+001, 5.2548600e+001, 7.5926300e+001, 4.2877300e+001, -5.8272000e+001,
        -1.6972800e+001, 7.8384500e+001, 7.5042700e+001, -1.6151300e+001, 7.0856900e+001,
        -7.9579500e+001, -2.6483700e+001, 5.6369900e+001, -8.8224900e+001, -6.4999600e+001,
        -5.3502200e+001, -5.4230000e+001, 1.8682600e+001, -4.1006100e+001, -5.4213400e+001
    };
    
    // Shifted Rastrigin (CEC 2005)
    static double offset_3[30] = {
        1.9005000e+000, -1.5644000e+000, -9.7880000e-001, -2.2536000e+000, 2.4990000e+000,
        -3.2853000e+000, 9.7590000e-001, -3.6661000e+000, 9.8500000e-002, -3.2465000e+000,
        3.8060000e+000, -2.6834000e+000, -1.3701000e+000, 4.1821000e+000, 2.4856000e+000,
        -4.2237000e+000, 3.3653000e+000, 2.1532000e+000, -3.0929000e+000, 4.3105000e+000,
        -2.9861000e+000, 3.4936000e+000, -2.7289000e+000, -4.1266000e+000, -2.5900000e+000,
        1.3124000e+000, -1.7990000e+000, -1.1890000e+000, -1.0530000e-001, -3.1074000e+000
    };
    
    // Shifted Schwefel (F2 CEC 2005. Also for F4)
    static double offset_4[30] = {
        3.5626700e+001, -8.2912300e+001, -1.0642300e+001, -8.3581500e+001, 8.3155200e+001,
        4.7048000e+001, -8.9435900e+001, -2.7421900e+001, 7.6144800e+001, -3.9059500e+001,
        4.8885700e+001, -3.9828000e+000, -7.1924300e+001, 6.4194700e+001, -4.7733800e+001,
        -5.9896000e+000, -2.6282800e+001, -5.9181100e+001, 1.4602800e+001, -8.5478000e+001,
        -5.0490100e+001, 9.2400000e-001, 3.2397800e+001, 3.0238800e+001, -8.5094900e+001,
        6.0119700e+001, -3.6218300e+001, -8.5883000e+000, -5.1971000e+000, 8.1553100e+001
    };
    
    // Shifted Griewank (CEC 2005)
    static double offset_5[30] = {
        -2.7626840e+002, -1.1911000e+001, -5.7878840e+002, -2.8764860e+002, -8.4385800e+001,
        -2.2867530e+002, -4.5815160e+002, -2.0221450e+002, -1.0586420e+002, -9.6489800e+001,
        -3.9574680e+002, -5.7294980e+002, -2.7036410e+002, -5.6685430e+002, -1.5242040e+002,
        -5.8838190e+002, -2.8288920e+002, -4.8888650e+002, -3.4698170e+002, -4.5304470e+002,
        -5.0658570e+002, -4.7599870e+002, -3.6204920e+002, -2.3323670e+002, -4.9198640e+002,
        -5.4408980e+002, -7.3445600e+001, -5.2690110e+002, -5.0225610e+002, -5.3723530e+002
    };
    
    // Shifted Ackley (CEC 2005)
    static double offset_6[30] = {
        -1.6823000e+001, 1.4976900e+001, 6.1690000e+000, 9.5566000e+000, 1.9541700e+001,
        -1.7190000e+001, -1.8824800e+001, 8.5110000e-001, -1.5116200e+001, 1.0793400e+001,
        7.4091000e+000, 8.6171000e+000, -1.6564100e+001, -6.6800000e+000, 1.4543300e+001,
        7.0454000e+000, -1.8621500e+001, 1.4556100e+001, -1.1594200e+001, -1.9153100e+001,
        -4.7372000e+000, 9.2590000e-001, 1.3241200e+001, -5.2947000e+000, 1.8416000e+000,
        4.5618000e+000, -1.8890500e+001, 9.8008000e+000, -1.5426500e+001, 1.2722000e+000
    };
    
    // Set the known optimum for each function
    switch (functionCode) {
        case 0: // Parabola (Sphere)
            // Optimal at origin
            break;
            
        case 1: // Griewank
            // Optimal at origin
            break;
            
        case 2: // Rosenbrock
            // Optimal at (1,1,1,...)
            for (int d = 0; d < SS.D; d++) {
                pos.x[d] = 1.0;
            }
            break;
            
        case 3: // Rastrigin
            // Optimal at origin
            break;
            
        case 4: // Tripod
            pos.x[0] = 0.0;
            pos.x[1] = -50.0;
            break;
            
        case 5: // Ackley
            // Optimal at origin
            break;
            
        case 6: // Schwefel
            // Optimal at (420.9687,...,420.9687)
            for (int d = 0; d < SS.D; d++) {
                pos.x[d] = 420.9687;
            }
            break;
            
        case 7: // Schwefel 1.2
            // Optimal at origin
            break;
            
        case 8: // Schwefel 2.22
            // Optimal at origin
            break;
            
        case 9: // Neumaier 3
            // Optimal at (i+1) for each i
            for (int d = 0; d < SS.D; d++) {
                pos.x[d] = d + 1;
            }
            break;
            
        case 10: // G3 (constrained)
            // Optimal evenly distributed
            for (int d = 0; d < SS.D; d++) {
                pos.x[d] = 1.0 / sqrt((double)SS.D);
            }
            break;
            
        case 12: // Schwefel
            // Optimal at (420.9687,...,420.9687)
            for (int d = 0; d < SS.D; d++) {
                pos.x[d] = 420.9687;
            }
            break;
            
        case 13: // 2D Goldstein-Price
            pos.x[0] = 0.0;
            pos.x[1] = -1.0;
            break;
            
        case 14: // Schaffer f6
            pos.x[0] = 0.0;
            pos.x[1] = 0.0;
            break;
            
        case 15: // Step
            // Optimal at origin
            break;
            
        case 16: // Schwefel 2.21
            // Optimal at origin
            break;
            
        case 17: // Lennard-Jones
            // Complex 3D structure - optimum depends on number of atoms
            // Cannot be easily represented in a simple way
            break;
            
        case 18: // Gear train
            // Known optimal solution (16,19,43,49) or equivalent
            if (SS.D == 4) {
                pos.x[0] = 16.0;
                pos.x[1] = 19.0;
                pos.x[2] = 43.0;
                pos.x[3] = 49.0;
            }
            break;
            
        case 21: // Compression Spring
            // Known optimal solution from literature
            if (SS.D == 3) {
                pos.x[0] = 7.0;      // Number of active coils (integer)
                pos.x[1] = 1.3;      // Wire diameter
                pos.x[2] = 0.35;     // Mean coil diameter
            }
            break;
            
        case 25: // Pressure Vessel
            // Known optimal solution from literature
            if (SS.D == 4) {
                pos.x[0] = 0.8125;    // x1 = 0.8125 (thickness of shell)
                pos.x[1] = 0.4375;    // x2 = 0.4375 (thickness of head)
                pos.x[2] = 42.0;      // x3 = 42.0 (inner radius)
                pos.x[3] = 176.0;     // x4 = 176.0 (length of cylindrical section)
            }
            break;
            
        // CEC 2005 Benchmark Functions with known shift values
        case 100: // F1 CEC 2005 (shifted sphere)
            // Optimum is at the offset_0 values
            for (int d = 0; d < SS.D && d < 30; d++) {
                pos.x[d] = offset_0[d];
            }
            break;
            
        case 102: // F6 CEC 2005 (shifted Rosenbrock)
            // Optimum is at offset_2 for Rosenbrock
            // Note: In Rosenbrock, the true optimum is at offset + [1,1,1...]
            for (int d = 0; d < SS.D && d < 30; d++) {
                pos.x[d] = offset_2[d]; // The +1 is handled in perf.c
            }
            break;
            
        case 103: // F9 CEC 2005 (shifted Rastrigin)
            // Optimum is at offset_3 values
            for (int d = 0; d < SS.D && d < 30; d++) {
                pos.x[d] = offset_3[d];
            }
            break;
            
        case 104: // F2 CEC 2005 (shifted Schwefel)
        case 107: // F4 CEC 2005 (shifted Schwefel with noise)
            // Optimum is at offset_4 values
            for (int d = 0; d < SS.D && d < 30; d++) {
                pos.x[d] = offset_4[d];
            }
            break;
            
        case 105: // F7 CEC 2005 (shifted Griewank, not rotated)
            // Optimum is at offset_5 values
            for (int d = 0; d < SS.D && d < 30; d++) {
                pos.x[d] = offset_5[d];
            }
            break;
            
        case 106: // F8 CEC 2005 (shifted Ackley, not rotated)
            // Optimum is at offset_6 values
            for (int d = 0; d < SS.D && d < 30; d++) {
                pos.x[d] = offset_6[d];
            }
            break;
            
        case 999: // Test function
            // Depends on the specific test
            break;
            
        default:
            // For other functions, use origin as default
            break;
    }
    
    return pos;
}

// Function to save position and fitness value to CSV
void save_to_csv(FILE* file, struct position pos, double fitness, const char* label) {
    // Print label for the test point
    fprintf(file, "%s,", label);
    
    // Print position values separated by commas
    for (int d = 0; d < pos.size; d++) {
        fprintf(file, "%.15f,", pos.x[d]);
    }
    
    // Print fitness value and newline
    fprintf(file, "%.15f\n", fitness);
}

// Function to generate validation data for a specific function
void generate_validation_data(int functionCode) {
    struct problem pb = problemDef(functionCode);
    char filename[100];
    FILE* file;
    
    // Create filename with function number and dimension
    sprintf(filename, "build/data/F%d-D%d.csv", functionCode, pb.SS.D);
    
    // Open file for writing
    file = fopen(filename, "w");
    if (!file) {
        printf("Error: Could not open file %s for writing\n", filename);
        return;
    }
    
    // Write header information
    fprintf(file, "# Function: %d, Dimension: %d\n", functionCode, pb.SS.D);
    fprintf(file, "# Objective value: %.15f\n", pb.objective);
    
    // Add note about CEC functions
    if (functionCode >= 100 && functionCode <= 107) {
        fprintf(file, "# Note: This is a CEC 2005 benchmark function with shift/rotation. Optimum position represents the global optimum.\n");
    }
    
    fprintf(file, "# Search space: ");
    for (int d = 0; d < pb.SS.D; d++) {
        fprintf(file, "[%.6f,%.6f]", pb.SS.min[d], pb.SS.max[d]);
        if (d < pb.SS.D - 1) fprintf(file, " × ");
    }
    fprintf(file, "\n");
    
    // Write column headers
    fprintf(file, "point_type,");
    for (int d = 0; d < pb.SS.D; d++) {
        fprintf(file, "x%d,", d+1);
    }
    fprintf(file, "fitness\n");
    
    // Generate and evaluate min bound position
    struct position min_pos = generate_min_position(pb.SS);
    double min_fitness = perf(min_pos, functionCode, pb.SS, pb.objective);
    save_to_csv(file, min_pos, min_fitness, "min_bounds");
    
    // Generate and evaluate max bound position
    struct position max_pos = generate_max_position(pb.SS);
    double max_fitness = perf(max_pos, functionCode, pb.SS, pb.objective);
    save_to_csv(file, max_pos, max_fitness, "max_bounds");
    
    // Generate and evaluate global optimal position
    struct position opt_pos = generate_optimal_position(functionCode, pb.SS);
    double opt_fitness = perf(opt_pos, functionCode, pb.SS, pb.objective);
    save_to_csv(file, opt_pos, opt_fitness, "global_optimum");
    
    // Generate and evaluate 10 random positions
    for (int i = 0; i < 10; i++) {
        char label[20];
        sprintf(label, "random_%d", i+1);
        struct position rand_pos = generate_random_position(pb.SS);
        double rand_fitness = perf(rand_pos, functionCode, pb.SS, pb.objective);
        save_to_csv(file, rand_pos, rand_fitness, label);
    }
    
    fclose(file);
    printf("Generated validation data for function %d with dimension %d\n", functionCode, pb.SS.D);
}

int main() {
    // Seed the random number generator with a fixed seed for reproducibility
    srand(12345);
    
    // Create an array of function codes to test
    int function_codes[] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
        11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
        21, 22, 23, 24, 25, 26,
        100, 102, 103, 104, 105, 106, 107
    };
    int num_functions = sizeof(function_codes) / sizeof(function_codes[0]);
    
    printf("Generating validation data for %d functions...\n", num_functions);
    
    // Generate validation data for each function
    for (int i = 0; i < num_functions; i++) {
        generate_validation_data(function_codes[i]);
    }
    
    printf("Validation data generation complete!\n");
    return 0;
} 