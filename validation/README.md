# Standard PSO Validation

This directory contains code to generate validation data for the Standard PSO benchmark functions.

## Overview

The validation system generates CSV files containing test points and fitness values for each benchmark function. These files can be used to validate implementations of the benchmark functions in other languages or frameworks.

## Building and Running

To build and run the validation system:

```bash
./build.sh
```

This script will:
1. Clean any previous build artifacts
2. Compile the validation program
3. Run the validation program to generate the CSV files
4. Output the validation data to the `build/data` directory

## Output Format

The validation system generates one CSV file for each benchmark function with the naming format `F<function_number>-D<dimension>.csv`. 

Each CSV file contains:
- Header information (function number, dimension, objective value, search space bounds)
- Column headers for position coordinates and fitness
- Test points:
  - Minimum bounds point
  - Maximum bounds point
  - Global optimum position
  - 10 random points with their fitness values

## CEC 2005 Benchmark Functions

For the CEC 2005 benchmark functions (function codes 100-107), the global optimum positions are set according to known shift vectors for specific dimensions:
- For dimension 30 (F1 F9): Representative shift vectors are provided
- For dimension 10 (F2, F4, F6, F7, F8): Representative shift vectors are provided
- For other dimensions: Pattern-based values are used to approximate the shifted optima

The fitness values for all points, including the global optima, are computed using the actual function evaluation to ensure they are accurate.

## Using the Validation Data

The validation data can be used to verify that your implementation of the benchmark functions produces the correct fitness values for the given test points.

For example, to validate an implementation of function 0 (Sphere function):

1. Read the test points from `build/data/F0-D30.csv`
2. For each test point, compute the fitness value using your implementation
3. Compare your computed fitness with the value in the CSV file
4. If the values match (within a small tolerance), your implementation is correct

## Available Functions

The validation system generates data for all 34 benchmark functions defined in `problemDef.c`, including:
- Classical test functions (Sphere, Griewank, Rosenbrock, etc.)
- Engineering test functions
- CEC 2005 benchmark functions (F1, F2, F4, F6, F7, F8, F9)

## License

See the project's main license file.