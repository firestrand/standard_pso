#!/bin/bash

# Function to run a test with given parameters
run_test() {
    local func_num=$1
    local num_runs=$2
    local seed=$3
    
    echo "Running PSO for function $func_num..."
    
    # Create results directory if it doesn't exist
    mkdir -p results/f${func_num}
    
    # Run the PSO algorithm
    ./build/standard_pso -f ${func_num} -R ${num_runs} -S 40 -sd ${seed}
}

# Build the project
./build.sh

# Check if function number is provided
if [ $# -eq 1 ]; then
    # Run specific function
    run_test $1 100 1294404794
else
    # Run all CEC 2005 benchmark functions
    #run_test 103 100 1294404794  # F9 Shifted Rastrigin
    #run_test 104 100 1294404794  # F2 Shifted Schwefel
    #run_test 105 100 1294404794  # F7 Shifted Griewank
    run_test 106 100 1294404794  # F8 Shifted Ackley
fi

echo "All tests completed. Results are in the results directory." 