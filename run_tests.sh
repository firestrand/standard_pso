#!/bin/bash

# Get the absolute path of this script's directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Create base results directory
mkdir -p "$SCRIPT_DIR/results"

# Function to run a test
run_test() {
    local runs=$1
    local seed=$2
    
    
    echo "Running with $runs runs and seed $seed..."
    
    # Run the PSO program with appropriate arguments
    "$SCRIPT_DIR/build/standard_pso" -f 13 -R $runs -S 40 -sd $seed
}

# Run all test functions with 100 runs and different seeds
for seed in 1294404794 1234567890; do
    run_test 100 $seed
done

echo "All tests completed. Results are in the results directory." 