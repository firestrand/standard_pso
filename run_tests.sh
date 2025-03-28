#!/bin/bash

# Get the absolute path of this script's directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Create base results directory
mkdir -p "$SCRIPT_DIR/results"

# Function to run a test
run_test() {
    local func=$1
    local runs=$2
    local seed=$3
    local func_dir="$SCRIPT_DIR/results/f${func}"
    
    # Create function directory if it doesn't exist
    mkdir -p "$func_dir"
    
    echo "Running function $func with $runs runs and seed $seed..."
    
    # Run the PSO program with appropriate arguments
    "$SCRIPT_DIR/build/standard_pso" -R $runs -S 40 -f $func -sd $seed
    
    echo "Results saved to $func_dir"
}

# Run all test functions with 100 runs and different seeds
for seed in 1294404794 1234567890; do
    for func in 4 11 15 17 18 20 21 100 102 103 104 105 106; do
        run_test $func 100 $seed
    done
done

echo "All tests completed. Results are in the results directory." 