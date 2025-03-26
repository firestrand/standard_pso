#!/bin/bash

# Create results directory
mkdir -p results

# Function to run a test
run_test() {
    local func=$1
    local runs=$2
    local seed=$3
    local output_file="results/function_${func}_${runs}runs_seed${seed}.txt"
    echo "Running function $func with $runs runs and seed $seed..."
    ./build/standard_pso -R $runs -S 40 -f $func -sd $seed > $output_file
    echo "Results saved to $output_file"
}

# Run all test functions with 100 runs and different seeds
for seed in 1294404794 1234567890; do
    for func in 4 11 15 17 18 20 21 100 102 103 104 105 106; do
        run_test $func 100 $seed
    done
done

echo "All tests completed. Results are in the results directory." 