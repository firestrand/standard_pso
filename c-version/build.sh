#!/bin/bash

# Create build directory if it doesn't exist
mkdir -p build

# Create results directory
mkdir -p results

# Navigate to build directory
cd build

# Configure with CMake
cmake ..

# Build the project
cmake --build .

# Check if build was successful
if [ $? -eq 0 ]; then
    echo "Build successful!"
    # Return to project root directory
    cd ..
    echo "Running PSO with first function (f4 Tripod)..."
    ./build/standard_pso -f 1 -R 1
else
    echo "Build failed!"
    exit 1
fi 