#!/bin/bash

# Build the validation program
make clean
make

# Run the validation program
./build/validate

echo "Validation data has been generated in build/data directory." 