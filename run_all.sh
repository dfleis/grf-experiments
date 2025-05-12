#!/bin/bash

# Print a message to indicate which script is running
echo "Running benchmark and MSE scripts for HTE and VCM models..."
echo

# Make sure all scripts are executable
chmod +x run_bench_tree_hte.sh
chmod +x run_bench_tree_vcm.sh
chmod +x run_bench_forest_hte.sh
chmod +x run_bench_forest_vcm.sh
chmod +x run_mse_hte.sh
chmod +x run_mse_vcm.sh

# Run all scripts in sequence
echo "Running tree VCM benchmarks..."
./run_bench_tree_vcm.sh
echo

echo "Running tree HTE benchmarks..."
./run_bench_tree_hte.sh
echo

echo "Running forest HTE benchmarks..."
./run_bench_forest_hte.sh
echo

echo "Running forest VCM benchmarks..."
./run_bench_forest_vcm.sh
echo

echo "Running HTE MSE evaluations..."
./run_mse_hte.sh
echo

echo "Running VCM MSE evaluations..."
./run_mse_vcm.sh
echo

echo "All scripts completed."
