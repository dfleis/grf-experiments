#!/bin/bash

# Print a message to indicate which script is running
echo "Running MSE scripts for HTE and VCM models..."
echo

# Make sure all scripts are executable
chmod +x run_mse_hte.sh
chmod +x run_mse_vcm.sh

# Run all scripts in sequence
echo "Running HTE MSE evaluations..."
./run_mse_hte.sh
echo

echo "Running VCM MSE evaluations..."
./run_mse_vcm.sh
echo

echo "All scripts completed."
