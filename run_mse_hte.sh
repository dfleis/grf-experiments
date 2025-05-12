#! /bin/bash

MODEL=hte
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Create logs directory if it doesn't exist
mkdir -p logs/

# Define function to run script and log output
run_and_log() {
    local script=$1
    local setting=$2
    local script_base=$(basename "$script" .R | sed 's/run-//')
    
    local logfile="logs/${script_base}_${MODEL}_setting${setting}_${TIMESTAMP}.log"
    echo "Running $script_base model=$MODEL setting=$setting (logging to $logfile)"
    Rscript "$script" $MODEL $setting > "$logfile" 2>&1
}

# Run all MSE benchmarks with settings 1-5
for setting in {1..5}; do
    run_and_log "run-mse-forest-small.R" $setting
done

for setting in {1..5}; do
    run_and_log "run-mse-forest.R" $setting
done
