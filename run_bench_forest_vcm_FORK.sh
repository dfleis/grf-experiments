#! /bin/bash

MODEL=vcm
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

## Run all benchmarks with settings 1-4
#for setting in {1..4}; do
#    run_and_log "run-bench-forest-small-FORK.R" $setting
#done

for setting in {1..4}; do
    run_and_log "run-bench-forest-FORK.R" $setting
done
