#!/bin/bash

# Slurm parameters
MEMO=30G       # memory allocation
TIME=0-05:00:00 # time allocation
CORE=1          # core allocation

CHDIR=/home/yyang/scratch/grf-experiments

# Assemble order prefix
ORDP="sbatch --account=def-yyang --mem="$MEMO" -n 1 -c "$CORE" --time="$TIME" --chdir="$CHDIR""

# Create log directory if it doesn't exist
LOGS=logs
mkdir -p $LOGS

# Loop over both model types
for MODEL in vcm hte; do
	# Assemble components for the slurm order for this job
	JOBN="forest_timing_${MODEL}"
	SCRIPT="run_forest_timing_${MODEL}.sh"

	OUTF=$LOGS"/"$JOBN".out"
	ERRF=$LOGS"/"$JOBN".err"

	# Assemble slurm order for this job
	ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT

	# Print order
	echo $ORD

	# Submit order
	$ORD
done
