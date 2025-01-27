#!/bin/bash

# Slurm parameters
MEMO=249G		# memory allocation
TIME=0-08:00:00		# time allocation
CORE=55			# core allocation

CHDIR=/home/yyang/scratch/grf-experiments

# Assemble order prefix
ORDP="sbatch --account=def-yyang --mem="$MEMO" -n 1 -c "$CORE" --time="$TIME" --chdir="$CHDIR""

# Create log directory if it doesn't exist
LOGS=logs
mkdir -p $LOGS

# Assemble components for the slurm order for this job
JOBN="bench_NODE"
SCRIPT="bench-NODE.sh"

OUTF=$LOGS"/"$JOBN".out"
ERRF=$LOGS"/"$JOBN".err"

# Assemble slurm order for this job
ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT

# Print order
echo $ORDP

# Submit order
$ORD
