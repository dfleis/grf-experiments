#!/bin/bash

# Slurm parameters
MEMO=10G		# memory allocation
TIME=0-01:00:00		# time allocation
CORE=10			# core allocation

CHDIR=/home/yyang/scratch/grf-experiments

# Assemble order prefix
ORDP="sbatch --account=def-yyang --mem="$MEMO" -n 1 -c "$CORE" --time="$TIME" --chdir="$CHDIR""



echo $ORDP
