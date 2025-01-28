#! /bin/bash

module load r/4.4.0

MODEL=vcm

for setting in {1..4}; do
	Rscript forest-timing.R $MODEL $setting
done
