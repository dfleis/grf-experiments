#! /bin/bash

module load r/4.4.0

MODEL=hte

for setting in {1..5}; do
	Rscript forest-mse.R $MODEL $setting
done
