#! /bin/bash

MODEL=vcm

for setting in {1..4}; do
	Rscript run-mse-forest-small.R $MODEL $setting
done

for setting in {1..4}; do
	Rscript run-mse-forest.R $MODEL $setting
done
