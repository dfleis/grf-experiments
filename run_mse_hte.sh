#! /bin/bash

MODEL=hte

for setting in {1..5}; do
	Rscript run-mse-forest-small.R $MODEL $setting
done

for setting in {1..5}; do
	Rscript run-mse-forest.R $MODEL $setting
done
