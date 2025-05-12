#! /bin/bash

MODEL=hte

for setting in {1..5}; do
	Rscript run-bench-tree-small.R $MODEL $setting
done

for setting in {1..5}; do
	Rscript run-bench-tree.R $MODEL $setting
done

for setting in {1..5}; do
	Rscript run-bench-forest-small.R $MODEL $setting
done

for setting in {1..5}; do
	Rscript run-bench-forest.R $MODEL $setting
done

