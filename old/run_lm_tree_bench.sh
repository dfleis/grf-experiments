#! /bin/bash

Kvals=(2 4 8 16)
pvals=(2)
nvals=(1000)

nrep=5
niter=5

Rscript grf-lm-forest-single-tree-bench.R $Kvals $pvals $nvals $nrep $niter

