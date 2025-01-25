rm(list = ls())

suppressMessages(library(tidyverse))
library(grf)
source("utils-bench-tree.R")

args <- commandArgs(TRUE)
# model <- args[1]
# stumps <- args[2]

#model <- "lm"
#stumps <- "trees"

stopifnot(model %in% c("lm", "multi_arm"))
stopifnot(stumps %in% c("stumps", "trees"))
  
methods <- c("grad", "fpt1", "fpt2") 
Kvals <- c(2, 4, 8, 16, 32, 64, 128)
pvals <- c(1, 5, 25)
nvals <- c(2500, 5000, 10000)

nrep <- 10 # reps using new samples under the same parameter regime
niter <- 50 # iterations using the same samples within each rep
seed_r <- 1
seed_cpp <- 1

path <- "data/bench-tree"
filename_head <- sprintf("bench-%s-%s", model, stumps)

bench_tree(model = model,
           stumps = stumps, 
           methods = methods, 
           Kvals = Kvals, 
           pvals = pvals, 
           nvals = nvals, 
           nrep = nrep, 
           niter = niter, 
           seed_r = seed_r, 
           seed_cpp = seed_cpp,
           path = path,
           filename_head = filename_head,
           bench_quiet = FALSE)


