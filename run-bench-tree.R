source("utils/utils-bench-tree.R")
seed <- 1
FILENAME_HEAD <- "bench-tree"

args <- commandArgs(TRUE)
model_type <- args[1] # "vcm" or "hte"
setting_id <- args[2] # "1", "2", "3", "4", "5"

# model_type <- "vcm"
# setting_ids <- c(1, 2, 3, 4)
stumps <- c(TRUE, FALSE)
Kvals <- c(4, 16, 64, 256)
pvals <- c(2, 5)
nvals <- as.integer(c(1e4, 5e4, 1e5, 5e5))

stopifnot(model_type %in% c("vcm", "hte"))
stopifnot(setting_id %in% c("1", "2", "3", "4", "5"))

methods <- c("grad", "fpt1", "fpt2")
nrep <- 3
niter <- 5
.quiet <- FALSE

bench_tree(
  methods = methods,
  model_type = model_type,
  setting_ids = setting_ids,
  stumps = stumps,
  Kvals = Kvals,
  pvals = pvals,
  nvals = nvals,
  nrep = nrep,
  niter = niter,
  seed = seed, 
  filename_head = FILENAME_HEAD,
  .quiet = .quiet
)
