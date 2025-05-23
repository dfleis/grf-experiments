source("utils/utils-bench-forest.R")
seed <- 1
FILENAME_HEAD <- "bench-forest-FORK"
NUM_THREADS <- 1

args <- commandArgs(TRUE)
model_type <- args[1] # "vcm" or "hte"
setting_id <- args[2] # "1", "2", "3", "4", "5"

# model_type <- "vcm"
# setting_id <- c(1, 2, 3, 4)
stumps <- c(TRUE, FALSE)
Kvals <- c(4, 16, 64, 256)
pvals <- 5
nvals <- as.integer(c(1e4, 2e4, 1e5))
numtreevals <- c(100)

stopifnot(model_type %in% c("vcm", "hte"))
stopifnot(setting_id %in% c("1", "2", "3", "4", "5"))

methods <- c("grad", "fpt1", "fpt2")
center_data <- TRUE
nrep <- 1
niter <- 3
.quiet <- FALSE

grf_args_global <- list(
  min.node.size = 5,
  sample.fraction = 0.5,
  honesty = TRUE,
  honesty.fraction = 0.5,
  honesty.prune.leaves = TRUE,
  alpha = 0.05,
  ci.group.size = 1,
  compute.oob.predictions = FALSE,
  num.threads = NUM_THREADS,
  seed = seed
)

bench_forest(
  methods = methods,
  model_type = model_type,
  setting_id = setting_id,
  numtreevals = numtreevals,
  Kvals = Kvals,
  pvals = pvals,
  nvals = nvals,
  center_data = center_data,
  grf_args = grf_args_global,
  nrep = nrep,
  niter = niter,
  seed = seed,
  filename_head = FILENAME_HEAD,
  .quiet = .quiet
)

