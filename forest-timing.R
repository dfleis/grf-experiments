suppressMessages(library(tidyverse))
suppressMessages(library(bench))
library(grf)
source("utils-data.R")
source("utils-forest.R")

BENCH_QUIET <- TRUE
NUM_THREADS <- 1
FILENAME_HEAD <- "data/forest-timing"

args <- commandArgs(TRUE)
model_type <- args[1] # "vcm" or "hte"
setting_id <- args[2] # "1", "2", "3", "4", "5"
##model_type <- "hte"
##setting_id <- "5"

stopifnot(model_type %in% c("vcm", "hte"))
stopifnot(setting_id %in% c("1", "2", "3", "4", "5"))

settings <- validate_setting(model_type = model_type, setting_id = setting_id)

#--------------------------------------------------
#----- DATA/FOREST SETTINGS
#--------------------------------------------------
seed <- 1

# Algorithm/model params
#methods <- c("grad", "fpt1", "fpt2")
methods <- c("grad", "fpt2")
num_trees <- c(1, 50, 100)

# Data params
Kvals <- c(256, 64, 16, 4)
pvals <- 5
nvals <- c(10000, 20000)

# Global GRF arguments
args_grf_global <- list(
  min.node.size = 5,
  compute.oob.predictions = FALSE, # just fit the forest
  ci.group.size = 1, # disable CI tree-wise batching
  num.threads = NUM_THREADS,
  seed = seed)

# Note on ci.group.size = 1.
# Each subsample is used to train ci.group.size trees. The conventional
# presentation of the GRF algorithm corresponds to ci.group.size = 1
# because every tree uses its own subsample. However, for CI estimation,
# GRF will compute trees in small batches were an individual subsample 
# is used to train ci.group.size trees.
# This mechanism affects the random number generator used in the
# C++ backend to do the subsampling in a way that makes it difficult
# to guarantee that two forests will use the same collection of 
# subsamples for training their trees. Specifically, in order for the
# drawn subsamples to be identical across the different methods 
# ("grad", "fpt1", "fpt2"), we disable CI batch sampling by setting
# ci.group.size = 1. This is fine for this simulation because we're 
# studying the forest fit times and not CI coverage.

#--------------------------------------------------
#----- SIMULATIONS
#--------------------------------------------------
set.seed(seed)
par_grid <- expand.grid(K = Kvals, p = pvals, n = nvals, nt = num_trees)
FUN_grf <- get_FUN_grf(model_type) 

t0 <- Sys.time()
sim_times <- vector(mode = "list", length = nrow(par_grid))
for (i in 1:nrow(par_grid)) {
  pars <- par_grid[i,]
  K <- pars$K
  p <- pars$p
  n <- pars$n
  nt <- pars$nt
  
  cat(format(Sys.time()), "::", i, "of", nrow(par_grid), "::",
      paste(model_type, setting_id), "::", FUN_grf$name, "::",
      paste(names(pars), pars, sep = " = ", collapse = ", "), "\n")
  
  # Generate data
  data_train <- generate_data(n = n, p = p, K = K, 
                              model_type = settings$model_type,
                              effect_type = settings$effect_type, 
                              prob_type = settings$prob_type)

  if (nt == 1) {
    disable_centering <- TRUE
  } else {
    disable_centering <- FALSE
  }
  
  # Common centering step as per the recommendations of grf::lm_forest and grf::multi_arm_causal_forest
  centered <- center_forests(data = data_train, 
                             model_type = model_type,
                             disable_centering = disable_centering,
                             num_trees = nt, 
                             min_node_size = args_grf_global$min.node.size, 
                             num.threads = args_grf_global$num.threads,
                             seed = args_grf_global$seed)   
  
  args_grf <- modifyList(args_grf_global, 
                         list(X = data_train$X, 
                              Y = data_train$Y, 
                              W = data_train$W,
                              Y.hat = centered$Y_hat, 
                              W.hat = centered$W_hat,
                              num.trees = nt))
  if (nt == 1) { 
    # disable subsampling and honesty if we're evaluating a single tree
    args_grf <- modifyList(args_grf, list(sample.fraction = 1, honesty = F))
  }
  make_args <- function(mm) modifyList(args_grf, list(method = mm))
  
  bp <- bench::press(
    method = methods,
    {
      bench::mark(
        do.call(FUN_grf$FUN, make_args(method)),
        iterations = 1,
        filter_gc = FALSE # large forests appear to almost always trigger GC
      )
    },
    .quiet = BENCH_QUIET
  )
  
  sim_times[[i]] <- bp %>% 
    select(method, median) %>%
    mutate(K = K, p = p, n = n, 
           nt = nt,
           model_type = settings$model_type,
           setting_id = settings$setting_id,
           effect_type = settings$effect_type,
           prob_type = settings$prob_type,
           FUN = FUN_grf$name)
  rm(bp)
}
t1 <- Sys.time()

#--------------------------------------------------
#----- SAVE DATA
#--------------------------------------------------
df_times <- bind_rows(sim_times) %>%
	mutate(median = as.numeric(median))
df_times_wide <- df_times %>%
  pivot_wider(
    names_from = method, 
    values_from = median
  ) %>%
  mutate(#fpt1_FACTOR = as.numeric(grad)/as.numeric(fpt1),
         fpt2_FACTOR = as.numeric(grad)/as.numeric(fpt2))

datetime <- format(Sys.time(), format = "%Y%m%d-%H%M")

make_filename <- function(label) {
  sprintf("%s-%s-%s-%s-%s.csv", FILENAME_HEAD, model_type, setting_id, label, datetime)
}

write.csv(df_times, file = make_filename("long"), row.names = FALSE)
write.csv(df_times_wide, file = make_filename("wide"), row.names = FALSE)

cat("Total elapsed time:", format(difftime(t1, t0), digits = 6), "\n")
