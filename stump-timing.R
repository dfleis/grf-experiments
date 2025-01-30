suppressMessages(library(tidyverse))
suppressMessages(library(bench))
library(grf)
source("utils-data.R")
source("utils-forest.R")

BENCH_QUIET <- TRUE
NUM_THREADS <- 1
FILENAME_HEAD <- "data/stump-timing"

args <- commandArgs(TRUE)
model_type <- args[1] # "vcm" or "hte"
setting_id <- args[2] # "1", "2", "3", "4", "5"
# model_type <- "hte"
# setting_id <- "1"

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

# Data params
Kvals <- c(4096, 1024, 256, 64, 16, 4)
pvals <- c(2, 5)
nvals <- c(10000, 20000)


prep_args_common <- function(X, method, stump, seed) {
  method <- match.arg(tolower(method), choices = c("grad", "fpt1", "fpt2"))
  
  # specify arguments for single tree tests
  num.trees = 1L
  sample.fraction = 1.0
  mtry = ncol(X)
  honesty = FALSE
  honesty.fraction = 1.0
  honesty.prune.leaves = FALSE
  ci.group.size = 1L
  compute.oob.predictions = FALSE
  num.threads = NUM_THREADS
  alpha = 0.0
  
  if (isTRUE(stump)) min.node.size = nrow(X) - 1
  else min.node.size = 5L
  
  # other (default) arguments
  sample.weights = NULL
  clusters = numeric(0)
  equalize.cluster.weights = FALSE
  imbalance.penalty = 0
  stabilize.splits = FALSE
  
  num.threads <- grf:::validate_num_threads(NULL)
  samples.per.cluster <- grf:::validate_equalize_cluster_weights(
    equalize.cluster.weights,
    clusters,
    sample.weights)
  
  args <- list(num.trees = num.trees,
               clusters = clusters,
               samples.per.cluster = samples.per.cluster,
               sample.fraction = sample.fraction,
               mtry = mtry,
               min.node.size = min.node.size,
               honesty = honesty,
               honesty.fraction = honesty.fraction,
               honesty.prune.leaves = honesty.prune.leaves,
               alpha = alpha,
               imbalance.penalty = imbalance.penalty,
               stabilize.splits = stabilize.splits,
               ci.group.size = ci.group.size,
               compute.oob.predictions = compute.oob.predictions,
               method.flag = switch(method, grad = 1, fpt1 = 2, fpt2 = 3),
               num.threads = num.threads,
               seed = seed,
               legacy.seed = grf:::get_legacy_seed(),
               sample.weights = sample.weights)
  return (args)
}
prep_data_hte <- function(data, method, stump, seed) {
  X <- data$X
  Y <- grf:::validate_observations(data$Y, data$X, allow.matrix = TRUE)
  Wmm <- stats::model.matrix(~data$W - 1)
  
  args <- prep_args_common(X, method, stump, seed)
  args$gradient.weights <- numeric(0)
  
  data <- grf:::create_train_matrices(
    X = X,
    outcome = Y,
    treatment = Wmm[,-1],
    sample.weights = args$sample.weights)
  
  args$sample.weights <- NULL
  
  return (list(data = data, args = args))
}
prep_data_vcm <- function(data, method, stump, seed) {
  X <- data$X
  Y <- grf:::validate_observations(data$Y, data$X, allow.matrix = TRUE)
  W <- data$W
  
  args <- prep_args_common(X, method, stump, seed)
  args$gradient.weights <- rep(1, NCOL(W) * NCOL(Y))
  
  data <- grf:::create_train_matrices(
    X = X,
    outcome = Y,
    treatment = W,
    sample.weights = args$sample.weights)
  
  args$sample.weights <- NULL
  
  return (list(data = data, args = args))
}

multi_causal_tree_wrapper <- function(prep_data) {
  # grf:::multi_causal_train underlies both grf::lm_forest and grf::multi_arm_causal_forest
  forest <- grf:::do.call.rcpp(grf:::multi_causal_train,
                               c(prep_data$data, prep_data$args))
  (length(forest$`_split_vars`[[1]]) - 1)/2
}

#--------------------------------------------------
#----- SIMULATIONS
#--------------------------------------------------
set.seed(seed)
par_grid <- expand.grid(K = Kvals, p = pvals, n = nvals)
FUN_grf <- get_FUN_grf(model_type)

t0 <- Sys.time()
sim_times <- vector(mode = "list", length = nrow(par_grid))
for (i in 1:nrow(par_grid)) {
  pars <- par_grid[i,]
  K <- pars$K
  p <- pars$p
  n <- pars$n
  
  cat(format(Sys.time()), "::", i, "of", nrow(par_grid), "::",
      paste(model_type, setting_id), "::", FUN_grf$name, "::",
      paste(names(pars), pars, sep = " = ", collapse = ", "), "\n")
  
  # Generate data
  data_train <- generate_data(n = n, p = p, K = K,
                              model_type = settings$model_type,
                              effect_type = settings$effect_type,
                              prob_type = settings$prob_type)
  
  bp <- bench::press(
    method = methods,
    {
      prep_data <- if (model_type == "vcm") {
        prep_data_vcm(data = data_train, method = method, stump = TRUE, seed = seed)
      } else if (model_type == "hte") {
        prep_data_hte(data = data_train, method = method, stump = TRUE, seed = seed)
      }
      
      bench::mark(
        multi_causal_tree_wrapper(prep_data),
        min_iterations = 50
      )
    },
    .quiet = BENCH_QUIET
  )
  
  sim_times[[i]] <- bp %>%
    select(method, median, n_itr, nsplits = result) %>%
    mutate(K = K, p = p, n = n,
           nt = nt,
           model_type = settings$model_type,
           setting_id = settings$setting_id,
           effect_type = settings$effect_type,
           prob_type = settings$prob_type,
           FUN = FUN_grf$name)
}
t1 <- Sys.time()


#--------------------------------------------------
#----- SAVE DATA
#--------------------------------------------------
df_times <- bind_rows(sim_times) %>%
  mutate(median = as.numeric(median)) %>%
  unnest_longer(nsplits) %>%
  unnest(nsplits)

df_times_wide <- df_times %>%
  select(-n_itr) %>%
  pivot_wider(
    names_from = method,
    values_from = median
  ) %>%
  mutate(fpt2_FACTOR = as.numeric(grad)/as.numeric(fpt2))

datetime <- format(Sys.time(), format = "%Y%m%d-%H%M")

make_filename <- function(label) {
  sprintf("%s-%s-%s-%s-%s.csv", FILENAME_HEAD, model_type, setting_id, label, datetime)
}

write.csv(df_times, file = make_filename("long"), row.names = FALSE)
write.csv(df_times_wide, file = make_filename("wide"), row.names = FALSE)
  
cat("Total elapsed time:", format(difftime(t1, t0), digits = 6), "\n")

