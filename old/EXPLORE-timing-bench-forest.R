###################################################################################################################
# forest-grid-VCM.R
# Comparison methods "grad", "fpt1", and "fpt2" for GRF estimates of heterogeneous effects 
# under a varying coefficient/conditionally linear model.
#
# Use grf::lm_forest to target estimates of \theta^*(x) from the varying coefficient model
#   E[Y_i | X_i = x] = W_i^\top \theta^*(x).
# Samples (X_i, Y_i, W_i) such that
#   * X_i = p-dimensional auxiliary covariates,
#   * Y_i = scalar-valued response/outcome,
#   * W_i = K-dimensional primary regressors.
# We target estimates of the coefficient function \theta^*:\mathbb R^p\to\mathbb R^K
# such that \theta^*(x) has coordinates/coordinate functions
#   \theta^*(x) = (\theta^*_1(x), ..., \theta^*_K(x)),
# with component/coordinate functions \theta^*_k:\mathbb R\to\mathbb R.
#
# Benchmarks for grf::lm_forest over a grid of parameters:
#   * minimum node size (min.node.size) x number of trees (num.trees)
#   * K (regressor dimension) x p (covariate dimension) x n (training sample size)
#
# Generated 2025-Jan-10 23:16 EST
###################################################################################################################
# library(rlang)
library(tidyverse)
library(purrr)
library(bench)
library(grf)

source("utils-data.R")

#--------------------------------------------------
#----- HELPERS/WRAPPERS
#--------------------------------------------------
center_forests <- function(data, num_trees = 2000, min_node_size = 5, ...) {
  # Robinson-style centering step to produce centered Y - Ybar and W - Wbar.
  # This is the default mechanism that underlies grf::lm_forest. We want to
  # use the same centered Y and W matrices across all methods since FPT does
  # not affect this step, and we're only interested in comparing the differences 
  # between FPT and grad
  num_trees_centering <- max(50, round(num_trees/4))
  with(data, {
    forest_Y <- grf::multi_regression_forest(X = X, Y = Y, 
                                             num.trees = num_trees_centering,
                                             min.node.size = min_node_size, 
                                             compute.oob.predictions = TRUE, ...)
    forest_W <- grf::multi_regression_forest(X = X, Y = W, 
                                             num.trees = num_trees_centering, 
                                             min.node.size = min_node_size, 
                                             compute.oob.predictions = TRUE, ...)
    Y_hat <- predict(forest_Y)$predictions
    W_hat <- predict(forest_W)$predictions
    list(Y_hat = Y_hat, W_hat = W_hat)
  })
}

# process_forest <- function(forest, test) {
#   X_test <- test$X
#   thetaX_test <- test$thetaX
#   
#   preds <- predict(forest, newdata = X_test, estimate.variance = FALSE, drop = TRUE)
#   est <- preds$predictions
#   colnames(est) <- test$effect_labels
#   err <- est - thetaX_test
#   cmse <- colMeans(err^2) # component-wise MSE
#   amse <- mean(cmse)      # average MSE
#     
#   df_cmse <- data.frame(effect = names(cmse), cmse = cmse, row.names = NULL)
#     
#   # final nb. of splits per tree, after honest pruning (`_split_vars` is the total nodes)
#   nsplits <- (sapply(forest$`_split_vars`, length) - 1)/2 
#   df_nsplits <- data.frame(tree = 1:length(nsplits), nsplits = nsplits)
#     
#   list(amse = amse, cmse = df_cmse, nsplits = df_nsplits)
# }
process_forest <- function(forest, test) {
  X_test <- test$X
  thetaX_test <- test$thetaX
  
  preds <- predict(forest, newdata = X_test, estimate.variance = FALSE, drop = TRUE)
  est <- preds$predictions
  colnames(est) <- colnames(thetaX_test)
  err <- est - thetaX_test
  cmse <- colMeans(err^2) # component-wise MSE
  amse <- mean(cmse)      # average MSE
  
  df_cmse <- data.frame(effect = names(cmse), cmse = cmse, row.names = NULL)
  
  # final nb. of splits per tree, after honest pruning (`_split_vars` is the total nodes)
  nsplits <- (sapply(forest$`_split_vars`, length) - 1)/2 
  df_nsplits <- data.frame(tree = 1:length(nsplits), nsplits = nsplits)
  
  list(amse = amse, cmse = df_cmse, nsplits = df_nsplits)
}

process_bench_results <- function(bb, test) {
  df_stats <- bb$result %>%
    map(process_forest, test = test) %>% # same as lapply(bp$result, process_forest, test = test)
    setNames(nm = bb$method) %>%
    purrr::list_transpose()
  
  bb %>% 
    select(-expression, -result, -memory) %>%
    mutate(amse = unname(df_stats$amse), 
           cmse = df_stats$cmse, 
           nsplits = df_stats$nsplits)
}

#--------------------------------------------------
#----- DATA/FOREST SETTINGS
#--------------------------------------------------
seed <- 1
model_type <- "vcm"
sigma_eps <- 1 # response/outcome noise 

# Global arguments
args_grf_global <- list(
  compute.oob.predictions = FALSE, 
  num.threads = 1,
  #honesty.prune.leaves = FALSE,
  ci.group.size = 1, # DISABLE CI SAMPLING 
  seed = seed)
# We want to disable the CI/variance estimation procedure
# because 1) we're not interested in evaluating CI coverage
# in this simulation, and 2) the CI batch sampling mechanism
# operates on the same RNG as the observation sampler, and 
# so in order to keep the samples drawn identical across
# the "grad", "fpt1", and "fpt2" methods, we must disable
# CI tree batching.

# GRF params
methods <- c("grad", "fpt1", "fpt2")
#num_trees <- rev(seq(100, 1000, by = 100))
#min_node_size <- 5*2^(0:7)#seq(10, 20, by = 10)
num_trees <- 1
min_node_size <- 5

# Data params
Kvals <- rev(c(2, 4, 8, 16, 32, 64, 128, 256))
pvals <- 1#c(1, 5)
nvals <- 20000 #c(2500, 10000)


# Parameter grid
par_grid <- expand.grid(K = Kvals, p = pvals, n = nvals, nt = num_trees, mns = min_node_size)

#--------------------------------------------------
#----- SIMULATIONS
#--------------------------------------------------
set.seed(seed)

rlang::local_options(bench.press_quiet = TRUE)
Sys.sleep(1)

t0 <- Sys.time()
sim_times <- vector(mode = "list", length = nrow(par_grid))
for (i in 1:nrow(par_grid)) {
  pars <- par_grid[i,]
  K <- pars$K
  p <- pars$p
  n <- pars$n
  nt <- pars$nt
  mns <- pars$mns
  
  cat(format(Sys.time()), "::", i, "of", nrow(par_grid), "::",
      paste(names(pars), pars, sep = " = ", collapse = ", "), "\n")
  
  # Generate data
  data_train <- generate_data(n = n, p = p, K = K, model_type = model_type, sigma_eps = sigma_eps)
    
  # Common centering step as per the recommendations of grf::lm_forest (common for all three)
  centered <- center_forests(data = data_train, 
                             num_trees = nt, 
                             min_node_size = mns, 
                             seed = args_grf_global$seed,
                             num.threads = args_grf_global$num.threads)   
  args_grf <- modifyList(args_grf_global, 
                         list(X = data_train$X, Y = data_train$Y, W = data_train$W,
                              Y.hat = centered$Y_hat, W.hat = centered$W_hat,
                              num.trees = nt, min.node.size = mns))
  make_args <- function(mm) modifyList(args_grf, list(method = mm))
  
  # bt_times <- sapply(methods, function(method) {
  #   gc()
  #   bt <- bench::bench_time(do.call(grf::lm_forest, make_args(method)))
  #   return (bt)
  # })
  # 
  gc()
  
  bp <- bench::press(
    method = methods,
    {
      bench::mark(
        do.call(grf::lm_forest, make_args(method)),
        iterations = 1,
        filter_gc = TRUE # large forests appear to almost always trigger GC
      )
    }
  )
  
  sim_times[[i]] <- bp %>% 
    select(method, median) %>%
    mutate(K = K, p = p, n = n)
  
  rm(bp)
  
  #sim_times[[i]] <- data.frame(grad = bt_times["real","grad"],
  #           fpt1 = bt_times["real","fpt1"],
  #           fpt2 = bt_times["real","fpt2"]) %>%
  #  mutate(K = K, p = p, n = n)
}
t1 <- Sys.time()
df_times <- bind_rows(sim_times) %>% 
  mutate(trees = num_trees,
         node  = min_node_size)

df_wide <- df_times %>%
  pivot_wider(
    names_from = method, 
    values_from = median
  ) %>%
  mutate(fpt1_FACTOR = as.numeric(grad)/as.numeric(fpt1),
         fpt2_FACTOR = as.numeric(grad)/as.numeric(fpt2))

df_wide

cat("Total elapsed time:", format(difftime(t1, t0), digits = 4))

