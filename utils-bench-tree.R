suppressMessages(library(tidyverse))
library(bench)

generate_data_lm_tree <- function(K, p, n) {
  eps <- rnorm(n, 0, 1)
  X   <- matrix(rnorm(n * p), nrow = n)
  W   <- matrix(rnorm(n * K), nrow = n)
  
  beta <- rnorm(K)
  hX <- sapply(beta, function(b) b * X[,1])
  
  Y <- rowSums(hX * W) + eps 
  list(X = X, Y = Y, W = W)
}

generate_data_multi_arm_tree <- function(K, p, n) {
  stopifnot("Must have K >= 2 treatment levels (including the baseline treatment)." = K >= 2)
  stopifnot("It must be possible to observe all treatment levels for this test (n >= K)." = K <= n)
  # Note: Given n samples of a K-level treatment assignment, if K is so large 
  # that we don't observe every treatment level as part of the n samples, then GRF
  # will drop any unobserved levels. For the sake of this comparison, we don't
  # want this to happen since it's effectively fitting a different sized model.
  # In fact, for our degenerate test of fitting a single split, GRF will 
  # return the root node rather than performing any split.
  
  eps <- rnorm(n, 0, 1)
  X   <- matrix(rnorm(n * p), nrow = n)

  prob_W <- rep(1, K)/K # treatment level probabilities
  W_levels <- factor(1:K)
  # W <- sample(W_levels, size = n, replace = T, prob = prob_W)
  # force observation of each level at least once
  W <- c(W_levels, sample(W_levels, size = n - K, replace = T, prob = prob_W))
  W <- sample(W)
  Wmm <- stats::model.matrix(~-1 + W) # dummy indicators

  beta <- c(0, rnorm(K - 1))
  hX <- sapply(beta, function(b) b * X[,1])
  
  Y <- rowSums(hX * Wmm) + eps 
  list(X = X, Y = Y, W = W)
}

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
  num.threads = 1L
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

prep_data_multi_arm_tree <- function(data, method, stump, seed) {
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

prep_data_lm_tree <- function(data, method, stump, seed) {
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

prep_data_tree <- function(model = c("multi_arm", "lm"), pars, seed) {
  model <- match.arg(model)
  
  if (model == "multi_arm") {
    data <- generate_data_multi_arm_tree(pars$K, pars$p, pars$n)
    prep_data <- prep_data_multi_arm_tree(data, pars$method, pars$stump, seed)
  }
  else if (model == "lm") {
    data <- generate_data_lm_tree(pars$K, pars$p, pars$n)
    prep_data <- prep_data_lm_tree(data, pars$method, pars$stump, seed)
  }
  else stop("Something went wrong in prep_data_tree... Model should either be 'multi_arm' or 'lm'.")
  
  return (prep_data)
}

multi_causal_tree_wrapper <- function(prep_data) {
  # grf:::multi_causal_train underlies both grf::lm_forest and grf::multi_arm_causal_forest
  forest <- grf:::do.call.rcpp(grf:::multi_causal_train, 
                               c(prep_data$data, prep_data$args))
  (length(forest$`_split_vars`[[1]]) - 1)/2
}

parse_bench_output <- function(bench_data) {
  bench_data %>%
    mutate(splits = result) %>%
    unnest_longer(c(gc, time, splits), indices_include = T) %>%
    unnest_wider(gc) %>%
    rename(gc0 = level0, gc1 = level1, gc2 = level2, iter = time_id) %>%
    mutate(gc   = gc0 + gc1 + gc2 > 0,
           time = as.numeric(time)) %>%
    select(method, rep, iter, gc, gc0, gc1, gc2, splits, time)
}

bench_par_grid <- function(model, stumps, methods, pars, nrep, niter, seed_r, seed_cpp) {
  stopifnot(stumps %in% c("stumps", "trees"))
  stump <- stumps == "stumps" # fit stump or full trees
  
  # benchmark over parameter grid
  set.seed(seed_r)
  res <- bench::press(
    rep = 1:nrep,
    method = methods,
    {
      pars$method <- method
      pars$stump <- stump
      prep_data <- prep_data_tree(model, pars, seed_cpp) 
      
      bench::mark(
        multi_causal_tree_wrapper(prep_data),
        iterations = niter,
        min_time = Inf,
        filter_gc = F) # filter iterations with GC later
    }
  )
  return (parse_bench_output(res) %>% mutate(pars))
}

bench_tree <- function(model, stumps = c(TRUE, FALSE), 
                       methods, Kvals, pvals, nvals, 
                       nrep, niter, seed_r, seed_cpp, 
                       path, filename_head, bench_quiet = TRUE) {
  rlang::local_options(bench.press_quiet = bench_quiet)

  # benchmarks are more sensible when parameters are in descending order? memory allocation?
  pars_grid <- expand.grid(K = sort(Kvals, decreasing = T), 
                           p = sort(pvals, decreasing = T), 
                           n = sort(nvals, decreasing = T))
  
  t0 <- Sys.time()
  for (i in 1:nrow(pars_grid)) {
    pars <- pars_grid[i,]
    pars_str <- paste(sapply(pars, as.character), collapse = " ")
    
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), model, stumps, 
        "\t", i, "of", nrow(pars_grid), "\t", pars_str, "\n")
    gc(); Sys.sleep(1)
    
    res <- bench_par_grid(model = model, 
                          stumps = stumps, 
                          methods = methods, 
                          pars = pars, 
                          nrep = nrep, 
                          niter = niter, 
                          seed_r = seed_r, 
                          seed_cpp = seed_cpp)
    
    filename <- sprintf("%s/%s-%s-nrep%s-niter%s.csv", path, filename_head, 
                        gsub(" ", "-", pars_str), nrep, niter)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    write.csv(res, file = filename, row.names = F)
  }
  
  cat(sprintf("Benchmark complete for %s %s with parameters\n", model, stumps))
  for (j in 1:ncol(pars_grid)) {
    cat("  ", colnames(pars_grid)[j], ":", unique(as.character(pars_grid[,j])), "\n")
  }
  cat("Total elapsed time:", format(difftime(Sys.time(), t0)), "\n")
}


