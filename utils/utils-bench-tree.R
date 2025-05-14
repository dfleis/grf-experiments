suppressMessages(library(tidyverse))
library(bench)
source("utils/utils-data.R")

train_grf_wrapper <- function(data, args) {
  # Wrapper for the core C++ function called by GRF for both 
  # grf::multi_arm_causal_forest and grf::lm_forest.
  #
  # WARNING: The underlying C++ function expects the data and arguments to
  # be formatted in a very specific way. For this reason, much of the code
  # below is dedicated to pre-processing the data for this function.
  #
  # The motivation for calling GRF in this way is that we want to time the
  # core GRF algorithm and not the pre-processing steps GRF does prior
  # to calling the core training function.
  
  forest <- grf:::do.call.rcpp(grf:::multi_causal_train, c(data, args))
  
  return((length(forest$`_split_vars`[[1]]) - 1)/2) # return number of splits for validation
}

#-------------------------------------------------------
#----- Pre-processing tools
#-------------------------------------------------------
prep_args_common <- function(X, method, stump, seed) {
  # GRF doesn't offer a maximum depth argument when fitting its trees.
  # Instead, GRF controls its stopping criteria/tree complexity via the 
  # arguments:
  #   * min.node.size: (target) minimum number of samples in a node (default = 5
  #     for multi_arm_causal_forest and lm_forest). This is only a "target" minimum
  #     node size, similar to randomForest, such that we can observe nodes with
  #     even fewer numbers of samples in a trained tree.
  #   * alpha: imbalance proportion of a split (default = 0.05).
  #   * imbalance.penalty: the degree to which imbalance penalized  (default = 0).
  #
  # Therefore, we can force GRF to fit stumps by setting alpha = 0 and
  # min.node.size = "the number of observations passed to each tree" - 1.
  # Specifically, if we disable subsampling and honest splitting, we can force
  # GRF to fit stumps by setting alpha = 0 and min.node.size = nrow(X) - 1.
  #
  # NOTE: This is only guaranteed to work when centering is disabled (Y.hat = 0, W.hat = 0),
  # otherwise GRF will tend to fit trees with 0 splits (i.e. only root nodes). 
  # This is fine for timing tests since we're only interested in comparing the core 
  # algorithmic differences between GRF-grad and GRF-FPT, and not the common
  # centering calculations.
  method <- match.arg(tolower(method), choices = c("grad", "fpt1", "fpt2"))
  
  if (isTRUE(stump)) min.node.size = nrow(X) - 1 # force GRF to fit stumps
  else min.node.size = 5L
  alpha = 0.0
  
  # Other arguments we modify for single tree/stump timing tests
  num.trees = 1L
  mtry = ncol(X)
  sample.fraction = 1.0
  honesty = FALSE
  honesty.fraction = 1.0
  honesty.prune.leaves = FALSE
  ci.group.size = 1L
  compute.oob.predictions = FALSE
  num.threads = 1L # only relevant when num.trees > 1
  
  # Remaining default arguments
  sample.weights = NULL
  clusters = numeric(0)
  equalize.cluster.weights = FALSE
  imbalance.penalty = 0
  stabilize.splits = FALSE
  
  # Process arguments (internal steps done in grf::lm_forest and grf::multi_arm_causal_forest)
  num.threads <- grf:::validate_num_threads(num.threads)
  samples.per.cluster <- grf:::validate_equalize_cluster_weights(
    equalize.cluster.weights, clusters, sample.weights)
  method.flag <- switch(method, grad = 1, fpt1 = 2, fpt2 = 3)
  
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
               method.flag = method.flag,
               num.threads = num.threads,
               seed = seed, 
               legacy.seed = grf:::get_legacy_seed(),
               sample.weights = sample.weights)
  return (args)
}

prep_data_common <- function(data, method, stump, seed) {
  args <- prep_args_common(X = data$X, method = method, stump = stump, seed = seed)
  
  W <- stats::model.matrix(~0 + data$W)
  Y <- grf:::validate_observations(data$Y, data$X, allow.matrix = TRUE)
  
  # Disable centering for single split/tree timing tests
  W.centered <- W # W - data$W.hat
  Y.centered <- Y # Y - data$Y.hat
  
  return(list(X = data$X, W.centered = W.centered, Y.centered = Y.centered, args = args))
}

prep_data_hte <- function(data, method, stump, seed) {
  # Internal pre-processing steps specific to grf::multi_arm_causal_forest.
  prep_init <- prep_data_common(data = data, method = method, stump = stump, seed = seed)
  prep_init$args$gradient.weights <- numeric(0)
  
  data <- grf:::create_train_matrices(X = prep_init$X, 
                                      outcome = prep_init$Y.centered, 
                                      treatment = prep_init$W.centered[,-1], 
                                      sample.weights = prep_init$args$sample.weights)
  prep_init$args$sample.weights <- NULL
  return (list(data = data, args = prep_init$args))
}

prep_data_vcm <- function(data, method, stump, seed) {
  # Internal pre-processing steps specific to grf::lm_forest.
  prep_init <- prep_data_common(data = data, method = method, stump = stump, seed = seed)
  prep_init$args$gradient.weights <- rep(1, NCOL(data$W) * NCOL(data$Y))
  
  data <- grf:::create_train_matrices(X = prep_init$X, 
                                      outcome = prep_init$Y.centered, 
                                      treatment = prep_init$W.centered, 
                                      sample.weights = prep_init$args$sample.weights)
  prep_init$args$sample.weights <- NULL
  return (list(data = data, args = prep_init$args))
}

prep_data <- function(data_args, method, stump, forest_seed) {
  # TODO: Add option to enable centering rather than force disabling it?
  # Might be useful for forest tests (multiple trees), but if we're just
  # doing timing tests then it makes more sense to just disable centering.
  model_type <- validate_model_type(data_args$model_type)
  data <- do.call(generate_data, data_args)
  prep_args <- list(data = data, method = method, stump = stump, seed = forest_seed)
  
  prep_FUN <- if (model_type == "vcm") {
    prep_data_vcm
  } else if (model_type == "hte") {
    prep_data_hte
  } else {
    # This condition should never hit since we call validate_model_type earlier
    stop("Something went wrong in prep_data. ",
         "Found model_type = ", model_type, ", but expected one of ", 
         "\"vcm\" or \"hte\".")
  }
  
  return (do.call(prep_FUN, prep_args))
}

#-------------------------------------------------------
#----- Post-processing tools
#-------------------------------------------------------
process_bench_tree <- function(bench_data, data_pars) {
  bench_data %>%
    mutate(splits = result) %>%
    unnest_longer(c(gc, time, splits), indices_include = T) %>%
    unnest_wider(gc) %>%
    rename(gc0 = level0, gc1 = level1, gc2 = level2, iter = time_id) %>%
    mutate(gc   = gc0 + gc1 + gc2 > 0,
           time = as.numeric(time)) %>%
    select(method, rep, iter, gc, gc0, gc1, gc2, splits, time) %>%
    mutate(data_pars)
}

#-------------------------------------------------------
#----- Benchmarking tools
#-------------------------------------------------------
bench_tree_pars <- function(methods, model_type, stump, 
                            data_pars, nrep, niter, 
                            seed = NULL, .quiet = TRUE) {
  # model_type: "vcm", "hte"
  # methods: "grad", "fpt1", "fpt2"
  # stump: TRUE, FALSE
  # data_pars: setting_id, K, p, n
  # nrep: bench replications (see ?bench::mark and ?bench::press)
  # niter: bench iterations (see ?bench::mark and ?bench::press)
  res <- bench::press(
    rep = 1:nrep, # dummy variable to do replications over new data samples
    method = methods,
    {
      data_args <- data_pars
      data_args$seed <- if (!is.null(seed)) seed + rep - 1 else rep - 1
      data_args$model_type <- model_type
      
      prepped <- prep_data(data_args = data_args, 
                           method = method, 
                           stump = stump, 
                           forest_seed = seed) 

      bench::mark(
        train_grf_wrapper(data = prepped$data, args = prepped$args),
        iterations = niter, 
        min_time = Inf,
        filter_gc = FALSE) # filter iterations with GC later
    },
    .quiet = .quiet
  )

  pars <- bind_cols(stump = stump, model_type = model_type, data_pars)
  return (process_bench_tree(res, pars))
}

#--------------------------------------------------
#----- Core bench function
#--------------------------------------------------
bench_tree <- function(methods, model_type, setting_id, 
                       stumps = c(TRUE, FALSE), Kvals, pvals, nvals, 
                       nrep, niter, seed = NULL,
                       path = "data", filename_head = "bench-tree", 
                       .quiet = TRUE) {
  setting_id <- validate_setting(model_type, setting_id)$setting_id
  
  pars_grid <- expand.grid(
    K = sort(Kvals, decreasing = T),
    p = sort(pvals, decreasing = T), 
    n = sort(nvals, decreasing = T),
    stump = stumps,
    setting_id = setting_id
  )
  
  file_dir <- sprintf("%s/%s-%s-%s", path, filename_head, nrep, niter)
  dir.create(file_dir, recursive = TRUE, showWarnings = FALSE)
  
  t0 <- Sys.time()
  for (i in 1:nrow(pars_grid)) {

    pars <- pars_grid[i,]
    pars_str <- paste(c(model_type, setting_id, pars$stump, as.integer(pars$K), 
                        as.integer(pars$p), as.integer(pars$n)), collapse = "-")
    # file_dir/filename_head-[NREPS]-[NITERS]-[MODEL]-[SETTING]-[STUMP]-[K]-[p]-[n].csv
    filename_str <- sprintf("%s/%s-%s-%s-%s.csv", file_dir, filename_head, 
                            nrep, niter, pars_str)
    
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "::", 
        i, "of", nrow(pars_grid), "::", filename_str, "\n")
    
    #gc(); Sys.sleep(1)
    
    res <- bench_tree_pars(
      methods = methods, 
      model_type = model_type, 
      stump = pars$stump,
      data_pars = modifyList(pars, list(stump = NULL)), 
      nrep = nrep, 
      niter = niter, 
      seed = seed,
      .quiet = .quiet
    )
    
    write.csv(res, file = filename_str, row.names = FALSE)
  }
  
  cat(sprintf("Benchmark complete for model_type = %s with\n", model_type))
  for (j in 1:ncol(pars_grid)) {
    cat("  ", colnames(pars_grid)[j], ":", unique(as.character(pars_grid[,j])), "\n")
  }
  cat(sprintf("Benchmark replications nrep = %s\n", nrep))
  cat(sprintf("Benchmark iterations niter = %s\n", niter))
  cat("Total elapsed time:", format(difftime(Sys.time(), t0)), "\n")
}

