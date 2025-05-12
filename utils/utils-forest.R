library(grf)
source("utils/utils-data.R")

get_FUN_grf <- function(model_type) {
  model_type <- validate_model_type(model_type)
  
  # Capture the expression before assignment
  expr <- if (model_type == "hte") {
    quote(multi_arm_causal_forest)
  } else if (model_type == "vcm") {
    quote(lm_forest)
  } else { 
    # This condition should never hit since we call validate_model_type earlier
    stop("Something went wrong in get_FUN_grf. ",
         "Found model_type = ", model_type, ", but expected one of ", 
         "\"vcm\" or \"hte\".")
  }
  list(FUN = eval(expr), name = deparse(expr))
}

make_centering_args <- function(model_type, X, args) {
  arg_names <- c("num.trees", "sample.weights", "clusters", "equalize.cluster.weights", 
                 "sample.fraction", "mtry", "min.node.size", "honesty", "honesty.fraction",
                 "honesty.prune.leaves", "alpha", "imbalance.penalty", "num.threads", "seed")
  
  model_type <- validate_model_type(model_type)
  FUN_grf <- get_FUN_grf(model_type)
  
  args <- args[arg_names]
  default_args <- as.list(formals(FUN_grf$FUN))[arg_names]
  
  args$num.trees <- max(50, args$num.trees/4)
  args$min.node.size <- 5
  args$honesty <- TRUE
  args$honesty.fraction <- 0.5
  args$compute.oob.predictions <- TRUE
  args$X <- X
  
  centering_args <- modifyList(default_args, args)
  centering_args$mtry <- eval(centering_args$mtry, envir = environment(list(X = X)))
  return (centering_args)
}

make_centering_data <- function(data, model_type, center_data = TRUE, grf_args = NULL) {
  # By default, GRF uses a Robinson-style centering step and fits its forest based on
  # centered values Y - Y.hat and W - W.hat.
  # 
  # Y.hat: 
  #   * Both VCM and HTE forests (lm_forest and multi_arm_causal_forest, respectively)
  #     center the response/outcome Y by training a multi_regression_forest on the 
  #     model Y ~ X, i.e.
  #       forest.Y <- multi_regression_forest(X = X, Y = Y, ...)
  #   * Following this, Y.hat is computed via in-sample predictions
  #       Y.hat <- predict(forest.Y)
  #   * Note multi_regression_forest is equivalent to regression_forest when Y is a
  #     scalar response, but not numerically identical (grf uses multi_regression_forest
  #     for centering in the case of both scalar and vector responses)
  #
  # W.hat: 
  #   * VCM forests (lm_forest) centers the regressors by training a multi_regression_forest
  #     on the model W ~ X, i.e.
  #       forest.W <- multi_regression_forest(X = X, Y = W, ...)
  #   * HTE forests (multi_arm_causal_forest) centers the regressors (treatment indicators) by 
  #     training a probability_forest on the model W ~ X, i.e.
  #       forest.W <- probability_forest(X = X, Y = W, ...)
  #   * Following this, W.hat is computed via in-sample predictions
  #       W.hat <- predict(forest.W)
  #
  # GRF-grad and GRF-FPT are identical with respect to regression_forest, 
  # multi_regression_forest, and probability_forest since the form of the
  # pseudo-outcomes are the same (up to a scalar multiple). Therefore, we
  # want to make sure GRF-grad and GRF-FPT see the same centered data 
  # Y - Y.hat and W - W.hat (since we're only interested in comparing the
  # ways in which GRF differs under the choice of grad/FPT). 
  model_type <- validate_model_type(model_type)
  
  FUN_forest_W <- if (model_type == "vcm") {
    multi_regression_forest
  } else if (model_type == "hte") {
    probability_forest
  } else {
    # This condition should never hit since we call validate_model_type earlier
    stop("Something went wrong in center_forest. ",
         "Found model_type = ", model_type, ", but expected one of ", 
         "\"vcm\" or \"hte\".")
  }
  
  out <- if (isFALSE(center_data)) {
    Y.hat <- matrix(0, nrow = NROW(data$Y), ncol = NCOL(data$Y))
    W.hat <- 0 * stats::model.matrix(~0+data$W)
    if (is.factor(data$W)) colnames(W.hat) <- levels(data$W)
    
    list(Y.hat = Y.hat, W.hat = W.hat)
  } else {
    args.orthog <- make_centering_args(model_type = model_type, X = data$X, args = grf_args)
    
    forest.Y <- do.call(multi_regression_forest, modifyList(args.orthog, list(Y = data$Y)))
    
    if (model_type == "hte") args.orthog <- modifyList(args.orthog, list(ci.group.size = 1))
    forest.W <- do.call(FUN_forest_W, modifyList(args.orthog, list(Y = data$W)))
    
    Y.hat <- predict(forest.Y)$predictions
    W.hat <- predict(forest.W)$predictions
    
    list(Y.hat = Y.hat, W.hat = W.hat)
  }
  out$X <- data$X
  out$Y <- data$Y
  out$W <- data$W
  
  return (out)
}

