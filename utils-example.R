suppressMessages(library(tidyverse))
library(grf)

source("utils-data.R")

generate_data_example <- function(n_test, n_train, p, effect_funs, model_type, ...) {
  # Training data
  train <- generate_data(n = n_train, p = p, 
                         effect_funs = effect_funs, 
                         model_type = model_type, ...)
  
  # Test data, equally space along X_1 (visualization purposes)
  x_range <- c(min(train$X[,1]), max(train$X[,1]))
  X_test <- 0 * generate_X(n = n_test, p = p)
  X_test[,1] <- seq(x_range[1], x_range[2], length.out = n_test)
  
  pars <- list(model_type = model_type)
  
  list(train = train, test = list(X = X_test, 
                                  thetaX = train$theta(X_test)))
}

make_preds_example <- function(train, test, method = c("grad", "fpt1", "fpt2"), args_grf) {
  # Specify model
  method <- match.arg(method)
  FUN_grf <- get_FUN_grf(train$pars$model_type)
  args <- modifyList(args_grf, list(X = train$X, Y = train$Y, W = train$W, method = method))
  
  if (isTRUE(args[["Y.hat"]] == 0)) args <- modifyList(args, list(Y.hat = 0 * train$Y))
  if (isTRUE(args[["W.hat"]] == 0)) {
    if (identical(FUN_grf, grf::multi_arm_causal_forest)) {
      W.hat <- 0 * model.matrix(~0+train$W)
      colnames(W.hat) <- levels(train$W)
    } else {
      W.hat <- 0 * train$W
    }
    args <- modifyList(args, list(W.hat = W.hat))
  }
  
  # Train forest
  forest <- do.call(FUN_grf, args)
  
  # Make new estimates
  est_var <- isFALSE(forest$ci.group.size == 1)
  preds <- predict(forest, newdata = test$X, estimate.variance = est_var, drop = TRUE)
  
  # Restructure
  est <- preds$predictions
  se <- if (!is.null(preds$variance.estimates)) { 
    sqrt(preds$variance.estimates)
  } else {
    matrix(NA, nrow = nrow(est), ncol = ncol(est))
  } 

  res <- list(truth = test$thetaX, 
              est = est, 
              se = se) %>%
    discard(is.null) %>%
    lapply(function(mat) {
      colnames(mat) <- colnames(test$thetaX)
      mat
    })
  
  #list(preds = res, X = test$X, method = method, effect_labels = colnames(test$thetaX))
  list(preds = res, X = test$X, method = method, effect_labels = train$effect_labels)
}

process_preds_example <- function(preds_data) {
  method <- attr(preds_data, "method")
  effect_labels <- preds_data$effect_labels

  preds_data$preds %>%
    lapply(bind_cols, data.frame(preds_data$X)) %>%
    imap(~.x %>% mutate(stat_type = .y)) %>%
    bind_rows() %>%
    pivot_longer(cols = any_of(names(effect_labels)),
                 names_to = "effect",
                 values_to = "value") %>%
    mutate(method = method) %>%
    mutate(effect = factor(effect, labels = effect_labels))
}

run_example <- function(n_test, n_train, p, effect_funs = NULL, model_type, args_grf = NULL, 
                        methods = NULL, ci_level = 0.95, ...) {
  if (!is.null(args_grf[["method"]])) methods <- unique(c(methods, args_grf[["method"]]))
  if (is.null(methods)) methods <- c("grad", "fpt1", "fpt2")
  
  data <- generate_data_example(n_test = n_test, n_train = n_train, p = p, 
                                effect_funs = effect_funs, model_type = model_type, ...)

  res <- lapply(methods, function(method) {
    preds_data <- make_preds_example(train = data$train, test = data$test, 
                                     method = method, args_grf = args_grf)
    process_preds_example(preds_data) %>%
      mutate(method = method)
  }) %>% bind_rows() %>%
    mutate(method = factor(method, levels = methods)) %>%
    mutate(method = recode(method, fpt1 = "FPT1", fpt2 = "FPT2"))
  
  df_list <- split(res, res$stat_type)

  Z_CI <- qnorm(1 - (1 - ci_level)/2)
  
  df_list$est$ci_lo <- df_list$est$value - Z_CI * df_list$se$value
  df_list$est$ci_hi <- df_list$est$value + Z_CI * df_list$se$value

  #df_list$train <- data$train
  
  structure(df_list, ci_level = ci_level)
}

