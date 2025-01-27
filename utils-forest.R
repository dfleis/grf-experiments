library(grf)

center_forests <- function(data, model_type, 
                           disable_centering = FALSE, 
                           num_trees = 2000, min_node_size = 5, ...) {
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
    
    forest_W <- if (model_type == "vcm") {
      grf::multi_regression_forest(X = X, Y = W, 
                                   num.trees = num_trees_centering, 
                                   min.node.size = min_node_size, 
                                   compute.oob.predictions = TRUE, ...)
    } else { # hte
      grf::probability_forest(X = X, Y = W, 
                              num.trees = num_trees_centering,
                              min.node.size = min_node_size,
                              compute.oob.predictions = TRUE, ...)
    }
    
    Y_hat <- predict(forest_Y)$predictions
    W_hat <- predict(forest_W)$predictions
    
    if (isTRUE(disable_centering)) {
      Y_hat <- 0 * Y_hat
      W_hat <- 0 * W_hat
    }
    return (list(Y_hat = Y_hat, W_hat = W_hat))
  })
}