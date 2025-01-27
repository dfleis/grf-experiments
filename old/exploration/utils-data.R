library(mvtnorm)
library(grf)

source("utils.R")

GLOBAL_TYPES_ENV <- new.env()
GLOBAL_TYPES_ENV$MODEL_TYPES <- c("hte", "vcm")
lockEnvironment(GLOBAL_TYPES_ENV, bindings = TRUE)
.get_model_types <- function() GLOBAL_TYPES_ENV$MODEL_TYPES


validate_model_type <- function(model_type) {
  match.arg(model_type, .get_model_types())
}


get_FUN_grf <- function(model_type) {
  model_type <- validate_model_type(model_type)
  FUN <- if (model_type == "hte") {
    grf::multi_arm_causal_forest
  } else if (model_type == "vcm") {
    grf::lm_forest
  } else {
    stop("Something went wrong. Invalid model_type in grf_forest")
  }
  return (FUN)
}




#----- Effect functions system
set_effect_labels <- function(effect_funs) {
  K <- length(effect_funs)
  effect_names <- paste0("theta", 1:K)
  effect_labs <- paste0("theta[", 1:K, "](x)")

  effect_funs <- lapply(1:K, function(k) structure(effect_funs[[k]], label = effect_labs[k]))
  names(effect_funs) <- effect_names
  return (effect_funs)
}
get_effect_labels <- function(effect_funs) {
  sapply(effect_funs, attr, "label")
}

create_linear_effect_funs <- function(K, model_type, sampler, ...) {
  model_type <- validate_model_type(model_type)
  
  if (missing(sampler)) sampler <- rnorm#runif
  beta <- rsampler(n = K, FUN = sampler, ...)
  
  if (model_type == "hte") beta[1] <- 0
  
  effect_funs <- vector(mode = "list", length = K)
  for (k in 1:K) {
    effect_funs[[k]] <- local({
      kk <- k
      function(x) beta[kk] * x[,1]
    })
  }
  effect_funs <- set_effect_labels(effect_funs)
  
  structure(effect_funs, model_type = model_type, beta = beta, pars = attributes(beta))
}

validate_effect_funs <- function(K = NULL, effect_funs = NULL, 
                                 model_type, warn = TRUE, ...) {
  if (is.null(K) & is.numeric(effect_funs)) {
    K <- effect_funs
    effect_funs <- NULL
  }
  if (is.null(effect_funs)) {
    effect_funs <- create_linear_effect_funs(K = K, model_type = model_type, ...) 
  } else {
    if (!is.list(effect_funs) || !all(sapply(effect_funs, is.function))) {
      stop("effect_funs must be a list of functions.")
    } else {
      if (!is.null(K) & warn)
        warning("Valid effect_funs supplied, ignoring argument K.")
    }
    # TODO: Do we want to forcibly set the first effect/treatment function
    # to the zero function (with the appropriate dimension) if the target
    # is heterogeneous treatment effects/contrasts?
    # if (model_type == "hte") {
    #   fun1 <- effect_funs[[1]]
    #   effect_funs[[1]] <- function(x) 0 * fun1(x)
    # }
  }
  return (effect_funs)
}


#----- Generate data from a varying coefficient model
generate_X <- function(n, p) {
  X <- matrix(runif(n * p, 0, 1), nrow = n)
  colnames(X) <- paste0("X", 1:p)
  return (X)
}
generate_W_hte <- function(n, K, probs = NULL, ...) {
  # Regressors W for heterogeneous treatment effect models.
  # For HTE models, W is a K-level factor that represents 
  # assignment to discretely many treatment levels
  
  if (is.null(probs)) probs <- rep(1, K)/K
  
  sample(factor(1:K), size = n, replace = TRUE, prob = probs)
}
generate_W_vcm <- function(n, K, mu = NULL, SIGMA = NULL, ...) {
  # Regressors W for (general) varying coefficient models. 
  # For VCM models, W is a K-dimensional random variable,
  # as in, e.g., a set of K continuous treatments.
  
  if (is.null(mu)) mu <- rep(0, K)
  if (is.null(SIGMA)) SIGMA <- diag(1, nrow = K, ncol = K)
  
  mvtnorm::rmvnorm(n, mean = mu, sigma = SIGMA)
}
generate_W <- function(n, K, model_type, ...) {
  model_type <- validate_model_type(model_type)
  
  W <- if (model_type == "hte") {
    generate_W_hte(n = n, K = K, ...)
  } else if (model_type == "vcm") {
    generate_W_vcm(n = n, K = K, ...)  
  }
  
  return (W)
}

generate_data <- function(n, p, K = NULL, effect_funs = NULL, model_type, 
                          sigma_eps = 1, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  # NOTE: This doesn't verify whether the first function in effect_funs is the zero
  # function, as is assumed for the heterogeneous treatment effect/contrast model
  effect_funs <- validate_effect_funs(K = K, effect_funs = effect_funs, model_type = model_type)
  K <- length(effect_funs)
  
  X <- generate_X(n = n, p = p)
  eps <- rnorm(n, 0, sigma_eps)
  W <- generate_W(n = n, K = K, model_type = model_type, ...)
  Wdummy <- model.matrix(~0+W)
  
  # Varying coefficient/effect functions
  efuns <- set_effect_labels(effect_funs)
  if (model_type == "hte") {
    efuns <- efuns[-1]
    Wdummy <- Wdummy[,-1]
  }
  theta <- function(xmat) sapply(efuns, function(tf) tf(xmat))
  # Effects
  thetaX <- theta(X)
  colnames(thetaX) <- names(efuns)

  Y <- rowSums(Wdummy * thetaX) + eps
  
  pars <- list(effect_funs = effect_funs, 
               n = n, p = p, K = K,
               sigma_eps = sigma_eps,
               model_type = model_type)
  
  list(X = X, Y = Y, W = W, # observable data
       thetaX = thetaX, # unobservable true/underlying effects
       theta = theta,   # unobservable true/underlying effect functions
       effect_labels = get_effect_labels(efuns), 
       pars = pars)
}
