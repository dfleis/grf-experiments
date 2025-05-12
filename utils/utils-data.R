library(mvtnorm) 
library(grf)
library(rfg)

#-------------------------------------------------------
#----- Model systems
#-------------------------------------------------------
GLOBAL_TYPES_ENV <- new.env()
GLOBAL_TYPES_ENV$MODEL_TYPES <- c("hte", "vcm")
lockEnvironment(GLOBAL_TYPES_ENV, bindings = TRUE)
.get_model_types <- function() GLOBAL_TYPES_ENV$MODEL_TYPES

validate_model_type <- function(model_type) {
  match.arg(model_type, .get_model_types())
}

#-------------------------------------------------------
#----- Data-generating settings systems (validation)
#-------------------------------------------------------
validate_setting <- function(model_type, setting_id) {
  # HTE Setting ID: Treatment effect function ID, treatment probability ID.
  #   HTE 1: Effect 1, probability 1 -- linear, dense, unconfounded.
  #   HTE 2: Effect 1, probability 2 -- linear, dense, confounded.
  #   HTE 3: Effect 2, probability 1 -- nonlinear, sparse when dim(X) > 2, unconfounded.
  #   HTE 4: Effect 3, probability 2 -- nonlinear, dense, confounded.
  #   HTE 5: Effect 4, probability3 -- RFG.
  # VCM Setting ID: Effect function ID
  #   VCM 1: Effect 1 -- linear, dense.
  #   VCM 2: Effect 2 -- nonlinear, sparse when dim(X) > 2.
  #   VCM 3: Effect 3 -- nonlinear, dense.
  #   VCM 4: Effect 4 -- RFG.
  model_type <- validate_model_type(model_type)
  
  effect_type <- setting_to_effect_type(model_type = model_type, setting_id = setting_id)
  prob_type <- setting_to_prob_type(model_type = model_type, setting_id = setting_id)
  
  list(model_type = model_type, 
       setting_id = setting_id, 
       effect_type = effect_type, 
       prob_type = prob_type)
}

setting_to_effect_type <- function(model_type, setting_id) {
  model_type <- validate_model_type(model_type)
  setting_id <- as.character(setting_id)
  
  effect_type <- switch(
    model_type,
    "vcm" = { 
      if (setting_id %in% c("1", "2", "3", "4")) {
        as.integer(setting_id)
      } else {
        stop("Invalid setting_id for model_type \"vcm\". Found setting_id = ", 
             setting_id, ", but expected one of 1, ..., 4.")
      }
    },
    "hte" = {
      switch(
        setting_id,
        "1" = "1", # linear sparse
        "2" = "1", # linear sparse
        "3" = "2", # logistic interaction
        "4" = "3", # logistic dense
        "5" = "4", # RFG
        stop("Invalid setting_id for model_type \"hte\". ",
             "Found setting_id = ", setting_id,
             ", but expected one of 1, ..., 5."))
    }
  )
  return (effect_type)
}
setting_to_prob_type <- function(model_type, setting_id) {
  model_type <- validate_model_type(model_type)
  setting_id <- as.character(setting_id)
  
  prob_type <- switch(
    model_type,
    "vcm" = NA,
    "hte" = {
      switch(
        setting_id,
        "1" = "1", # uniform
        "2" = "2", # linear
        "3" = "1", # uniform
        "4" = "2", # linear
        "5" = "3", # softmax
        stop("Invalid setting_id for model_type \"hte\". ",
             "Found setting_id = ", setting_id, 
             ", but expected one of 1, 2, 3, 4, 5."))
    }
  )
  return (prob_type)
}


#-------------------------------------------------------
#----- Effect functions systems
#-------------------------------------------------------
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

# Logistic-type function that appears in the causal forest simulations of GRF
varsigma <- function(u) 1 + 1/(1 + exp(-20 * (u - 1/3)))

# Effect function templates
effect_linear <- function(x, b) b * x[,1]
effect_logistic_interaction <- function(x, b) varsigma(b[1] * x[,1]) * varsigma(b[2] * x[,2])
effect_logistic_dense <- function(x, b) varsigma(x %*% b)
effect_RFG <- function(x, RFG_FUN) RFG_FUN(x)

create_effect_funs_factory <- function(K, par_list, effect_template, model_type) {
  model_type <- validate_model_type(model_type)
  stopifnot(length(par_list) == K)
  effect_funs <- vector(mode = "list", length = K)
  for (k in 1:K) {
    effect_funs[[k]] <- local({
      kk <- k
      function(x) effect_template(x, par_list[[kk]])
    })
  }
  if (model_type == "hte") 
    effect_funs[[1]] <- local(function(x) 0 * x[,1])
  
  effect_funs <- set_effect_labels(effect_funs)
  structure(effect_funs, model_type = model_type)
}

create_effect_funs_linear <- function(K, model_type) {
  # Effect setting type 1: theta^*_k(x) = beta_{k1} x_1 for beta_{k1} ~ N(0, 1)
  # Linear, sparse.
  beta <- rnorm(n = K)
  create_effect_funs_factory(K = K, par_list = beta, 
                             effect_template = effect_linear, model_type = model_type)
}
create_effect_funs_logistic_interaction <- function(K, model_type) {
  # Effect setting type 2: theta^*_k(x) = varsigma(beta_{k1} x_1) varsigma(beta_{k2} x_2) for beta_{k1},beta_{k2} ~ N(0, 1)
  # Logistic, interaction, sparse.
  betas <- replicate(K, rnorm(2), simplify = FALSE)
  create_effect_funs_factory(K = K, par_list = betas, 
                             effect_template = effect_logistic_interaction, model_type = model_type)
}
create_effect_funs_logistic_dense <- function(K, p, model_type) {
  # Effect setting type 3: theta^*_k(x) = varsigma(beta^T_k x) for beta_k ~ N_p(0, I)
  # Logistic, dense.
  betas <- replicate(K, rnorm(p), simplify = FALSE)
  create_effect_funs_factory(K = K, par_list = betas, 
                             effect_template = effect_logistic_dense, model_type = model_type)
}
create_effect_funs_RFG <- function(K, p, model_type) {
  # Effect setting type 4: theta^*_k(x) = RFG(x)
  # Random function generator.
  #RFG_pars <- replicate(K, generate_RFG_pars(p), simplify = FALSE)
  RFG_FUNS <- replicate(K, rfg::rfg(p))
  create_effect_funs_factory(K = K, par_list = RFG_FUNS, 
                             effect_template = effect_RFG, model_type = model_type)
}

create_effect_funs <- function(effect_type, K, p, model_type) {
  # Effect type 1: theta^*_k(x) = beta_{k1} x_1 for beta_{k1} ~ N(0, 1)
  # Effect type 2: theta^*_k(x) = varsigma(beta_{k1} x_1) varsigma(beta_{k2} x_2) for beta_{k1},beta_{k2} ~ N(0, 1)
  # Effect type 3: theta^*_k(x) = varsigma(beta^T_k x) for beta_k ~ N_p(0, I)
  # Effect type 4: theta^*_k(x) = RFG(x)
  effect_type <- as.character(effect_type)
  effect_funs <- switch(effect_type, 
         "1" = create_effect_funs_linear(K = K, model_type = model_type),
         "2" = create_effect_funs_logistic_interaction(K = K, model_type = model_type),
         "3" = create_effect_funs_logistic_dense(K = K, p = p, model_type = model_type),
         "4" = create_effect_funs_RFG(K = K, p = p, model_type = model_type),
         stop("Invalid effect_type passed to create_effect_funs."))
  return (effect_funs)
}

#-------------------------------------------------------
#----- Treatment probability system
#-------------------------------------------------------
# Treatment probability function pi(x) = (pi_1(x), ..., pi_K(x)) templates
prob_uniform <- function(x, K, ...) rep(1, K)/K
prob_linear <- function(x, K, ...) c(x[1], rep(1 - x[1], times = K - 1)/(K - 1))
prob_softmax <- function(x, K, gamma) {
  expXgamma <- exp(x %*% t(gamma))
  expXgamma/sum(expXgamma)
}

create_prob_funs_factory <- function(K, pars, prob_template) {
  prob_funs <- function(x) {
    t(apply(x, 1, FUN = function(xx) prob_template(x = xx, K = K, pars)))
  }
  return (prob_funs)
}

create_prob_funs_uniform <- function(K) {
  create_prob_funs_factory(K = K, prob_template = prob_uniform)
}
create_prob_funs_linear <- function(K) {
  create_prob_funs_factory(K = K, prob_template = prob_linear)
}
create_prob_funs_softmax <- function(K, p) {
  gamma <- matrix(rnorm(K * p), nrow = K)
  create_prob_funs_factory(K = K, pars = gamma, prob_template = prob_softmax)
}

create_probs_funs <- function(prob_type, K, p) {
  # Probability type 1: pi_k(x) = 1/K (uniform)
  # Probability type 2: pi_k(x) = x_1 for k = 1, and (1 - x_1)/(K - 1) for k = 2, ..., K (we assume x ~ U([0,1]^p))
  # Probability type 3: pi_k(x) = exp(gamma_k^T x) / sum^K_{j=1} exp(gamma_j^T x) for gamma_j ~ N_K(0, I)
  prob_type <- as.character(prob_type)
  prob_funs <- switch(prob_type,
                      "1" = create_prob_funs_uniform(K = K),
                      "2" = create_prob_funs_linear(K = K),
                      "3" = create_prob_funs_softmax(K = K, p = p),
                      stop("Invalid prob_type passed to create_prob_funs."))
  return (prob_funs)
}

#-------------------------------------------------------
#----- Data generation systems
#-------------------------------------------------------
generate_X <- function(n, p, theta = 0.3) {
  #--- Generate correlated normal quantiles Z ~ N_p(0, \Theta) where \Theta_{i,j} = \theta^{|i - j|}
  #--- then form covariates X from the normal copula such that X = (F(Z_1), F(Z_2), ..., F(Z_p))
  #X <- matrix(runif(n * p, 0, 1), nrow = n)
  Theta <- theta^abs(outer(1:p, 1:p, function(i, j) {i - j}))
  X <- apply(mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Theta), 2, pnorm)
  colnames(X) <- paste0("X", 1:p)
  return (X)
}

generate_W_vcm <- function(n, K, x, mu = NULL, SIGMA = NULL, ...) {
  # Regressors W for (general) varying coefficient models. 
  # For VCM models, W is a K-dimensional random variable,
  # as in, e.g., a set of K continuous treatments.
  if (is.null(mu)) mu <- rep(0, K)
  if (is.null(SIGMA)) SIGMA <- diag(1, nrow = K, ncol = K)
  mvtnorm::rmvnorm(n, mean = mu, sigma = SIGMA)
}
generate_W_hte <- function(n, K, x, prob_funs = NULL, ...) {
  # Regressors W for heterogeneous treatment effect models.
  # For HTE models, W is a K-level factor that represents 
  # assignment to one of discretely many treatment levels.
  treatment_levels <- factor(1:K)
  W <- if (is.null(prob_funs)) {
    probs <- rep(1, K)/K
    sample(treatment_levels, size = n, replace = TRUE, prob = probs)
  } else {
    probs <- prob_funs(x)
    apply(probs, 1, function(pr) {
      sample(treatment_levels, size = 1, prob = pr)
    })
  }
  return (W)
}
generate_W <- function(n, K, x, model_type, prob_type = NULL, ...) {
  model_type <- validate_model_type(model_type)
  stopifnot(nrow(x) == n)
  p <- ncol(x)
  
  W <- if (model_type == "hte") {
    prob_funs <- create_probs_funs(prob_type = prob_type, K = K, p = p)
    generate_W_hte(n = n, K = K, x = x, prob_funs = prob_funs, ...)
  } else { # vcm
    generate_W_vcm(n = n, K = K, x = x, ...)  
  }
  return (W)
}

generate_data <- function(n, p, K, model_type, 
                          setting_id = NULL,
                          effect_type = NULL, 
                          prob_type = NULL, 
                          sigma_eps = 1, 
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  model_type <- validate_model_type(model_type)
  if (!is.null(setting_id)) {
    # If setting_id is supplied, then it takes precedence to specify effect_type and prob_type.
    setting_info <- validate_setting(model_type = model_type, setting_id = setting_id)
    effect_type <- setting_info$effect_type
    prob_type <- setting_info$prob_type
  }
  
  # Generate effect functions
  effect_funs <- create_effect_funs(effect_type = effect_type, K = K, p = p, model_type = model_type)
  
  # Generate covariates, noise
  X <- generate_X(n = n, p = p)
  eps <- rnorm(n, 0, sigma_eps)
  
  # Generate regressors/treatments
  W <- generate_W(n = n, K = K, model_type = model_type, prob_type = prob_type, x = X)
  Wdummy <- model.matrix(~0+W)
  
  # Format effect functions
  efuns <- set_effect_labels(effect_funs)
  if (model_type == "hte") {
    efuns <- efuns[-1]
    Wdummy <- Wdummy[,-1]
  }
  theta <- function(xmat) sapply(efuns, function(tf) tf(xmat))
  
  # Compute effects
  thetaX <- theta(X)
  colnames(thetaX) <- names(efuns)
  
  # Compute response/outcome
  Y <- rowSums(Wdummy * thetaX) + eps
  
  pars <- list(
    effect_funs = effect_funs, 
    n = n, p = p, K = K,
    sigma_eps = sigma_eps,
    model_type = model_type,
    effect_type = effect_type,
    prob_type = prob_type
  )
  
  list(
    X = X, Y = Y, W = W, # observable data
    thetaX = thetaX, # unobservable true/underlying effects
    theta = theta,   # unobservable true/underlying effect functions
    effect_labels = get_effect_labels(efuns), 
    pars = pars
  )
}