library(mvtnorm)
library(grf)

GLOBAL_TYPES_ENV <- new.env()
GLOBAL_TYPES_ENV$MODEL_TYPES <- c("hte", "vcm")
lockEnvironment(GLOBAL_TYPES_ENV, bindings = TRUE)
.get_model_types <- function() GLOBAL_TYPES_ENV$MODEL_TYPES

validate_model_type <- function(model_type) {
  match.arg(model_type, .get_model_types())
}

validate_setting <- function(model_type, setting_id) {
  # HTE Setting ID: Treatment effect function ID, treatment probability ID.
  #   HTE 1: 1, 1
  #   HTE 2: 1, 2
  #   HTE 3: 2, 1
  #   HTE 4: 3, 2
  #   HTE 5: 4, 3
  # VCM Setting ID: Effect function ID
  #   VCM 1: 1
  #   VCM 2: 2
  #   VCM 3: 3
  #   VCM 4: 4
  model_type <- validate_model_type(model_type)
  
  effect_type <- validate_effect_type(model_type = model_type, setting_id = setting_id)
  prob_type <- validate_prob_type(model_type = model_type, setting_id = setting_id)
  
  list(model_type = model_type, 
       setting_id = setting_id, 
       effect_type = effect_type, 
       prob_type = prob_type)
}

validate_effect_type <- function(model_type, setting_id) {
  model_type <- validate_model_type(model_type)
  
  effect_type <- switch(
    model_type,
    "vcm" = { 
      if (setting_id %in% c("1", "2", "3", "4")) {
        as.integer(setting_id)
      } else {
        stop("Invalid setting_id for model_type \"vcm\". Found setting_id = ", 
             setting_id, ", but expected one of 1, 2, 3, 4.")
      }
    },
    "hte" = {
      switch(
        setting_id,
        "1" = 1,
        "2" = 1,
        "3" = 2,
        "4" = 3,
        "5" = 4,
        stop("Invalid setting_id for model_type \"hte\". Found setting_id = ", 
             setting_id, ", but expected one of 1, 2, 3, 4, 5.")
      )
    },
    stop("Invalid model_type. Found ", model_type, ", but expected one of \"vcm\" or \"hte\".")
  )
  return (effect_type)
}
validate_prob_type <- function(model_type, setting_id) {
  model_type <- validate_model_type(model_type)
  
  prob_type <- switch(
    model_type,
    "vcm" = NA,
    "hte" = {
      switch(
        setting_id,
        "1" = 1,
        "2" = 2,
        "3" = 1,
        "4" = 2,
        "5" = 3,
        stop("Invalid setting_id for model_type \"hte\". Found setting_id = ", 
             setting_id, ", but expected one of 1, 2, 3, 4, 5.")
      )
    },
    stop("Invalid model_type. Found ", model_type, ", but expected one of \"vcm\" or \"hte\".")
  )
  return (prob_type)
}

get_FUN_grf <- function(model_type) {
  model_type <- validate_model_type(model_type)
  
  # Capture the expression before assignment
  expr <- if (model_type == "hte") {
    quote(multi_arm_causal_forest)
  } else if (model_type == "vcm") {
    quote(lm_forest)
  } else {
    stop("Something went wrong. Invalid model_type in get_FUN_grf.")
  }
  
  list(FUN = eval(expr), name = deparse(expr))
}

#-------------------------------------------------------
#----- Effect functions system
#-------------------------------------------------------
# logistic-type function that appears in GRF for simulations
varsigma <- function(u) 1 + 1/(1 + exp(-20 * (u - 1/3)))

# Random function generator internal g_l(z_l) function parameters
generate_RFG_pars <- function(p) {
  a <- 0.1
  b <- 2.0
  
  r <- rexp(20, rate = 0.5)
  p_l <- floor(2.5 + r)
  p_l <- pmin(p_l, p)
  
  par_list <- vector(mode = "list", length = 20)
  for (i in 1:20) {
    tstMat <- array(rnorm(p_l[i]), dim=c(p_l[i], p_l[i]))
    U <- qr.Q(qr(tstMat))     # generate random orthogonal matrix
    d <- runif(p_l[i], a, b)  # eigenvalues
    D <- diag(d^2, p_l[i], p_l[i])
    V <- U %*% D %*% t(U)     # get V_l
    
    mu <- rnorm(p_l[i], 0, 1) # get mu_l
    oo <- sample(1:p, p_l[i])
    
    par_list[[i]] <- list(V = V, mu = mu, oo = oo)
  }
  
  return (par_list)
}

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

create_effect_funs_type1 <- function(K, model_type) {
  # Effect setting type 1: theta^*_k(x) = beta_{k1} x_1 for beta_{k1} ~ N(0, 1)
  model_type <- validate_model_type(model_type)
  
  beta <- rnorm(n = K)
  if (model_type == "hte") beta[1] <- 0
  
  effect_funs <- vector(mode = "list", length = K)
  for (k in 1:K) {
    effect_funs[[k]] <- local({
      kk <- k
      function(x) beta[kk] * x[,1]
    })
  }
  effect_funs <- set_effect_labels(effect_funs)
  structure(effect_funs, model_type = model_type)
}
create_effect_funs_type2 <- function(K, model_type) {
  # Effect setting type 2: theta^*_k(x) = varsigma(beta_{k1} x_1) varsigma(beta_{k2} x_2) for beta_{k1},beta_{k2} ~ N(0, 1)
  model_type <- validate_model_type(model_type)
  
  beta1 <- rnorm(n = K)
  beta2 <- rnorm(n = K)
  
  effect_funs <- vector(mode = "list", length = K)
  for (k in 1:K) {
    effect_funs[[k]] <- local({
      kk <- k
      function(x) varsigma(beta1[kk] * x[,1]) * varsigma(beta2[kk] * x[,2])
    })
  }
  if (model_type == "hte") effect_funs[[1]] <-
    local({
      function(x) 0 * x[,1]
    })
  
  effect_funs <- set_effect_labels(effect_funs)
  structure(effect_funs, model_type = model_type)
}
create_effect_funs_type3 <- function(K, p, model_type) {
  # Effect setting type 3: theta^*_k(x) = varsigma(beta^T_k x) for beta_k ~ N_p(0, I)
  model_type <- validate_model_type(model_type)
  
  beta <- matrix(rnorm(K * p), nrow = K, ncol = p)

  effect_funs <- vector(mode = "list", length = K)
  for (k in 1:K) {
    effect_funs[[k]] <- local({
      kk <- k
      function(x) varsigma(x %*% beta[kk,])
    })
  }
  if (model_type == "hte") effect_funs[[1]] <-
    local({
      function(x) 0 * x[,1]
    })
  
  effect_funs <- set_effect_labels(effect_funs)
  structure(effect_funs, model_type = model_type)
}
create_effect_funs_type4 <- function(K, p, model_type) {
  # Effect setting type 4: theta^*_k(x) = RFG(x) (random function generator)
  model_type <- validate_model_type(model_type)
  
  RFG_pars <- lapply(1:K, function(k) generate_RFG_pars(p))
  
  effect_funs <- vector(mode = "list", length = K)
  for (k in 1:K) {
    a_l <- runif(20, -1, 1)
    effect_funs[[k]] <- local({
      kk <- k
      function(x) {
        g <- sapply(RFG_pars[[kk]], function(pars) {
          z <- x[, pars$oo, drop=F]
          zmu <- t(z) - pars$mu
          apply(zmu, 2, function(zm) exp(-0.5 * zm %*% pars$V %*% zm))
        })
        return (g %*% a_l)
      }
    })
  }
  if (model_type == "hte") effect_funs[[1]] <-
    local({
      function(x) 0 * x[,1]
    })
  
  effect_funs <- set_effect_labels(effect_funs)
  structure(effect_funs, model_type = model_type)
}

create_effect_funs <- function(effect_type, K, p, model_type) {
  # Effect type 1: theta^*_k(x) = beta_{k1} x_1 for beta_{k1} ~ N(0, 1)
  # Effect type 2: theta^*_k(x) = varsigma(beta_{k1} x_1) varsigma(beta_{k2} x_2) for beta_{k1},beta_{k2} ~ N(0, 1)
  # Effect type 3: theta^*_k(x) = varsigma(beta^T_k x) for beta_k ~ N_p(0, I)
  # Effect type 4: theta^*_k(x) = RFG(x) (random function generator, see "Fx_generator")
  
  model_type <- validate_model_type(model_type)
  
  effect_funs <- list()
  if (effect_type == 1) {
    effect_funs <- create_effect_funs_type1(K = K, model_type = model_type)  
  } else if (effect_type == 2) {
    effect_funs <- create_effect_funs_type2(K = K, model_type = model_type)  
  } else if (effect_type == 3) {
    effect_funs <- create_effect_funs_type3(K = K, p = p, model_type = model_type)
  } else if (effect_type == 4) {
    effect_funs <- create_effect_funs_type4(K = K, p = p, model_type = model_type)
  } else {
    stop("Invalid 'effect_type' passed to create_effect_funs.")
  }
  
  return (effect_funs)
}

#-------------------------------------------------------
#----- Treatment probability system
#-------------------------------------------------------

create_treatment_probs_type1 <- function(x, K) {
  matrix(1, nrow = nrow(x), ncol = K)/K
}
create_treatment_probs_type2 <- function(x, K) {
  cbind(x[,1], replicate(K - 1, 1 - x[,1])/(K - 1))
}
create_treatment_probs_type3 <- function(x, K) {
  p <- ncol(x)
  gamma <- matrix(rnorm(K * p), nrow = K, ncol = p)
  expXgamma <- exp(x %*% t(gamma))
  expXgamma/rowSums(expXgamma)
}

create_treatment_probs <- function(prob_type, K, x) {
  # Prob type 1: pi_k(x) = 1/K
  # Prob type 2: pi_k(x) = x_1 for k = 1, and (1 - x_1)/(K - 1) for k = 2, ..., K (recall x is Uniform(0,1))
  # Prob type 3: pi_k(x) = exp(gamma_k^T x) / sum^K_{j=1} exp(gamma_j^T x)
  
  treatment_probs <- matrix(NA, nrow = nrow(x), ncol = K)
  if (prob_type == 1) {
    treatment_probs <- create_treatment_probs_type1(x = x, K = K)
  } else if (prob_type == 2) {
    treatment_probs <- create_treatment_probs_type2(x = x, K = K)
  } else if (prob_type == 3) {
    treatment_probs <- create_treatment_probs_type3(x = x, K = K)
  } else {
    stop("Invalid 'prob_type' passed to create_treatment_probs.")
  }
  
  return (treatment_probs)
}

#-------------------------------------------------------
#----- Data generation systems
#-------------------------------------------------------
generate_X <- function(n, p) {
  X <- matrix(runif(n * p, 0, 1), nrow = n)
  colnames(X) <- paste0("X", 1:p)
  return (X)
}
generate_W_vcm <- function(n, K, mu = NULL, SIGMA = NULL, ...) {
  # Regressors W for (general) varying coefficient models. 
  # For VCM models, W is a K-dimensional random variable,
  # as in, e.g., a set of K continuous treatments.
  
  if (is.null(mu)) mu <- rep(0, K)
  if (is.null(SIGMA)) SIGMA <- diag(1, nrow = K, ncol = K)
  
  mvtnorm::rmvnorm(n, mean = mu, sigma = SIGMA)
}
generate_W_hte <- function(n, K, probs = NULL, ...) {
  # Regressors W for heterogeneous treatment effect models.
  # For HTE models, W is a K-level factor that represents 
  # assignment to discretely many treatment levels
  
  if (is.null(probs)) {
    probs <- rep(1, K)/K
    W <- sample(factor(1:K), size = n, replace = TRUE, prob = probs)
  } else {
    W <- apply(probs, 1, function(pr) {
      sample(factor(1:K), size = 1, prob = pr)
    })
  }
  return (W)
}

generate_W <- function(n, K, model_type, prob_type = NULL, x, ...) {
  model_type <- validate_model_type(model_type)
  
  if (model_type == "hte") {
    probs <- create_treatment_probs(prob_type = prob_type, K = K, x = x)
    W <- generate_W_hte(n = n, K = K, probs = probs, ...)
  } else { # VCM
    W <- generate_W_vcm(n = n, K = K, ...)  
  }
  
  return (W)
}


generate_data <- function(n, p, K, model_type, 
                          effect_type = NULL, prob_type = NULL, 
                          sigma_eps = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  X <- generate_X(n = n, p = p)
  eps <- rnorm(n, 0, sigma_eps)
  
  # Regressors/treatments
  W <- generate_W(n = n, K = K, model_type = model_type, prob_type = prob_type, x = X)
  Wdummy <- model.matrix(~0+W)
  
  # Effect functions
  effect_funs <- create_effect_funs(effect_type = effect_type, K = K, p = p, model_type = model_type)
  efuns <- set_effect_labels(effect_funs)
  if (model_type == "hte") {
    efuns <- efuns[-1]
    Wdummy <- Wdummy[,-1]
  }
  theta <- function(xmat) sapply(efuns, function(tf) tf(xmat))
  
  # Effects
  thetaX <- theta(X)
  colnames(thetaX) <- names(efuns)
  
  # Response/outcome
  Y <- rowSums(Wdummy * thetaX) + eps
  
  pars <- list(effect_funs = effect_funs, 
               n = n, p = p, K = K,
               sigma_eps = sigma_eps,
               model_type = model_type,
               effect_type = effect_type,
               prob_type = prob_type)
  
  list(X = X, Y = Y, W = W, # observable data
       thetaX = thetaX, # unobservable true/underlying effects
       theta = theta,   # unobservable true/underlying effect functions
       effect_labels = get_effect_labels(efuns), 
       pars = pars)
}


# set.seed(124)
# n <- 10000
# p <- 5
# K <- 4
# model_type <- "hte"
# effect_type <- 4 # 1, 2, 3, 4
# prob_type <- 3 # 1, 2, 3
# 
# data <- generate_data(n = n, p = p, K = K, effect_type = effect_type, prob_type = prob_type)
# data$X
# data$Y
# table(data$W)
