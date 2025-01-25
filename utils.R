rsampler <- function(n, FUN = rnorm, pars = NULL, seed = NULL, ...) {
  # Example usage:
  # rsample(10, runif, min = -4, pars = list(max = 4))
  # rsample(10, rpois, lambda = 5)
  # rsample(10, rmultinom, size = 1, prob = rep(1, 5), seed = 1)
  
  # Input validation
  stopifnot(is.function(FUN))
  stopifnot("pars must NULL or a list of parameters to be passed to FUN"
            ={is.null(pars) || is.list(pars)})
  if (!is.null(seed)) set.seed(seed)
  
  # Check for parameter collisions between pars and dot
  pars_dot <- list(...)
  if (!is.null(pars) && length(pars_dot) > 0) {
    # TODO: warn/ignore/error if "n" appears in pars list or ... args
    collisions <- intersect(names(pars), names(pars_dot))
    if (length(collisions) > 0) {
      stop(sprintf("Parameter %s is named in both pars and ...", 
                   paste(collisions, collapse = ", ")))
    }
  }
  pars_all <- c(pars, pars_dot)
  stopifnot("pars must be a named list of parameters"={names(pars_all) != ""})
  
  if (length(pars_all) == 0) {
    pars_all <- formals(FUN)
  } else {
    args_FUN <- formalArgs(FUN)
    
    # Check for unused named parameters and warn
    pars_unused <- setdiff(names(pars_all), args_FUN)
    
    if (length(pars_unused) > 0) {
      warning(sprintf("Unused parameter(s) %s", 
                      paste(pars_unused, collapse = ", ")))
      # Remove unused parameters
      pars_all <- pars_all[intersect(names(pars_all), args_FUN)]
    }
    pars_all <- modifyList(formals(FUN), pars_all)
  }

  make_args <- function(.n) modifyList(as.list(pars_all), list(n = .n))
  FUN_sampler <- function(.n) do.call(FUN, args = make_args(.n))
  samples <- FUN_sampler(n)
  
  pars_sampler <- list(FUN = FUN_sampler, 
                       pars = make_args(n))
  return (structure(samples, sampler = pars_sampler, seed = seed))
}
