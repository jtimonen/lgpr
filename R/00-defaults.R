#' Default priors and parameter values
#'
#' @description
#' \itemize{
#'   \item \code{default_prior} returns a named list that defines a prior
#'   \item \code{default_ppc_fun} returns a function to be used as argument
#'   of \code{\link{ppc}}
#' }
#' @param name parameter name
#' @param num_uncrt number of uncertain components
#' @param object an object of class \linkS4class{lgpfit} or \code{lgpmodel}
#' @return see description
#' @name defaults
NULL

#' @rdname defaults
default_prior <- function(name, num_uncrt = NULL) {
  if (name == "effect_time_info") {
    desc <- default_prior_effect_time_info(num_uncrt)
  } else {
    desc <- list(default_prior_common(name))
  }
  return(desc)
}

#' @rdname defaults
default_prior_common <- function(name) {
  allowed <- c(
    "alpha", # 1
    "ell", # 2
    "sigma", # 3
    "phi", # 4
    "wrp", # 5
    "beta", # 6
    "gamma", # 7
    "effect_time" # 8
  )
  idx <- check_allowed(name, allowed)
  if (idx == 1) {
    prior <- student_t(nu = 20) # alpha
  } else if (idx == 2) {
    prior <- log_normal(mu = 0, sigma = 1) # ell
  } else if (idx == 3) {
    prior <- igam(shape = 2, scale = 1, square = TRUE) # sigma
  } else if (idx == 4) {
    prior <- log_normal(mu = 1, sigma = 1, square = TRUE) # phi
  } else if (idx == 5) {
    prior <- igam(shape = 14, scale = 5) # wrp
  } else if (idx == 6) {
    prior <- bet(a = 0.2, b = 0.2) # beta
  } else if (idx == 7) {
    prior <- bet(a = 1, b = 1) # gamma
  } else {
    prior <- uniform() # effect_time
  }
  return(prior)
}

#' @rdname defaults
default_prior_effect_time_info <- function(num_uncrt) {
  check_not_null(num_uncrt)
  if (num_uncrt > 0) {
    # There is no default prior for effect time
    # TODO: more informative message?
    stop("you must specify 'effect_time_info' in <prior>!")
  } else {
    # Will not be used
    desc <- list(backwards = FALSE, lower = NaN, upper = NaN, zero = NaN)
    desc <- list(desc)
  }
  return(desc)
}

#' @rdname defaults
default_ppc_fun <- function(object) {
  likelihood <- get_obs_model(object)
  f1 <- bayesplot::ppc_dens_overlay
  f2 <- bayesplot::ppc_hist
  fun <- if (likelihood == "gaussian") f1 else f2
  check_type(fun, "function")
  return(fun)
}
