#' Default priors and parameter values
#'
#' @description
#' \itemize{
#'   \item \code{default_vm_params} returns variance mask function parameters
#'   (two numbers)
#'   \item \code{default_prior} returns a named list that defines a prior
#' }
#' @param name parameter name
#' @param num_uncrt number of uncertain components
#' @return see description
#' @name defaults
NULL

#' @rdname defaults
default_vm_params <- function() {
  c(0.025, 1)
}

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
  if (name == "alpha") {
    prior <- student_t(nu = 20)
  } else if (name == "ell") {
    prior <- log_normal(mu = 0, sigma = 1)
  } else if (name == "sigma") {
    prior <- igam(shape = 2, scale = 1, square = TRUE)
  } else if (name == "phi") {
    prior <- log_normal(mu = 1, sigma = 1, square = TRUE)
  } else if (name == "wrp") {
    prior <- igam(shape = 14, scale = 5)
  } else if (name == "beta") {
    prior <- bet(a = 0.2, b = 0.2)
  } else if (name == "gamma") {
    prior <- bet(a = 3, b = 3)
  } else if (name == "effect_time") {
    prior <- uniform()
  } else {
    stop("invalid parameter name '", name, "'")
  }
  return(prior)
}

#' @rdname defaults
default_prior_effect_time_info <- function(num_uncrt) {
  check_not_null(num_uncrt)
  if (num_uncrt > 0) {
    # There is no default prior for effect time
    stop("you must specify 'effect_time_info' in <prior>!")
  } else {
    # Will not be used
    desc <- list(backwards = FALSE, lower = NaN, upper = NaN, zero = NaN)
    desc <- list(desc)
  }
  return(desc)
}
