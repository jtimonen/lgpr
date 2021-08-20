# Names of allowed prior types
prior_type_names <- function(idx = NULL) {
  names <- c(
    "Uniform",
    "Normal",
    "Student-t",
    "Gamma",
    "Inv-Gamma",
    "Log-Normal"
  )
  names <- tolower(names)
  out <- if (!is.null(idx)) names[idx] else names
  return(out)
}

# Fill a partially defined prior
fill_prior <- function(prior, num_uncrt, par_names) {
  defaulted <- c()
  specified <- c()
  for (name in par_names) {
    if (name %in% names(prior)) {
      # User-specified prior
      specified <- c(specified, name)
      pr <- dollar(prior, name)
      prior[[name]] <- list_if_named(pr)
    } else {
      # Default prior
      defaulted <- c(defaulted, name)
      prior[[name]] <- default_prior(name, num_uncrt)
    }
  }
  list(
    prior = prior,
    defaulted = defaulted,
    specified = specified
  )
}

# Information about which parameters have default priors
defaulting_info <- function(filled, verbose) {
  spec <- dollar(filled, "specified")
  dflt <- dollar(filled, "defaulted")
  str1 <- paste(spec, collapse = ", ")
  str2 <- paste(dflt, collapse = ", ")
  msg1 <- paste0("User-specified priors found for: {", str1, "}.")
  msg2 <- paste0(
    "If any of the following parameters are included in the",
    " model, default priors are used for them: {", str2, "}."
  )
  info <- paste0(msg1, "\n", msg2, "\n")
  log_info(info, verbose)
  return(dflt)
}

# Parse common parameters of a prior
parse_prior_common <- function(prior, si) {
  DIM_BETA <- as.numeric(dollar(si, "num_het") > 0)
  DIM_TEFF <- as.numeric(dollar(si, "num_unc") > 0)
  num_teff <- dollar(si, "num_teff")
  alpha <- parse_prior_pos(prior, dollar(si, "J"), "alpha")
  ell <- parse_prior_pos(prior, dollar(si, "num_ell"), "ell")
  wrp <- parse_prior_pos(prior, dollar(si, "num_wrp"), "wrp")
  beta <- parse_prior_unit(prior, DIM_BETA, "beta")
  teff <- parse_prior_teff(prior, DIM_TEFF, num_teff)
  c(alpha, ell, wrp, beta, teff)
}

# Parse prior of a positive parameter
parse_prior_pos <- function(prior, DIM, par_name) {
  desc <- dollar(prior, par_name) # no [[1]]
  parsed <- list()
  pp <- parse_prior_single(desc, DIM)
  f1 <- paste0("prior_", par_name)
  f2 <- paste0("hyper_", par_name)
  parsed[[f1]] <- dollar(pp, "prior")
  parsed[[f2]] <- dollar(pp, "hyper")
  return(parsed)
}

# Parse a unit interval parameter prior (beta prior)
parse_prior_unit <- function(prior, DIM, par_name) {
  parsed <- list()
  desc <- dollar(prior, par_name)[[1]]
  hyper <- c(dollar(desc, "alpha"), dollar(desc, "beta"))
  f <- paste0("hyper_", par_name)
  parsed[[f]] <- repvec(hyper, DIM)
  return(parsed)
}

# Parse effect time parameter prior
parse_prior_teff <- function(prior, DIM, num_teff) {
  effect_time_info <- dollar(prior, "effect_time_info")[[1]]
  is_backwards <- as.numeric(dollar(effect_time_info, "backwards"))
  lower <- dollar(effect_time_info, "lower")
  upper <- dollar(effect_time_info, "upper")
  zero <- dollar(effect_time_info, "zero")
  lower <- ensure_len(lower, num_teff)
  upper <- ensure_len(upper, num_teff)
  zero <- ensure_len(zero, num_teff)

  desc <- dollar(prior, "effect_time")[[1]]
  out <- prior_to_num(desc)
  type <- dollar(out, "prior")
  hyper <- dollar(out, "hyper")
  prior <- c(type[1], is_backwards)

  # Return
  list(
    prior_teff = repvec(prior, DIM),
    hyper_teff = repvec(hyper, DIM),
    teff_zero = repvec(zero, DIM),
    teff_lb = repvec(lower, DIM),
    teff_ub = repvec(upper, DIM)
  )
}

# Parse given prior for a single parameter type
parse_prior_single <- function(desc, num) {
  check_not_named(desc)
  L <- length(desc)
  err_msg <- paste0("<desc> should have length 1 or ", num, "! found = ", L)
  if (L != num) {
    if (L == 1) {
      # The parameter type has the same prior in all components
      out <- prior_to_num(desc[[1]])
      prior <- repvec(dollar(out, "prior"), num)
      hyper <- repvec(dollar(out, "hyper"), num)
    } else {
      stop(err_msg)
    }
  } else {

    # The parameter type has possibly different prior in different components
    prior <- repvec(c(0, 0), L)
    hyper <- repvec(c(0, 0, 0), L)
    for (j in seq_len(L)) {
      desc_j <- desc[[j]]
      out <- prior_to_num(desc_j)
      pr <- repvec(dollar(out, "prior"), 1)
      hp <- repvec(dollar(out, "hyper"), 1)
      prior[j, ] <- pr
      hyper[j, ] <- hp
    }
  }
  list(
    prior = prior,
    hyper = hyper
  )
}

#' Convert given prior to numeric format
#'
#' @param desc Prior description as a named list, containing fields
#' \itemize{
#'   \item \code{dist} - Distribution name. Must be one of
#'   {'uniform', 'normal', 'student-t', 'gamma', 'inv-gamma', 'log-normal'}
#'   (case-insensitive)
#'   \item \code{square} - Is the prior for a square-transformed parameter.
#' }
#' Other list fields are interpreted as hyperparameters.
#' @return a named list of parsed options
prior_to_num <- function(desc) {
  types <- prior_type_names()
  distribution_name <- dollar(desc, "dist")
  dist_num <- check_allowed(distribution_name, types)
  fields <- names(desc)
  fields <- fields[!(fields %in% c("dist", "square"))]
  hyper <- position_hyper_params(desc[fields])
  sq <- dollar(desc, "square")
  list(
    prior = c(dist_num, as.numeric(sq)),
    hyper = hyper,
    hyper_names = fields
  )
}

# Position the hyper parameters from a list to a vector that goes to Stan
position_hyper_params <- function(desc) {
  hyper <- c(0, 0, 0)
  NAMES <- names(desc)
  H1 <- c("mu", "alpha", "nu")
  H2 <- c("sigma", "beta")
  for (name in NAMES) {
    check_allowed(name, c(H1, H2))
    val <- dollar(desc, name)
    idx <- if (name %in% H2) 2 else 1
    hyper[idx] <- val
  }
  return(hyper)
}

# Default prior, given parameter name and number of uncertain components
default_prior <- function(name, num_uncrt = NULL) {
  if (name == "effect_time_info") {
    desc <- default_prior_effect_time_info(num_uncrt)
  } else {
    desc <- list(default_prior_common(name))
  }
  return(desc)
}

# Default priors for common parameters
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

# Default effect time info
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
