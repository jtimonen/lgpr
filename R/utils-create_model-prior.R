#' Parse given prior
#'
#' @inheritParams create_model.likelihood
#' @param prior A named list, defining the prior distribution of model
#' (hyper)parameters. See the "Defining priors" section below
#' (\code{\link{lgp}}).
#' @param stan_input a list of stan input fields
#' @return a named list of parsed options
#' @family internal model creation functions
create_model.prior <- function(prior, stan_input, verbose) {
  log_progress("Parsing prior...", verbose)
  num_uncrt <- dollar(stan_input, "num_uncrt")
  num_ns <- dollar(stan_input, "num_ns")
  filled <- fill_prior(prior, num_uncrt)
  spec <- dollar(filled, "specified")
  dflt <- dollar(filled, "defaulted")
  str1 <- paste(spec, collapse = ", ")
  str2 <- paste(dflt, collapse = ", ")
  wrp_defaulted <- "wrp" %in% dflt
  if (num_ns > 0 && wrp_defaulted) {
    model_desc <- "involves a gp_ns() or gp_vm() expression"
    msg <- warn_msg_default_prior("input warping steepness", "wrp", model_desc)
    warning(msg)
  }
  msg1 <- paste0("User-specified priors found for: {", str1, "}.")
  msg2 <- paste0(
    "If any of the following parameters are included in the",
    " model, default priors are used for them: {", str2, "}."
  )
  info <- paste0(msg1, "\n", msg2, "\n")
  log_info(info, verbose)

  raw <- dollar(filled, "prior")
  to_stan <- parse_prior_full(raw, stan_input)
  list(
    to_stan = to_stan,
    raw = raw
  )
}

# Fill a partially defined prior
fill_prior <- function(prior, num_uncrt) {
  par_names <- c("alpha", "ell", "wrp", "sigma", "phi", "beta", "gamma")
  par_names <- c(par_names, "effect_time", "effect_time_info")
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

# Convert the Stan input encoding of a prior to a human-readable data frame
prior_to_df <- function(stan_input, digits = 3) {
  # Positive parameters
  check_positive(digits)
  pnames <- c("alpha", "ell", "wrp", "sigma", "phi")
  df <- NULL
  for (p in pnames) {
    df <- rbind(df, prior_to_df_pos(stan_input, p, digits))
  }

  # Beta
  num_heter <- dollar(stan_input, "num_heter")
  if (num_heter > 0) {
    num_bt <- dollar(stan_input, "num_bt")
    df_bet <- prior_to_df_unit(stan_input, "beta", num_bt, digits)
    df <- rbind(df, df_bet)
  }

  # Gamma
  obs_model <- dollar(stan_input, "obs_model")
  if (obs_model == 5) {
    df_gam <- prior_to_df_unit(stan_input, "gamma", 1, digits)
    df <- rbind(df, df_gam)
  }

  # Effect time
  num_uncrt <- dollar(stan_input, "num_uncrt")
  if (num_uncrt > 0) {
    df_p <- prior_to_df_teff(stan_input, digits)
    df <- rbind(df, df_p)
  }

  return(df)
}

# Helper function for converting prior representation to human readable df
prior_to_df_pos <- function(stan_input, parname, digits) {
  prior <- dollar(stan_input, paste0("prior_", parname))
  hyper <- dollar(stan_input, paste0("hyper_", parname))
  D <- dim(prior)[1]
  pnames <- rep("foo", D)
  dnames <- rep("foo", D)
  bounds <- rep("foo", D)
  for (j in seq_len(D)) {
    par <- paste0(parname, "[", j, "]")
    out <- prior_to_str(par, prior[j, ], hyper[j, ], digits)
    tpar <- dollar(out, "parname")
    pnames[j] <- par
    dnames[j] <- paste0(tpar, " ~ ", dollar(out, "distribution"))
    bounds[j] <- "[0, Inf)"
  }
  df <- data.frame(pnames, bounds, dnames)
  colnames(df) <- c("Parameter", "Bounds", "Prior")
  return(df)
}

# Helper function for converting prior representation to human readable df
prior_to_df_unit <- function(stan_input, parname, num, digits) {
  hyper <- dollar(stan_input, paste0("hyper_", parname))
  check_positive(digits)
  a <- round(hyper[1], digits = digits)
  b <- round(hyper[2], digits = digits)
  dist <- paste0("beta(", a, ", ", b, ")")
  bounds <- "[0, 1]"
  nam1 <- paste0(parname, "[1]")
  nam2 <- paste0(parname, "[1-", num, "]")
  par <- if (num > 1) nam2 else nam1
  dist <- paste0(par, " ~ ", dist)
  df <- data.frame(par, bounds, dist)
  colnames(df) <- c("Parameter", "Bounds", "Prior")
  return(df)
}

# Helper function for converting prior representation to human readable df
prior_to_df_teff <- function(stan_input, digits) {
  num_bt <- dollar(stan_input, "num_bt")
  prior <- dollar(stan_input, "prior_teff")
  type <- prior[1]
  backwards <- prior[2]
  hyper <- dollar(stan_input, "hyper_teff")
  zero <- dollar(stan_input, "teff_zero")
  lower <- dollar(stan_input, "teff_lb")
  upper <- dollar(stan_input, "teff_ub")
  pnames <- rep("foo", num_bt)
  dnames <- rep("foo", num_bt)
  bounds <- rep("foo", num_bt)
  for (j in seq_len(num_bt)) {
    par <- paste0("teff[", j, "]")
    tpar <- par
    tpar <- minus.append(tpar, zero[j])
    tpar <- minus.prepend(tpar, backwards)
    out <- prior_to_str(par, c(type, 0), hyper, digits)
    pnames[j] <- par
    dnames[j] <- paste0(tpar, " ~ ", dollar(out, "distribution"))
    bounds[j] <- paste0("[", lower[j], ", ", upper[j], "]")
  }
  df <- data.frame(pnames, bounds, dnames)
  colnames(df) <- c("Parameter", "Bounds", "Prior")
  return(df)
}

# Append minus and val to a string if val is not zero
minus.append <- function(str, val) {
  if (val != 0) str <- paste0(str, " - ", val)
  return(str)
}

# Prepend minus to a string
minus.prepend <- function(str, prepend) {
  if (prepend != 0) str <- paste0(" - (", str, ")")
  return(str)
}

# Human-readable prior statement
prior_to_str <- function(parname, prior, hyper, digits) {
  hyper <- round(hyper, digits)

  # Check distribution type
  tp <- prior[1]
  check_allowed(tp, seq_len(6))
  names <- prior_type_names()
  pname <- names[tp]

  # Check if there is a transform
  tf <- prior[2]
  check_allowed(tf, c(0, 1))
  if (tf == 1) parname <- paste0("(", parname, ")^2")

  # Get prior statement
  if (tp %in% c(2, 4, 5, 6)) {
    str <- paste0(pname, "(", hyper[1], ",", hyper[2], ")")
  } else if (tp == 3) {
    str <- paste0(pname, "(", hyper[1], ")")
  } else {
    str <- paste(pname, sep = "")
  }

  # Return
  list(parname = parname, distribution = str)
}


# Parse a fully defined prior
parse_prior_full <- function(prior, stan_input) {
  list_pos <- parse_prior_full_pos(prior, stan_input)
  list_unit <- parse_prior_full_unit(prior, stan_input)
  list_teff <- parse_prior_full_teff(prior, stan_input)
  c(list_pos, list_unit, list_teff)
}

# Parse a fully defined prior (positive parameters)
parse_prior_full_pos <- function(prior, stan_input) {
  obs_model <- dollar(stan_input, "obs_model")
  par_names <- c("alpha", "ell", "wrp", "sigma", "phi")
  nums <- stan_input[c("num_comps", "num_ell", "num_ns")]
  num_sigma <- as.numeric(obs_model == 1)
  num_phi <- as.numeric(obs_model == 3)
  nums <- c(unlist(nums), num_sigma, num_phi)

  parsed <- c()
  K <- length(par_names)
  for (k in seq_len(K)) {
    desc <- dollar(prior, par_names[k])
    pp <- parse_prior_single(desc, nums[k])
    f1 <- paste0("prior_", par_names[k])
    f2 <- paste0("hyper_", par_names[k])
    parsed[[f1]] <- dollar(pp, "prior")
    parsed[[f2]] <- dollar(pp, "hyper")
  }
  return(parsed)
}


# Parse a fully defined prior (unit interval parameters)
parse_prior_full_unit <- function(prior, stan_input) {
  obs_model <- dollar(stan_input, "obs_model")
  par_names <- c("beta", "gamma")
  num_gamma <- as.numeric(obs_model == 5)
  num_heter <- dollar(stan_input, "num_heter")
  num_heter <- as.numeric(num_heter > 0)
  nums <- c(num_heter, num_gamma)

  parsed <- c()
  K <- length(par_names)
  for (k in seq_len(K)) {
    desc <- dollar(prior, par_names[k])[[1]]
    hyper <- c(dollar(desc, "alpha"), dollar(desc, "beta"))
    f <- paste0("hyper_", par_names[k])
    parsed[[f]] <- repvec(hyper, nums[k])
  }
  return(parsed)
}

# Parse a fully defined prior (effect time parameters)
parse_prior_full_teff <- function(prior, stan_input) {
  num_bt <- dollar(stan_input, "num_bt")
  num_uncrt <- dollar(stan_input, "num_uncrt")

  effect_time_info <- dollar(prior, "effect_time_info")[[1]]
  is_backwards <- as.numeric(dollar(effect_time_info, "backwards"))
  lower <- dollar(effect_time_info, "lower")
  upper <- dollar(effect_time_info, "upper")
  zero <- dollar(effect_time_info, "zero")
  lower <- ensure_len(lower, num_bt)
  upper <- ensure_len(upper, num_bt)
  zero <- ensure_len(zero, num_bt)
  DIM <- as.numeric(num_uncrt > 0)

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
