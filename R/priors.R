#' Prior definitions
#'
#' @param square is prior for a square-transformed parameter?
#' @name priors
#' @aliases normal, log_normal, gam, igam, uniform, student_t, bet
#' @return a named list
#' @description These use the same parametrizations as defined in the Stan
#' documentation. See the docs for
#' \href{https://mc-stan.org/docs/2_26/functions-reference/gamma-distribution.html}{gamma}
#' and
#' \href{https://mc-stan.org/docs/2_26/functions-reference/inverse-gamma-distribution.html}{inverse gamma}
#' distributions.
#' @examples
#' # Log-normal prior
#' log_normal(mu = 1, sigma = 1)
#'
#' # Cauchy prior
#' student_t(nu = 1)
#'
#' # Exponential prior with rate = 0.1
#' gam(shape = 1, inv_scale = 0.1)
#'
#' # Create a similar priors as in LonGP (Cheng et al., 2019)
#' # Not recommended, because a lengthscale close to 0 is possible.
#' a <- log(1) - log(0.1)
#' log_normal(mu = 0, sigma = a / 2) # for continuous lengthscale
#' student_t(nu = 4) # for interaction lengthscale
#' igam(shape = 0.5, scale = 0.005, square = TRUE) # for sigma
NULL

#' @export
#' @rdname priors
uniform <- function(square = FALSE) {
  list(
    dist = "uniform",
    square = square
  )
}

#' @export
#' @rdname priors
#' @param mu mean
#' @param sigma standard deviation
normal <- function(mu, sigma, square = FALSE) {
  check_numeric(mu)
  check_positive(sigma)
  list(
    dist = "normal",
    square = square,
    mu = mu,
    sigma = sigma
  )
}

#' @export
#' @rdname priors
#' @param nu degrees of freedom
student_t <- function(nu, square = FALSE) {
  check_positive(nu)
  list(
    dist = "student-t",
    square = square,
    nu = nu
  )
}

#' @export
#' @rdname priors
#' @param shape shape parameter (alpha)
#' @param inv_scale inverse scale parameter (beta)
gam <- function(shape, inv_scale, square = FALSE) {
  check_positive(shape)
  check_positive(inv_scale)
  list(
    dist = "gamma",
    alpha = shape,
    beta = inv_scale,
    square = square
  )
}

#' @export
#' @rdname priors
#' @param shape shape parameter (alpha)
#' @param scale scale parameter (beta)
#' @family functions related to the inverse-gamma distribution
igam <- function(shape, scale, square = FALSE) {
  check_positive(shape)
  check_positive(scale)
  list(
    dist = "inv-gamma",
    alpha = shape,
    beta = scale,
    square = square
  )
}

#' @export
#' @rdname priors
#' @param mu mean
#' @param sigma standard deviation
log_normal <- function(mu, sigma, square = FALSE) {
  check_numeric(mu)
  check_positive(sigma)
  list(
    dist = "log-normal",
    square = square,
    mu = mu,
    sigma = sigma
  )
}


#' @export
#' @rdname priors
#' @param a shape parameter
#' @param b shape parameter
bet <- function(a, b) {
  check_positive(a)
  check_positive(b)
  list(
    dist = "beta",
    square = FALSE,
    alpha = a,
    beta = b
  )
}


#' Parse given prior
#'
#' @param prior A named list, defining the prior distribution of model
#' (hyper)parameters. See the "Defining priors" section below
#' (\code{\link{lgp}}).
#' @param stan_input a list of stan input fields
#' @param obs_model observation model as integer
#' @return a named list of parsed options
parse_prior <- function(prior, stan_input, obs_model) {
  num_uncrt <- dollar(stan_input, "num_uncrt")
  num_ns <- get_num_ns(stan_input)
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
  msg1 <- paste0("  * user-specified priors found for: {", str1, "}")
  msg2 <- paste0(
    "  * if any of the following parameters are included in the",
    " model, default priors are used for them: {", str2, "}"
  )
  info <- paste0(msg1, "\n", msg2, "\n")
  raw <- dollar(filled, "prior")
  to_stan <- parse_prior_full(raw, stan_input, obs_model)
  list(
    to_stan = to_stan,
    info = info,
    raw = raw
  )
}

#' Fill a partially defined prior
#'
#' @inheritParams parse_prior
#' @param num_uncrt number of uncertain components
#' @return a named list
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

#' Names of allowed  prior types
#'
#' @param idx an integer or NULL
#' @return a character vector
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

#' Convert the Stan input encoding of a prior to a human-readable format
#'
#' @param stan_input a list of stan input fields
#' @inheritParams prior_to_str
#' @name prior_df
#' @return a data frame
NULL

#' @rdname prior_df
prior_to_df <- function(stan_input, digits = 3) {

  # Positive parameters
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

#' @rdname prior_df
#' @param parname parameter name
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

#' @rdname prior_df
#' @param parname parameter name
prior_to_df_unit <- function(stan_input, parname, digits) {
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

#' @rdname prior_df
#' @param num number of parameters of the type
prior_to_df_unit <- function(stan_input, parname, num, digits) {
  hyper <- dollar(stan_input, paste0("hyper_", parname))
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

#' @rdname prior_df
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

#' Add minus to a string depending on options
#'
#' @param str a string
#' @param val a value to append with minus (will not be appended if value
#' is zero)
#' @param prepend should a minus be prepended (true if this is not zero)
#' @name minus

#' @rdname minus
minus.append <- function(str, val) {
  if (val != 0) str <- paste0(str, " - ", val)
  return(str)
}

#' @rdname minus
minus.prepend <- function(str, prepend) {
  if (prepend != 0) str <- paste0(" - (", str, ")")
  return(str)
}

#' Human-readable prior statement
#'
#' @param parname parameter name
#' @param prior two integers
#' @param hyper three real numbers
#' @param digits number of digits to show for floats
#' @return A list.
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
