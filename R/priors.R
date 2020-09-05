#' Prior definitions
#'
#' @param square is prior for a square-transformed parameter?
#' @name priors
#' @aliases normal, log_normal, gam, igam, uniform, student_t, bet
#' @return a named list
#' @description These use the same parametrizations as defined in the Stan
#' documentation. See the docs for
#' \href{https://mc-stan.org/docs/2_24/functions-reference/gamma-distribution.html}{gamma} and
#' \href{https://mc-stan.org/docs/2_24/functions-reference/inverse-gamma-distribution.html}{inverse gamma} distributions.
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

#' Default prior for a named parameter
#'
#' @param name parameter name
#' @return a named list
prior_get_default <- function(name) {
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
  } else if (name == "effect_time") {
    prior <- uniform()
  } else {
    stop("invalid parameter name '", name, "'")
  }
  return(prior)
}

#' Parse given prior
#'
#' @param prior A named list, defining the prior distribution of model
#' (hyper)parameters.
#' @param stan_input a list of stan input fields
#' @param obs_model observation model as integer
#' @return a named list of parsed options
parse_prior <- function(prior, stan_input, obs_model) {
  num_uncrt <- dollar(stan_input, "num_uncrt")
  filled <- fill_prior(prior, num_uncrt)
  spec <- dollar(filled, "specified")
  def <- dollar(filled, "defaulted")
  str1 <- paste(spec, collapse = ", ")
  str2 <- paste(def, collapse = ", ")
  msg1 <- paste0("User-specified priors found for: {", str1, "}")
  msg2 <- paste0("Default priors used for: {", str2, "}")
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
  names <- c("alpha", "ell", "wrp", "sigma", "phi", "beta")
  names <- c(names, "effect_time", "effect_time_info")
  defaulted <- c()
  specified <- c()
  for (name in names) {
    if (name %in% names(prior)) {
      # User-specified prior
      specified <- c(specified, name)
      pr <- prior[[name]]
      is_named <- !is.null(names(pr))
      if (is_named) {
        pr <- list(pr)
      }
      prior[[name]] <- pr
    } else {
      # Default prior
      defaulted <- c(defaulted, name)
      if (name == "effect_time_info") {
        if (num_uncrt > 0) {
          msg <- "you must specify 'effect_time_info' in the prior list!"
          stop(msg)
        } else {
          desc <- list(backwards = FALSE, lower = NaN, upper = NaN, zero = NaN)
          desc <- list(desc)
        }
      } else {
        desc <- list(prior_get_default(name))
      }
      prior[[name]] <- desc
    }
  }
  list(
    prior = prior,
    defaulted = defaulted,
    specified = specified
  )
}

#' Parse given fully defined prior
#'
#' @inheritParams parse_prior
#' @return a named list of parsed options
parse_prior_full <- function(prior, stan_input, obs_model) {

  # Count number of different parameters
  nums <- stan_input[c("num_comps", "num_ell", "num_ns")]
  num_sigma <- as.numeric(obs_model == 1)
  num_phi <- as.numeric(obs_model == 3)
  num_uncrt <- dollar(stan_input, "num_uncrt")
  num_heter <- dollar(stan_input, "num_heter")
  num_bt <- dollar(stan_input, "num_bt")
  pnames <- c("alpha", "ell", "wrp", "sigma", "phi")
  nums <- c(unlist(nums), num_sigma, num_phi)
  names(nums)[4:5] <- c("num_sigma", "num_phi")

  # Common parameters
  common <- c()
  K <- length(pnames)
  for (k in seq_len(K)) {

    # Ensure that prior definition exists for this parameter
    desc <- dollar(prior, pnames[k])

    # Prior type and hyper params to stan input format
    pp <- parse_prior_single(desc, nums[k])
    f1 <- paste0("prior_", pnames[k])
    f2 <- paste0("hyper_", pnames[k])
    common[[f1]] <- dollar(pp, "prior")
    common[[f2]] <- dollar(pp, "hyper")
  }

  # Beta and teff parameters
  common[["hyper_beta"]] <- create_hyper_beta(prior, num_heter)
  desc <- dollar(prior, "effect_time")
  info <- dollar(prior, "effect_time_info")
  teff <- parse_prior_teff(desc, info, num_uncrt, num_bt)

  # Return
  c(common, teff)
}

#' Create the teff_ inputs to Stan
#'
#' @inheritParams prior_to_num
#' @param info effect time info
#' @param num_uncrt the \code{num_uncrt} input created for Stan
#' @param num_bt the \code{num_bt} input created for Stan
#' @return a named list
parse_prior_teff <- function(desc, info, num_uncrt, num_bt) {
  effect_time_info <- info[[1]]
  is_backwards <- as.numeric(dollar(effect_time_info, "backwards"))
  lower <- dollar(effect_time_info, "lower")
  upper <- dollar(effect_time_info, "upper")
  zero <- dollar(effect_time_info, "zero")
  lower <- ensure_len(lower, num_bt)
  upper <- ensure_len(upper, num_bt)
  zero <- ensure_len(zero, num_bt)
  DIM <- as.numeric(num_uncrt > 0)

  out <- prior_to_num(desc[[1]])
  type <- dollar(out, "prior")
  hyper <- dollar(out, "hyper")
  prior <- c(type[1], is_backwards)

  # Return
  out <- list(
    prior_teff = repvec(prior, DIM),
    hyper_teff = repvec(hyper, DIM),
    teff_zero = repvec(zero, DIM),
    teff_lb = repvec(lower, DIM),
    teff_ub = repvec(upper, DIM)
  )
}

#' Create the hyper_beta input to Stan
#'
#' @inheritParams parse_prior_full
#' @param num_heter the \code{num_heter} input created for Stan
#' @return an array of size (0, 2) or (1, 2)
create_hyper_beta <- function(prior, num_heter) {
  desc <- dollar(prior, "beta")
  L <- length(desc)
  if (L != 1) {
    stop("length of prior$beta must be 1!")
  }
  a <- dollar(desc[[1]], "alpha")
  b <- dollar(desc[[1]], "beta")
  bet <- c(a, b)
  n_rows <- as.numeric(num_heter > 0)
  out <- repvec(bet, n_rows)
  return(out)
}

#' Parse the given prior for a single parameter type
#'
#' @param desc An unnamed list of prior descriptions.
#' @param num number of parameters of this type
#' @return a named list of parsed options
parse_prior_single <- function(desc, num) {
  nams <- names(desc)
  if (!is.null(nams)) {
    str <- paste(nams, collapse = ", ")
    stop("the list <desc> should not have names! found = {", str, "}")
  }
  L <- length(desc)
  if (L != num) {
    if (L == 1) {

      # The parameter type has the same prior in all components
      out <- prior_to_num(desc[[1]])
      prior <- repvec(dollar(out, "prior"), num)
      hyper <- repvec(dollar(out, "hyper"), num)
    } else {
      msg <- paste0(
        "<desc> should be a list of length 1 or ", num,
        "! found = ", L
      )
      stop(msg)
    }
  } else {

    # The parameter type has possibly different prior in different components
    prior <- repvec(c(0, 0), L)
    hyper <- repvec(c(0, 0, 0), L)
    for (j in seq_len(L)) {
      out <- prior_to_num(desc[[j]])
      prior[j, ] <- repvec(dollar(out, "prior"), num)
      hyper[j, ] <- repvec(dollar(out, "hyper"), num)
    }
  }
  list(
    prior = prior,
    hyper = hyper
  )
}

#' Parse the given prior to numeric format
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

#' Position the hyper parameters from a list to a vector that goes to Stan
#'
#' @param desc Hyperparameters as a named list
#' @return three real numbers
position_hyper_params <- function(desc) {
  hyper <- c(0, 0, 0)
  NAMES <- names(desc)
  for (name in NAMES) {
    val <- desc[[name]]
    H1 <- c("mu", "alpha", "nu")
    H2 <- c("sigma", "beta")
    if (name %in% H1) {
      hyper[1] <- val
    } else if (name %in% H2) {
      hyper[2] <- val
    } else {
      msg <- paste0(
        "invalid hyperparameter name '", name, "'! must be one of {",
        paste(c(H1, H2), collapse = ", "), "}"
      )
      stop(msg)
    }
  }
  return(hyper)
}


#' Names of allowed  prior types
#'
#' @param idx an integer or NULL
#' @return a character vector
prior_type_names <- function(idx = NULL) {
  names <- c(
    "Uniform", "Normal", "Student-t",
    "Gamma", "Inv-Gamma", "Log-Normal"
  )
  names <- tolower(names)
  if (!is.null(idx)) {
    return(names[idx])
  } else {
    return(names)
  }
}

#' Convert the Stan input encoding of a prior to a human-readable format
#'
#' @param stan_input a list of stan input fields
#' @inheritParams prior_to_char
#' @name prior_df
#' @return a data frame
NULL

#' @rdname prior_df
prior_to_df <- function(stan_input, digits = 3) {

  # Common parameters
  pnames <- c("alpha", "ell", "wrp", "sigma", "phi")
  df <- NULL
  for (p in pnames) {
    df_p <- prior_to_df_oneparam(stan_input, p, digits)
    df <- rbind(df, df_p)
  }

  # Beta
  num_heter <- dollar(stan_input, "num_heter")
  num_uncrt <- dollar(stan_input, "num_uncrt")
  if (num_heter > 0) {
    df_p <- prior_to_df_beta(stan_input)
    df <- rbind(df, df_p)
  }

  # Effect time
  if (num_uncrt > 0) {
    df_p <- prior_to_df_teff(stan_input, digits)
    df <- rbind(df, df_p)
  }

  return(df)
}

#' @rdname prior_df
#' @param parname parameter name
prior_to_df_oneparam <- function(stan_input, parname, digits) {
  prior <- dollar(stan_input, paste0("prior_", parname))
  hyper <- dollar(stan_input, paste0("hyper_", parname))
  D <- dim(prior)[1]
  pnames <- rep("foo", D)
  dnames <- rep("foo", D)
  bounds <- rep("foo", D)
  for (j in seq_len(D)) {
    par <- paste0(parname, "[", j, "]")
    out <- prior_to_char(par, prior[j, ], hyper[j, ], digits)
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
prior_to_df_beta <- function(stan_input) {
  num_bt <- dollar(stan_input, "num_bt")
  hyper <- dollar(stan_input, paste0("hyper_beta"))
  dist <- paste0("beta(", hyper[1], ", ", hyper[2], ")")
  bounds <- "[0, 1]"
  par <- paste0("beta[1-", num_bt, "]")
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
    if (zero[j] != 0) {
      tpar <- paste0(tpar, " - ", zero[j])
    }
    if (backwards == 1) {
      tpar <- paste0(" - (", tpar, ")")
    }
    out <- prior_to_char(par, c(type, 0), hyper, digits)
    pnames[j] <- par
    dnames[j] <- paste0(tpar, " ~ ", dollar(out, "distribution"))
    bounds[j] <- paste0("[", lower[j], ", ", upper[j], "]")
  }
  df <- data.frame(pnames, bounds, dnames)
  colnames(df) <- c("Parameter", "Bounds", "Prior")
  return(df)
}


#' Human-readable prior statement
#'
#' @param parname parameter name
#' @param prior two integers
#' @param hyper three real numbers
#' @param digits number of digits to show for floats
#' @return A list.
prior_to_char <- function(parname, prior, hyper, digits) {
  hyper <- round(hyper, digits)

  # Check distribution type
  tp <- prior[1]
  if (tp < 1 || tp > 6) {
    stop("Prior type must 1, 2, 3, 4, 5 or 6! Found = ", tp)
  }
  names <- prior_type_names()
  pname <- names[tp]

  # Check if there is a transform
  tf <- prior[2]
  if (tf == 1) {
    parname <- paste0("(", parname, ")^2")
  } else if (tf == 0) {
    parname <- parname
  } else {
    msg <- paste0(
      "Transform must be either 0 (identity) or 1 (squaring). ",
      "Found = ", tf
    )
    stop(msg)
  }

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
