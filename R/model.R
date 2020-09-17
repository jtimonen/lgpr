#' Create a model
#'
#' @inheritParams parse_formula
#' @inheritParams parse_response
#' @inheritParams parse_likelihood
#' @inheritParams parse_prior
#' @inheritParams parse_covs_and_comps
#' @inheritParams parse_options
#' @param prior_only Should this run in prior sampling mode,
#' where likelihood is ignored?
#' @param verbose Should more verbose output be printed?
#' @family main functions
create_model <- function(formula,
                         data,
                         likelihood = "gaussian",
                         prior = NULL,
                         c_hat = NULL,
                         num_trials = NULL,
                         options = NULL,
                         prior_only = FALSE,
                         verbose = FALSE,
                         sample_f = !(likelihood == "gaussian")) {

  # Parse formula and options
  model_formula <- parse_formula(formula, data)
  list_opts <- parse_options(options)

  # Parse response and likelihood
  parsed <- parse_response(data, likelihood, model_formula)
  list_y <- dollar(parsed, "to_stan")
  y_scaling <- dollar(parsed, "scaling")
  list_lh <- parse_likelihood(likelihood, c_hat, num_trials, list_y, sample_f)

  # Parse covariates and components
  parsed <- parse_covs_and_comps(data, model_formula)
  list_x <- dollar(parsed, "to_stan")
  x_cont_scalings <- dollar(parsed, "x_cont_scalings")
  x_cat_levels <- dollar(parsed, "x_cat_levels")
  caseid_map <- dollar(parsed, "caseid_map")

  # Group variable names
  var_names <- list(
    y = model_formula@y_name,
    x_cont = names(x_cont_scalings),
    x_cat = names(x_cat_levels)
  )

  # Parse the prior
  obm <- dollar(list_lh, "obs_model")
  parsed <- parse_prior(prior, list_x, obm)
  full_prior <- dollar(parsed, "raw")
  list_prior <- dollar(parsed, "to_stan")
  if (verbose) {
    cat(dollar(parsed, "info"))
  }

  # Other
  list_other <- list(
    is_verbose = as.numeric(verbose),
    is_likelihood_skipped = as.numeric(prior_only)
  )

  # Create slots of the 'lgpmodel' object
  stan_input <- c(
    list_opts,
    list_y,
    list_lh,
    list_x,
    list_prior,
    list_other
  )
  info <- list(
    created = date(),
    pkg_desc = get_pkg_description(),
    caseid_map = caseid_map
  )
  var_scalings <- list(y = y_scaling, x_cont = x_cont_scalings)
  var_info <- list(x_cat_levels = x_cat_levels)

  # Create the 'lgpmodel' object
  out <- new("lgpmodel",
    model_formula = model_formula,
    var_names = var_names,
    var_info = var_info,
    var_scalings = var_scalings,
    stan_input = stan_input,
    info = info,
    stan_model_name = "lgp",
    full_prior = full_prior
  )
  return(out)
}

#' Parse the given modeling options
#'
#' @param options A named list with the following possible fields:
#' \itemize{
#'   \item \code{delta} Amount of added jitter to ensure positive definite
#'   covariance matrices.
#' }
#' @return a named list of parsed options
parse_options <- function(options = NULL) {
  input <- options

  # Set defaults
  opts <- list(
    skip_generated = FALSE,
    delta = 1e-8
  )

  # Replace defaults if found from input
  for (opt_name in names(opts)) {
    if (opt_name %in% names(input)) {
      opts[[opt_name]] <- input[[opt_name]]
    }
  }

  # Format for Stan input
  list(
    is_generated_skipped = as.numeric(dollar(opts, "skip_generated")),
    delta = dollar(opts, "delta")
  )
}

#' Create the function that does a standardizing transform and its inverse
#'
#' @param y variable measurements
#' @param name variable name
#' @return an object of class \linkS4class{lgpscaling}
create_scaling <- function(y, name) {
  if (length(y) < 2) {
    stop("length of <y> must be at least 2!")
  }
  m <- mean(y)
  std <- stats::sd(y)
  if (std == 0) {
    stop("the varible measurements have zero variance!")
  }
  fun <- function(x) {
    (x - m) / std
  }
  fun_inv <- function(x) {
    x * std + m
  }
  new("lgpscaling", fun = fun, fun_inv = fun_inv, var_name = name)
}


#' Creating the idx_expand input for Stan
#'
#' @name idx_expand
#' @param components the \code{components} input for Stan
#' @param x_cat the \code{x_cat} input for Stan
#' @param x_cont_mask the \code{x_cont_mask} input for Stan
#' @param x_fac object returned by \code{\link{create_idx_expand_picker}}
#' @param map object returned by \code{\link{map_factor_to_caseid}}
#' @return the \code{idx_expand} input for Stan and a mapping data frame
NULL

#' @rdname idx_expand
create_idx_expand <- function(components, x_cat, x_cont_mask) {
  pick <- create_idx_expand_picker(components, x_cat)
  x_fac <- dollar(pick, "x_fac")
  factor_name <- dollar(pick, "factor_name")
  inds <- which(components[, 4] + components[, 7] > 0)
  inds <- as.numeric(inds) # remove names
  L <- length(inds)
  if (L == 0) {
    idx_expand <- x_fac
    map <- NULL
  } else {
    i_cont <- components[inds, 9]
    to_red <- x_cont_mask[i_cont, ]
    if (length(i_cont) == 1) {
      to_red <- repvec(to_red, 1)
    }
    idx_mask <- reduce_rows(to_red)
    map <- map_factor_to_caseid(x_fac, idx_mask, factor_name)
    idx_expand <- map_caseid_to_row(x_fac, map) + 1 # note the plus 1
  }

  # Return
  list(map = map, idx_expand = idx_expand)
}

#' @rdname idx_expand
create_idx_expand_picker <- function(components, x_cat) {
  n_obs <- dim(x_cat)[2]
  inds <- c(components[, 4], components[, 7])
  inds <- as.numeric(inds[inds != 0])
  J <- length(inds)
  if (J == 0) {
    x_fac <- rep(1, n_obs)
    name <- NULL
  } else {
    all_same <- all(inds == inds[1])
    if (!all_same) {
      str <- paste(inds, collapse = ", ")
      msg <- paste0(
        "The het() and unc() expressions must have the same ",
        "categorical covariate in every term! ",
        "Found inds = {", str, "}"
      )
      stop(msg)
    }
    i1 <- inds[1]
    x_fac <- as.numeric(x_cat[i1, ])
    name <- rownames(x_cat)[i1]
  }

  # Return
  list(x_fac = x_fac, factor_name = name)
}

#' @rdname idx_expand
#' @param factor_name factor name
map_factor_to_caseid <- function(x_fac, x_cont_mask, factor_name) {
  id <- c()
  case_id <- c()
  num_cases <- 0
  facs <- unique(x_fac)
  for (u in facs) {
    inds <- which(x_fac == u)
    vals <- x_cont_mask[inds]
    L <- length(vals)
    all0 <- all.equal(vals, rep(0, L)) == TRUE
    all1 <- all.equal(vals, rep(1, L)) == TRUE
    if (all0) {
      num_cases <- num_cases + 1
      id[num_cases] <- u
      case_id[num_cases] <- num_cases
    } else if (all1) {
      num_cases <- num_cases + 1
    } else {
      stop(paste0("inconsistent x_cont_mask values for x_fac = ", u))
    }
  }
  out <- data.frame(id, case_id)
  colnames(out) <- c(factor_name, "case_id")
  return(out)
}

#' @rdname idx_expand
map_caseid_to_row <- function(x_fac, map) {
  out <- rep(0, length(x_fac))
  num_levels <- dim(map)[1]
  for (j in seq_len(num_levels)) {
    id <- map[j, 1]
    caseid <- map[j, 2]
    inds <- which(x_fac == id)
    out[inds] <- caseid
  }
  return(out)
}

#' Parse the response variable from given data and formula
#'
#' @param data A \code{data.frame} where each column corresponds to one
#' variable, and each row is one observation. Continuous covariates and the
#' response variable must have type \code{"numeric"} and categorical covariates
#' must have type \code{"factor"}. Missing values should be indicated with
#' \code{NaN} or \code{NA}. The response variable cannot contain missing
#' values.
#' @param likelihood Determines the observation model. Must be either
#' \code{"gaussian"} (default), \code{"poisson"}, \code{"nb"} (negative
#' binomial), \code{"binomial"} or \code{"bb"} (beta binomial).
#' @param model_formula An object of class \code{lgpformula}.
#' @return a named list of parsed options
parse_response <- function(data, likelihood, model_formula) {
  obs_model <- likelihood_as_int(likelihood)

  # Check that data is a data.frame and contains the response variable
  y_name <- model_formula@y_name
  c_data <- class(data)
  if (c_data != "data.frame") {
    stop("<data> must be a data.frame! found = ", c_data)
  }
  check_in_data(y_name, data)
  Y_RAW <- data[[y_name]]

  # Check that the response is numeric and compatible with observation model
  check_response(Y_RAW, obs_model)

  # Create y_cont and y_disc inputs for Stan
  num_obs <- length(Y_RAW)
  if (obs_model != 1) {
    y_disc <- array(Y_RAW, dim = c(1, num_obs))
    y_cont <- array(Y_RAW, dim = c(0, num_obs))
    normalizer <- new("lgpscaling", var_name = y_name)
  } else {
    normalizer <- create_scaling(Y_RAW, y_name) # create scaling and inverse
    y_disc <- array(Y_RAW, dim = c(0, num_obs))
    Y_NORM <- call_fun(normalizer@fun, Y_RAW) # standardize the response
    y_cont <- array(Y_NORM, dim = c(1, num_obs))
  }
  to_stan <- list(
    num_obs = num_obs,
    y_disc = y_disc,
    y_cont = y_cont
  )

  # Return also the scaling and its inverse
  list(
    to_stan = to_stan,
    scaling = normalizer
  )
}


#' Parse the given observation model
#'
#' @inheritParams parse_response
#' @param c_hat The GP mean. Must be a vector of length \code{dim(data)[1]}, or
#' a real number defining a constant GP mean. If \code{NULL}, this is set to
#'  \itemize{
#'    \item \code{c_hat = 0}, if \code{likelihood} is \code{"gaussian"}, because
#'    with Gaussian likelihood the response variable is by default centered to
#'    have zero mean.
#'    \item \code{c_hat = } \code{log(mean(y))} if \code{likelihood} is
#'    \code{"poisson"} or \code{"nb"},
#'    \item \code{c_hat = } \code{log(p/(1-p))}, where
#'    \code{p = mean(y/num_trials)} if \code{likelihood} is \code{"binomial"},
#'  }
#' where \code{y} denotes the response variable. You can modify this vector to
#' account for normalization between data points. With Gaussian likelihood
#' though, do not modify this argument, normalize the data beforehand instead.
#' @param num_trials This argument (number of trials) is only needed when
#' likelihood is binomial or beta binomial. Must have length one or equal
#' to number of data points. Setting \code{num_trials=1} and
#' \code{likelihood="binomial"} corresponds to Bernoulli observation model.
#' @param list_y a list field returned by \code{\link{parse_response}}
#' @param sample_f Determines if the latent function values are be sampled
#' (must be \code{TRUE} if likelihood is not \code{"gaussian"}).
#' @return a list of parsed options
parse_likelihood <- function(likelihood, c_hat, num_trials, list_y, sample_f) {
  LH <- likelihood_as_int(likelihood)
  y <- if (LH != 1) dollar(list_y, "y_disc") else dollar(list_y, "y_cont")
  y <- as.numeric(y)
  num_obs <- length(y)
  num_trials <- set_num_trials(num_trials, num_obs, LH)
  c_hat <- set_c_hat(c_hat, y, LH, num_trials)
  sample_f_final <- if (LH != 1) TRUE else sample_f
  if (sample_f_final != sample_f) {
    stop("sample_f must be TRUE when likelihood is ", likelihood)
  }
  list(
    obs_model = LH,
    y_num_trials = num_trials,
    c_hat = c_hat,
    is_f_sampled = as.numeric(sample_f_final)
  )
}

#' Set c_hat (non-gaussian observation models)
#'
#' @param c_hat the \code{c_hat} argument given as input to
#' \code{\link{create_model}}
#' @param response response variable data
#' @param LH likelihood as int
#' @param num_trials the num_trials data (binomial likelihood)
#' @return a real number
set_c_hat <- function(c_hat, response, LH, num_trials) {
  nb_or_pois <- LH %in% c(2, 3)
  binomial <- LH %in% c(4, 5)
  if (is.null(c_hat)) {
    if (nb_or_pois) {
      c_hat <- log(mean(response))
    } else if (binomial) {
      check_all_leq(response, num_trials)
      p <- mean(response / num_trials)
      c_hat <- log(p / (1 - p))
    } else {
      c_hat <- 0
    }
  } else {
    if (LH == 1) {
      stop(
        "Only give the c_hat argument if observation model is not Gaussian!",
        " With Gaussian likelihood, you should use",
        " c_hat = NULL, in which case the GP mean will be set to zero",
        " and the response variable is standardized to have mean zero."
      )
    }
  }
  n <- length(response)
  L <- length(c_hat)
  if (L != 1 && L != n) {
    stop(
      "Invalid length of <c_hat>! Must be 1 or equal to number ",
      "of observartions (", n, "). Found = ", L
    )
  }
  if (L == 1) {
    c_hat <- rep(c_hat, n)
  }
  return(c_hat)
}

#' Set num_trials (binomial and beta-binomial observation models)
#'
#' @param num_trials the \code{num_trials} argument
#' @param num_obs number of observations
#' @param LH likelihood as integer encoding
#' @return a numeric vector
set_num_trials <- function(num_trials, num_obs, LH) {
  is_binom <- LH %in% c(4, 5)
  if (is.null(num_trials)) {
    num_trials <- rep(1, num_obs)
  } else {
    if (!is_binom) {
      msg <- paste0(
        "Only give the <num_trials> argument if likelihood is",
        " binomial or beta-binomial (bb)!"
      )
      stop(msg)
    }
    L <- length(num_trials)
    if (L == 1) {
      num_trials <- rep(num_trials, num_obs)
    } else if (L != num_obs) {
      stop(
        "Invalid length of <num_trials>! Must be 1 or equal to number ",
        "of observartions (", num_obs, "). Found = ", L
      )
    }
  }
  DIM <- as.numeric(is_binom)
  num_trials <- array(num_trials, dim = c(DIM, num_obs))
  return(num_trials)
}


#' Check that the response is numeric and compatibile with observation model
#'
#' @param y the response variable measurements
#' @param obs_model observation model (integer encoding)
check_response <- function(y, obs_model) {

  # Check that y is numeric
  c_y <- class(y)
  if (c_y != "numeric") {
    stop("the response variable must be numeric! found = ", c_y)
  }

  # Check compatibility with observation model
  if (obs_model != 1) {
    diff <- max(abs(y - round(y)))
    if (diff > 0) {
      stop(
        "the response variable should contain only integers ",
        "with this observation model"
      )
    }
    if (any(y < 0)) {
      stop(
        "the response variable measurements cannot be negative ",
        "with this observation model!"
      )
    }
  }
  return(TRUE)
}
