#' Parse the response variable and its likelihood model
#'
#' @param data A \code{data.frame} where each column corresponds to one
#' variable, and each row is one observation. Continuous covariates and the
#' response variable must have type \code{"numeric"} and categorical covariates
#' must have type \code{"factor"}. Missing values should be indicated with
#' \code{NaN} or \code{NA}. The response variable cannot contain missing
#' values. Column names should not contain trailing or leading underscores.
#' @param likelihood Determines the observation model. Must be either
#' \code{"gaussian"} (default), \code{"poisson"}, \code{"nb"} (negative
#' binomial), \code{"binomial"} or \code{"bb"} (beta binomial).
#' @param sample_f Determines if the latent function values are sampled
#' (must be \code{TRUE} if likelihood is not \code{"gaussian"}). If this is
#' \code{TRUE}, the response variable will be normalized to have zero mean
#' and unit variance.
#' @param c_hat The GP mean. This should only be given if \code{sample_f} is
#' \code{TRUE}, otherwise the GP will always have zero mean. If \code{sample_f}
#' is \code{TRUE}, the given \code{c_hat} can be a vector of length
#' \code{dim(data)[1]}, or a real number defining a constant GP mean. If not
#' specified and \code{sample_f} is \code{TRUE}, \code{c_hat} is set to
#'  \itemize{
#'    \item \code{c_hat = mean(y)}, if \code{likelihood} is \code{"gaussian"},
#'    \item \code{c_hat = } \code{log(mean(y))} if \code{likelihood} is
#'    \code{"poisson"} or \code{"nb"},
#'    \item \code{c_hat = } \code{log(p/(1-p))}, where
#'    \code{p = mean(y/num_trials)} if \code{likelihood} is \code{"binomial"}
#'    or \code{"bb"},
#'  }
#' where \code{y} denotes the response variable measurements.
#' @param num_trials This argument (number of trials) is only needed when
#' likelihood is \code{"binomial"} or \code{"bb"}. Must have length one or
#' equal to the number of data points. Setting \code{num_trials=1} and
#' \code{likelihood="binomial"} corresponds to Bernoulli observation model.
#' @param y_name Name of response variable
#' @param verbose Should some informative messages be printed?
#' @return a list of parsed options
#' @name parse_y
NULL

#' @rdname parse_y
parse_y <- function(data, likelihood, c_hat, num_trials,
                    y_name, sample_f, verbose) {
  if (verbose) cat("Parsing response and likelihood...\n")
  LH <- likelihood_as_int(likelihood)

  # Check that data contains the response variable,
  # which is numeric and compatible with observation model
  check_in_data(y_name, data)
  Y_RAW <- dollar(data, y_name)
  check_response(Y_RAW, LH)

  # Call different subroutines depending on sample_f
  if (!sample_f) {
    if (LH != 1) {
      stop("<sample_f> must be TRUE when <likelihood> is ", likelihood)
    }
    if (!is.null(c_hat)) {
      stop("<c_hat> must be NULL when <sample_f> is FALSE!")
    }
    if (!is.null(num_trials)) {
      stop("<num_trials> must be NULL when <sample_f> is FALSE!")
    }
    y_info <- parse_y.marginal(Y_RAW, y_name)
  } else {
    y_info <- parse_y.latent(Y_RAW, y_name, LH, c_hat, num_trials)
  }
  return(y_info)
}

#' @rdname parse_y
parse_y.marginal <- function(Y_RAW, y_name) {
  num_obs <- length(Y_RAW)
  normalizer <- create_scaling(Y_RAW, y_name) # create scaling and inverse
  y_norm <- call_fun(normalizer@fun, Y_RAW) # standardize the response

  # Return Stan inputs and also the scaling and its inverse
  list(
    to_stan = list(
      num_obs = num_obs,
      y_norm = y_norm,
      obs_model = 1
    ),
    scaling = normalizer
  )
}

#' @rdname parse_y
parse_y.latent <- function(Y_RAW, y_name, LH, c_hat, num_trials) {
  num_obs <- length(Y_RAW)
  if (LH == 1) {
    y_int <- array(Y_RAW, dim = c(0, num_obs))
    y_real <- array(Y_RAW, dim = c(1, num_obs))
  } else {
    y_int <- array(Y_RAW, dim = c(1, num_obs))
    y_real <- array(Y_RAW, dim = c(0, num_obs))
  }
  num_trials <- set_num_trials(num_trials, Y_RAW, LH)
  c_hat <- set_c_hat(c_hat, Y_RAW, LH, num_trials)

  # Create stan input parts
  to_stan <- list(
    num_obs = num_obs,
    y_int = y_int,
    y_real = y_real,
    obs_model = LH,
    y_num_trials = num_trials,
    c_hat = c_hat
  )

  # Return stan input and dummy normalizer, which is identity mapping
  list(to_stan = to_stan, scaling = new("lgpscaling", var_name = y_name))
}

#' Set c_hat
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
  gaussian <- LH == 1

  # Create  c_hat if it is NULL
  if (is.null(c_hat)) {
    if (nb_or_pois) {
      c_hat <- log(mean(response)) # Poisson or  NB
    } else if (binomial) {
      check_all_leq(response, num_trials)
      p <- mean(response / num_trials)
      c_hat <- log(p / (1 - p)) # Binomial or BB
    } else if (gaussian) {
      c_hat <- mean(response) # Gaussian
    } else {
      stop("unknown likelihood LH =", LH)
    }
  }
  n <- length(response)
  L <- length(c_hat)

  # Throw error if given c_hat is misspecified
  if (L != 1 && L != n) {
    stop(
      "Invalid length of <c_hat>! Must be 1 or equal to number ",
      "of observartions (", n, "). Found = ", L
    )
  }
  if (L == 1) c_hat <- rep(c_hat, n) # given c_hat is one number
  return(c_hat)
}

#' Set num_trials (binomial and beta-binomial observation models)
#'
#' @param num_trials the \code{num_trials} argument
#' @param y response variable observations
#' @param LH likelihood as integer encoding
#' @return a numeric vector
set_num_trials <- function(num_trials, y, LH) {
  num_obs <- length(y)
  is_binom <- LH %in% c(4, 5)
  if (is.null(num_trials)) {
    num_trials <- rep(1, num_obs)
  } else {
    if (!is_binom) {
      msg <- paste0(
        "Only give the <num_trials> argument if likelihood is",
        " 'binomial' or 'bb' (beta-binomial)!"
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
  response <- y
  check_type(response, "numeric")
  if (obs_model != 1) {
    check_non_negative_all(response)
    check_integer_all(response)
  }
  return(TRUE)
}


#' Create the function that does a standardizing transform and its inverse
#'
#' @param y variable measurements
#' @param name variable name
#' @return an object of class \linkS4class{lgpscaling}
create_scaling <- function(y, name) {
  check_length_geq(y, 2)
  m <- mean(y)
  std <- stats::sd(y)
  if (std == 0) {
    msg <- paste0("the variable <", name, "> has zero variance!")
    stop(msg)
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
