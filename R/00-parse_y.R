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

#' Subroutines for parsing the response variable and its likelihood model
#'
#' @param Y_RAW raw response variable taken from the user-supplied data frame
#' @param LH integer encoding of likelihood
#' @inheritParams parse_y
#' @name parse_y_helper
NULL

#' @rdname parse_y_helper
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

#' @rdname parse_y_helper
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

#' Check that the response is numeric and compatible with observation model
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
