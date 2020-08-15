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
#' binomial) or \code{"binomial"}. To use Bernoulli likelihood, use
#' \code{likelihood="binomial"} and set \code{num_trials} as a vector of ones.
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
    Y_NORM <- normalizer@fun(Y_RAW) # standardize the response
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
#' likelihood is binomial. Must have length one or equal to number of data
#' points. Setting \code{num_trials=1} corresponds to Bernoulli observation
#' model.
#' @param list_y a list field returned by \code{\link{parse_response}}
#' @return a list of parsed options
parse_likelihood <- function(likelihood, c_hat, num_trials, list_y) {
  LH <- likelihood_as_int(likelihood)
  y <- if (LH != 1) list_y$y_disc else list_y$y_cont
  y <- as.numeric(y)
  num_obs <- length(y)
  num_trials <- set_num_trials(num_trials, num_obs, LH)
  c_hat <- set_c_hat(c_hat, y, LH, num_trials)
  list(
    obs_model = LH,
    y_num_trials = num_trials,
    c_hat = c_hat,
    is_likelihood_skipped = 0
  )
}

#' Set c_hat (non-gaussian observation models)
#'
#' @param c_hat the \code{c_hat} argument given as input to
#' \code{\link{lgp_model}}
#' @param response response variable
#' @param LH likelihood as int
#' @param num_trials the num_trials data (binomial likelihood)
#' @return a real number
set_c_hat <- function(c_hat, response, LH, num_trials) {
  nb_or_pois <- LH %in% c(2, 3)
  binomial <- LH == 4
  if (is.null(c_hat)) {
    if (nb_or_pois) {
      c_hat <- log(mean(response))
    } else if (binomial) {
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

#' Set num_trials (binomial and Bernoulli observation models)
#'
#' @param num_trials the \code{num_trials} argument
#' @param num_obs number of observations
#' @param LH likelihood as integer encoding
#' @return a numeric vector
set_num_trials <- function(num_trials, num_obs, LH) {
  if (is.null(num_trials)) {
    num_trials <- rep(1, num_obs)
  } else {
    if (LH != 4) {
      stop("Only give the <num_trials> argument if likelihood is binomial!")
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
