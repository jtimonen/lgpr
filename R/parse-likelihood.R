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
