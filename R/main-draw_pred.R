#' Draw pseudo-observations from a predictive distribution
#'
#' @export
#' @param fit An object of class \linkS4class{lgpfit}.
#' @param pred An object of \linkS4class{Prediction}.
#' @family main functions
draw_pred <- function(fit, pred) {
  model <- get_model(fit)
  stan_fit <- get_stanfit(fit)
  draw_pred.subroutine(model, stan_fit, pred)
}

# Draw pseudo-observations from a predictive distribution
draw_pred.subroutine <- function(model, stan_fit, pred) {
  if (is(pred, "GaussianPrediction")) {
    msg <- paste0(
      "draw_pred() doesn't work for GaussianPrediction objects ",
      "because currently they don't hold information about the ",
      "covariance of the predicitive distribution."
    )
  }
  draw_pred.Prediction(model, stan_fit, pred)
}

# Draw pseudo-observations from a predictive distribution
draw_pred.Prediction <- function(model, stan_fit, pred) {
  stopifnot(is(pred, "Prediction"))
  obs_model <- get_obs_model(model)
  h <- pred@h
  S <- dim(h)[1] # number of draws
  N <- dim(h)[2] # number of points
  y_draws <- matrix(0.0, S, N)

  # Get possible observation model parameter draws
  if (obs_model == "gaussian") {
    theta_obs <- as.vector(get_draws(stan_fit, pars = "sigma"))
  } else if (obs_model == "nb") {
    theta_obs <- as.vector(get_draws(stan_fit, pars = "phi"))
  } else if (obs_model == "bb") {
    theta_obs <- as.vector(get_draws(stan_fit, pars = "gamma"))
  } else {
    theta_obs <- rep(NA, N)
  }

  # Get possible number of trials
  if (is_bin_or_bb(obs_model)) {
    num_trials <- as.vector(get_num_trials(model))
  } else {
    num_trials <- NULL
  }

  for (s in seq_len(S)) {
    h_s <- h[s, ] # vector of length N
    theta_s <- theta_obs[s] # single number
    if (obs_model == "gaussian") {
      y_s <- stats::rnorm(n = length(h_s), mean = h_s, sd = theta_s)
    } else if (obs_model == "poisson") {
      y_s <- stats::rpois(n = length(h_s), lambda = h_s)
    } else if (obs_model == "nb") {
      y_s <- stats::rnbinom(n = length(h_s), mu = h_s, size = theta_s)
    } else if (obs_model == "binomial") {
      y_s <- stats::rnbinom(n = length(h_s), prob = h_s, size = num_trials)
    } else if (obs_model == "bb") {
      gam_t <- (1 - gamma) / gamma
      alpha <- gam_t * h_s
      beta <- gam_t * (1.0 - h_s)
      prob <- stats::rbeta(n = length(alpha), alpha, beta)
      y_s <- stats::rbinom(n = length(alpha), prob = prob, size = num_trials)
    } else {
      stop("Unknown obs_model! Please report a bug.")
    }
    y_draws[s, ] <- y_s
  }
  return(y_draws)
}
