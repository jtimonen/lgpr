#' Draw pseudo-observations from a predictive distribution
#'
#' @description This doesn't work for \linkS4class{GaussianPrediction} objects
#' because for memory reasons they don't store the
#' covariance matrices of the predictive distributions.
#' @export
#' @param fit An object of class \linkS4class{lgpfit} that can been created
#' with \code{sample_f=TRUE}.
#' @param pred An object of \linkS4class{Prediction}.
#' @family main functions
draw_pred <- function(fit, pred = NULL) {
  if (!is_f_sampled(fit)) {
    stop("fit has been created with sample_f = FALSE")
  }
  if (is.null(pred)) {
    pred <- get_pred(fit)
  }
  model <- get_model(fit)
  stan_fit <- get_stanfit(fit)
  draw_pred.subroutine(model, stan_fit, pred)
}

# Draw pseudo-observations from a predictive distribution
draw_pred.subroutine <- function(model, stan_fit, pred) {
  if (!is(pred, "Prediction")) {
    stop("pred must be an object of class Prediction")
  }
  draw_pred.Prediction(model, stan_fit, pred)
}

# Draw pseudo-observations from a predictive distribution
draw_pred.Prediction <- function(model, stan_fit, pred) {
  stopifnot(is(pred, "Prediction"))
  obs_model <- get_obs_model(model)
  h <- pred@h
  S <- dim(h)[1] # number of draws
  L <- num_evalpoints(pred) # number of evaluation points
  y_draws <- matrix(0.0, S, L)

  # Get possible observation model parameter draws
  if (obs_model == "gaussian") {
    theta_obs <- as.vector(get_draws(stan_fit, pars = "sigma"))
  } else if (obs_model == "nb") {
    theta_obs <- as.vector(get_draws(stan_fit, pars = "phi"))
  } else if (obs_model == "bb") {
    theta_obs <- as.vector(get_draws(stan_fit, pars = "gamma"))
  } else {
    theta_obs <- rep(NA, S)
  }

  # Get possible number of trials
  if (is_bin_or_bb(obs_model)) {
    num_trials <- as.vector(get_num_trials(model))
    if (length(num_trials) != L) {
      stop("length(num_trials) not equal to num_evalpoints(pred)!")
    }
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
      gammma <- theta_s # number
      gam_t <- (1 - gamma) / gamma # number
      alpha <- gam_t * h_s # vector
      beta <- gam_t * (1.0 - h_s) # vector
      prob <- stats::rbeta(n = length(alpha), alpha, beta)
      y_s <- stats::rbinom(n = length(alpha), prob = prob, size = num_trials)
    } else {
      stop("Unknown obs_model! Please report a bug.")
    }
    y_draws[s, ] <- y_s
  }
  return(y_draws)
}
