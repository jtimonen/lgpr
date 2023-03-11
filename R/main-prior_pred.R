#' Prior (predictive) sampling
#'
#' @description
#' These functions take an \linkS4class{lgpmodel} object, and
#' \itemize{
#'   \item \code{prior_pred} samples from the prior predictive distribution of
#'   the model
#'   \item \code{sample_param_prior} samples only its parameter prior using
#'   \code{\link[rstan]{sampling}}
#' }
#' @inheritParams pred
#' @param model An object of class \linkS4class{lgpmodel}.
#' @param refresh Argument for \code{\link[rstan]{sampling}}.
#' @param quiet This forces \code{verbose} to be \code{FALSE}. If you want
#' to suppress also the output from Stan, give the additional argument
#' \code{refresh=0}.
#' @param ... Additional arguments for \code{\link[rstan]{sampling}}.
#' @name prior_pred
#' @family main functions
#' @return
#' \itemize{
#'   \item \code{prior_pred} returns a list with components
#'   \itemize{
#'      \item \code{y_draws}: A matrix containing the prior predictive draws
#'      as rows. Can be passed to \code{bayesplot::pp_check()} for
#'      graphical prior predictive checking.
#'      \item \code{pred_draws}: an object of class \linkS4class{Prediction},
#'      containing prior draws of each model component and their sum
#'      \item \code{param_draws}: a \code{stanfit} object of prior parameter
#'      draws (obtained by calling \code{sample_param_prior} internally)
#'    }
#'   \item \code{sample_param_prior} returns
#'   an object of class \code{\link[rstan]{stanfit}}
#' }
NULL

#' @rdname prior_pred
#' @export
prior_pred <- function(model,
                       verbose = TRUE,
                       quiet = FALSE,
                       refresh = 0,
                       STREAM = get_stream(),
                       ...) {
  check_type(model, "lgpmodel")
  if (quiet) verbose <- FALSE
  log_progress("Sampling from parameter prior...", verbose)
  stan_fit <- sample_param_prior(
    model,
    verbose = verbose,
    quiet = quiet,
    refresh = refresh,
    ...
  )
  log_progress("Drawing from GP priors conditional on parameters...", verbose)
  f_draws <- draw_f_prior(model, stan_fit, verbose, STREAM)
  log_progress("Drawing from prior predictive distribution...", verbose)
  draw_prior_pred(model, stan_fit, f_draws, verbose)
}


#' @rdname prior_pred
#' @export
sample_param_prior <- function(model, verbose = TRUE, quiet = FALSE, ...) {
  check_type(model, "lgpmodel")
  if (quiet) verbose <- FALSE
  object <- dollar(stanmodels, "parameter_prior")
  data <- model@stan_input

  # Run parameter prior sampling
  rstan::sampling(
    object = object,
    data = data,
    check_data = TRUE,
    ...
  )
}


# Draw from prior predictive distribution given prior draws of parameters
# and functions
draw_prior_pred <- function(model, stan_fit, f_draws, verbose) {
  # Create h, with shape num_draws x num_points
  f_sum_draws <- dollar(f_draws, "f_draws")
  c_hat <- get_chat(model)
  h <- map_f_to_h(model, f_sum_draws, c_hat, reduce = NULL)

  # Create Prediction
  pred_draws <- new("Prediction",
    f_comp = dollar(f_draws, "f_draws_comp"),
    f = f_sum_draws,
    h = h,
    x = dollar(f_draws, "x"),
    extrapolated = FALSE
  )

  # Draw from predictive distribution and return
  y_draws <- draw_pred.subroutine(model, stan_fit, pred_draws)
  list(
    y_draws = y_draws,
    pred_draws = pred_draws,
    param_draws = stan_fit
  )
}

# Draw from additive GP prior given kernel and other parameter draws
#
# @param model An object of class \linkS4class{lgpmodel}.
# @param stan_fit An object of class \code{stanfit}, containing
# prior parameter draws.
# @param x Data frame of points where to draw the GP values.
# @inheritParams pred
# @return A named list.
draw_f_prior <- function(model, stan_fit, verbose, STREAM) {
  kc <- create_kernel_computer(model, stan_fit, NULL, NULL, NULL, FALSE, STREAM)
  input <- kc@input
  K_input <- kc@K_input
  S <- num_paramsets(kc) # number of parameter sets
  P <- num_evalpoints(kc) # number of output points
  J <- num_components(kc) # number of components
  comp_names <- component_names(kc)
  delta <- dollar(input, "delta")

  # Create output arrays
  f_draws_comp <- array(0.0, c(S, P, J))
  f_draws <- array(0.0, c(S, P))

  # Setup
  progbar <- verbose && S > 1
  pb <- progbar_setup(L = S)
  hdr <- dollar(pb, "header")
  idx_print <- dollar(pb, "idx_print")
  log_progress(hdr, progbar)

  # Loop through parameter sets
  for (idx in seq_len(S)) {
    # Compute full kernel matrices for one parameter set
    K_prior <- kernel_all(K_input, input, idx, kc@STREAM)

    # Draw each component
    f_i <- draw_gp_components(K_prior, delta) # dim = (P, J)

    # Store result and update progress
    f_draws_comp[idx, , ] <- f_i
    f_draws[idx, ] <- rowSums(f_i)

    if (progbar) progbar_print(idx, idx_print)
  }
  log_progress(" ", progbar)

  # Return
  f_draws_comp <- aperm(f_draws_comp, c(3, 1, 2)) # dim (S, P, J) -> (J, S, P)
  out <- list(
    f_draws_comp = arr3_to_list(f_draws_comp, comp_names), # list with len J
    f_draws = f_draws, # dim (S, P)
    x = get_data(model)
  )
  return(out)
}

# Draw components from zero-mean GP priors
draw_gp_components <- function(K, delta) {
  N <- nrow(K[[1]])
  J <- length(K)
  mu0 <- rep(0.0, N)
  f_draws <- matrix(0.0, N, J)
  for (j in seq_len(J)) {
    Sigma <- K[[j]] + delta * diag(N)
    f_draws[, j] <- MASS::mvrnorm(n = 1, mu = mu0, Sigma = Sigma)
  }
  return(f_draws)
}
