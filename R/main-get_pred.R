#' Extract model predictions and function posteriors
#'
#' @export
#' @description \emph{NOTE:} It is not recommended for users to call this. Use
#' \code{\link{pred}} instead.
#' @inheritParams pred
#' @return an object of class \linkS4class{GaussianPrediction} or
#' \linkS4class{Prediction}
get_pred <- function(fit, draws = NULL, reduce = NULL, verbose = TRUE) {
  check_type(fit, "lgpfit")
  if (contains_postproc(fit)) {
    pr <- fit@postproc_results
    return(dollar(pr, "pred"))
  }
  in1 <- is.null(draws)
  in2 <- is.null(reduce)
  if (!in1 && !in2) stop("either draws or reduce (or both) must be NULL")
  if (is_f_sampled(fit)) {
    pred <- get_pred.sampled(fit, draws, reduce)
  } else {
    pred <- get_pred.gaussian(fit, draws, reduce, verbose)
  }
  return(pred)
}

# Get predictions for a model where the latent function f was marginalized
get_pred.gaussian <- function(fit, draws, reduce, verbose) {
  if (contains_postproc(fit)) {
    pred <- dollar(fit@postproc_results, "pred")
  } else {
    log_info(
      "No existing postprocessing information stored. Recomputing...",
      verbose
    )
  }
  pred <- pred(fit, reduce = reduce, draws = draws, verbose = verbose)
  return(pred)
}

# Get predictions for a model where the latent function f was sampled
get_pred.sampled <- function(fit, draws, reduce) {
  new("Prediction",
    f_comp = get_pred.sampled.f_comp(fit, draws, reduce),
    f = get_pred.sampled.f(fit, draws, reduce),
    h = get_pred.sampled.h(fit, draws, reduce),
    x = get_data(get_model(fit)),
    extrapolated = FALSE
  )
}

# Get the draws of each component of \code{f}
# Returns a list of arrays of shape \code{num_draws} x \code{num_obs}
get_pred.sampled.f_comp <- function(fit, draws, reduce) {
  nams <- component_names(fit)
  D <- length(nams)
  fp <- get_draws(fit, pars = "f_latent", draws = draws, reduce = reduce)
  fp <- array_to_arraylist(fp, D)
  names(fp) <- nams
  return(fp)
}

# Gets draws of the total \code{f}
# Returns an array of shape \code{num_draws} x \code{num_obs}
get_pred.sampled.f <- function(fit, draws, reduce) {
  f_comp <- get_pred.sampled.f_comp(fit, draws, reduce)
  f <- STAN_matrix_array_sum(f_comp, get_stream())
  return(f)
}

# Gets draws of the total \code{f}, adds \code{c_hat} to each draw, and
# then maps through the inverse link function. Applies \code{reduce} only
# after this transformation.
# Returns an array of shape \code{num_draws} x \code{num_obs}
get_pred.sampled.h <- function(fit, draws, reduce) {
  # Get f and add c_hat
  f <- get_pred.sampled.f(fit, draws, reduce = NULL)
  c_hat <- get_chat(fit)
  h <- map_f_to_h(fit, f, c_hat, reduce)
  return(h)
}
