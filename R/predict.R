#' Extract model predictions
#'
#' @inheritParams get_draws
#' @return an object of class \linkS4class{GaussianPrediction} or
#' \linkS4class{Prediction}
#' @name get_pred
#' @family prediction extraction functions
NULL

#' @export
#' @rdname get_pred
get_pred <- function(fit, draws = NULL, reduce = NULL) {
  check_type(fit, "lgpfit")
  in1 <- is.null(draws)
  in2 <- is.null(reduce)
  if (!in1 && !in2) stop("either draws or reduce (or both) must be NULL")
  if (is_f_sampled(fit)) {
    pred <- get_pred.sampled(fit, draws, reduce)
  } else {
    pred <- get_pred.gaussian(fit, draws, reduce)
  }
  return(pred)
}

#' @rdname get_pred
get_pred.gaussian <- function(fit, draws, reduce) {
  f_comp <- get_pred.gaussian.f_comp(fit, draws, reduce)
  f <- get_pred.gaussian.f(fit, draws, reduce)
  y <- get_pred.gaussian.y(fit, draws, reduce)
  new("GaussianPrediction",
    f_comp_mean = dollar(f_comp, "mean"),
    f_comp_std = dollar(f_comp, "std"),
    f_mean = dollar(f, "mean"),
    f_std = dollar(f, "std"),
    y_mean = dollar(y, "mean"),
    y_std = dollar(y, "std")
  )
}

#' @rdname get_pred
get_pred.sampled <- function(fit, draws, reduce) {
  new("Prediction",
    f_comp = get_pred.sampled.f_comp(fit, draws, reduce),
    f = get_pred.sampled.f(fit, draws, reduce),
    h = get_pred.sampled.h(fit, draws, reduce)
  )
}

#' Get analytic posterior or predictive distributions
#'
#' @description These are helper functions for \code{\link{get_pred.gaussian}}.
#' \itemize{
#'   \item Function \code{get_pred.gaussian.f_comp} gets the componentwise
#'   posterior distributions of each component of \code{f} on the normalized
#'   scale.
#'   \item Function \code{get_pred.gaussian.f} gets the total posterior
#'   distribution of \code{f} on the normalized scale.
#'   \item Function \code{get_pred.gaussian.y} gets the predictive distribution
#'   on the unnormalized original data scale. It applies \code{reduce} only
#'   after adding sigma^2 to the variance of \code{f}.
#' }
#' @inheritParams get_pred.gaussian
#' @return a list with names \code{mean} and \code{std}
#' @name get_pred_gaussian
#' @family prediction extraction functions
NULL

#' @rdname get_pred_gaussian
get_pred.gaussian.f_comp <- function(fit, draws, reduce) {
  nams <- c(get_component_names(fit), "total")
  R <- length(nams)
  fp <- get_draws(fit, pars = "f_post", draws = draws, reduce = reduce)
  alist <- array_to_arraylist(fp, 2 * R)
  mu <- alist[1:R] # means
  std <- alist[(R + 1):(2 * R)] # stds
  names(mu) <- nams
  names(std) <- nams

  # Return
  list(
    mean = mu[1:(R - 1)],
    std = std[1:(R - 1)]
  )
}

#' @rdname get_pred_gaussian
get_pred.gaussian.f <- function(fit, draws, reduce) {
  nams <- c(get_component_names(fit), "total")
  R <- length(nams)
  fp <- get_draws(fit, pars = "f_post", draws = draws, reduce)
  alist <- array_to_arraylist(fp, 2 * R)
  mu <- alist[1:R] # means
  std <- alist[(R + 1):(2 * R)] # stds
  names(mu) <- nams
  names(std) <- nams
  list(
    mean = dollar(mu, "total"),
    std = dollar(std, "total")
  )
}

#' @rdname get_pred_gaussian
get_pred.gaussian.y <- function(fit, draws, reduce) {

  # Get posterior of f and sigma on normalized scale without any reduction
  f <- get_pred.gaussian.f(fit, draws, reduce = NULL)
  sigma <- get_draws(fit, pars = "sigma", draws = draws, reduce = NULL)
  sigma <- as.vector(sigma)

  # Add sigma to get y_pred
  y_mean <- dollar(f, "mean")
  fstd <- dollar(f, "std")
  y_var <- add_to_columns(fstd^2, sigma^2)
  y_std <- sqrt(y_var)

  # Scale y_pred to original scale
  y_scl <- dollar(fit@model@var_scalings, "y")
  y_mean <- call_fun(y_scl@fun_inv, y_mean)
  y_std <- call_fun(y_scl@fun_inv, y_std)

  # Apply reduce and return
  list(
    mean = apply_reduce(y_mean, reduce),
    std = apply_reduce(y_std, reduce)
  )
}

#' Get predictions for a model where the latent function f was sampled
#'
#' @description These are helper functions for \code{\link{get_pred.sampled}}.
#' \itemize{
#'   \item Function \code{get_pred.sampled.f_comp} gets the draws of each
#'   component of \code{f} on the normalized scale..
#'   \item Function \code{get_pred.sampled.f} gets draws of the total \code{f}
#'   on the normalized scale.
#'   \item Function \code{get_pred.sampled.h} gets draws of the total \code{f}.
#'   maps them to unnormalized scale and through the inverse link function.
#'   It applies \code{reduce} only after this transformation.
#' }
#' @inheritParams get_pred.sampled
#' @return an array of shape \code{num_draws} x \code{num_obs} (actually a
#' list of such arrays for \code{get_pred.sampled.f_comp})
#' @name get_pred_sampled
#' @family prediction extraction functions
NULL

#' @rdname get_pred_sampled
get_pred.sampled.f_comp <- function(fit, draws, reduce) {
  nams <- get_component_names(fit)
  D <- length(nams)
  fp <- get_draws(fit, pars = "f_latent", draws = draws, reduce = reduce)
  fp <- array_to_arraylist(fp, D)
  names(fp) <- nams
  return(fp)
}

#' @rdname get_pred_sampled
get_pred.sampled.f <- function(fit, draws, reduce) {
  f_comp <- get_pred.sampled.f_comp(fit, draws, reduce)
  f <- STAN_matrix_array_sum(f_comp, get_stream())
  return(f)
}

#' @rdname get_pred_sampled
get_pred.sampled.h <- function(fit, draws, reduce) {

  # Get f
  f <- get_pred.sampled.f(fit, draws, reduce = NULL)

  # Scale to original unnormalized scale
  y_scl <- dollar(fit@model@var_scalings, "y")
  f <- call_fun(y_scl@fun_inv, f)

  # Apply inverse link function and reduction
  likelihood <- get_obs_model(fit)
  h <- link_inv(f, likelihood)
  h <- apply_reduce(h, reduce)
  return(h)
}
