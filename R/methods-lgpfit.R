#' Extract posterior draws
#'
#' @export
#' @description Uses \code{rstan::extract} with \code{permuted = FALSE} and
#' \code{inc_warmup = FALSE}.
#' @param fit an object of class \linkS4class{lgpfit}
#' @param ... other keyword arguments to \code{rstan::extract}
#' @return a named list
#' @family fit postprocessing functions
get_draws <- function(fit, ...) {
  check_type(fit, "lgpfit")
  rstan::extract(fit@stan_fit, permuted = FALSE, inc_warmup = FALSE, ...)
}

#' Extract posterior of the function f and its components
#'
#' @export
#' @param fit an object of class \linkS4class{lgpfit}
#' @param draws Indices of posterior draws for which to get \code{f}. This can
#' be a single integer, a vector of indices, or \code{NULL} (default). In the
#' latter case all draws are obtained.
#' @return Returns a list with names \code{num_draws} and \code{f}.
#' The latter is a named list of which has length equal to the number of
#' components plus one. Let \code{S = length(draws)}. Each list element is
#' \itemize{
#'   \item Array of size \code{S} x \code{num_obs}, where each row is one
#'   posterior draw of the function f, \code{is_sampled(model)} is \code{TRUE}.
#'   \item A list with fields \code{mean} and \code{std}, if
#'   \code{is_sampled(model)} is \code{FALSE}. Both fields are arrays of size
#'   \code{S} x \code{num_obs}. These are the analytically computed means and
#'   standard deviations for each posterior draw.
#' }
#' @family fit postprocessing functions
get_posterior_f <- function(fit, draws = NULL) {
  check_type(fit, "lgpfit")
  f_sampled <- is_f_sampled(fit)
  names <- get_component_names(fit)
  D <- length(names)
  R <- D + 1
  all_names <- c(names, "total")
  pars <- if (f_sampled) "f_latent" else "f_post"
  fp <- get_draws(fit, pars = pars)
  fp <- squeeze_second_dim(fp)
  S <- dim(fp)[1]
  if (is.null(draws)) {
    draws <- c(1:S)
  }
  if (!f_sampled) {
    alist <- array_to_arraylist(fp, 2 * R, draws)
    mean <- alist[1:R] # means
    std <- alist[(R + 1):(2 * R)] # stds
    f_out <- zip_lists(mean, std)
  } else {
    f_out <- array_to_arraylist(fp, R, draws)
  }
  names(f_out) <- all_names

  # Return
  list(
    f = f_out,
    num_draws = length(draws)
  )
}

#' Scale the function f posterior to original unnormalized scale
#'
#' @export
#' @description Can only be used with Gaussian observation model.
#' @param fit an object of class \linkS4class{lgpfit}
#' @param f_total a list with fields \code{mean} and \code{std}
#' @return a similar object as \code{f_total}
#' @family fit postprocessing functions
scale_f_total <- function(fit, f_total) {
  check_type(fit, "lgpfit")
  check_not_null(f_total)
  f_sampled <- is_f_sampled(fit)
  check_false(f_sampled)
  fun_inv <- fit@model@var_scalings$y@fun_inv
  f_total$mean <- scale_f_post_helper(fun_inv, f_total$mean)
  f_total$std <- scale_f_post_helper(fun_inv, f_total$std)
  return(f_total)
}

#' Posterior summary
#'
#' @export
#' @param fit an object of class \linkS4class{lgpfit}
#' @return a character representation
#' @family fit postprocessing functions
fit_summary <- function(fit) {
  check_type(fit, "lgpfit")
  print(fit@stan_fit, pars = c("f_post", "lp__"), include = FALSE)
}
