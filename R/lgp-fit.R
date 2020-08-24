#' Fit a model
#'
#' @export
#' @param model An object of class \linkS4class{lgpmodel}.
#' @param ... Optional arguments passed to \code{rstan::sampling}, for example
#' \code{iter}, \code{chains} or \code{control}. See
#' \code{\link[rstan]{sampling}} for the possible arguments.
lgp_fit <- function(model, ...) {
  rstan::sampling(object = stanmodels[["lgp"]], data = model@stan_input, ...)
}
