#' Mains function of the `lgpr` package
#'
#' @export
#' @description
#' \itemize{
#'   \item The \code{lgp} function is a wrapper for both \code{lgp_model}
#'   and \code{lgp_sampling}, which returns an object of class
#'   \linkS4class{lgpfit}.
#'   \item The \code{\link{lgp_model}} function creates an object of class
#'   \linkS4class{lgpmodel}.
#'   \item The \code{lgp_sampling} function takes an \linkS4class{lgpmodel}
#'   object and fits the model using \code{rstan::sampling}.
#'   \item The \code{lgp_optimizing} function takes an \linkS4class{lgpmodel}
#'   object and fits the model using \code{rstan::optimizing}.
#' }
#' @inheritParams lgp_model
#' @inheritParams lgp_sampling
#' @inheritParams lgp_optimizing
#' @name lgp
NULL

#' @rdname lgp
lgp <- function(formula,
                data,
                likelihood = "gaussian",
                prior = NULL,
                c_hat = NULL,
                num_trials = NULL,
                options = NULL,
                threshold = 0.95,
                relevance_method = "f_mean",
                verbose = FALSE,
                ...) {
  
  # Create and fit the model
  model <- lgp_model(
    formula, data, likelihood, prior, c_hat, num_trials, options,
    verbose
  )
  stan_fit <- lgp_sampling(model = model, ...)
  new('lgpfit', model = model, stan_fit = stan_fit)
}

#' @rdname lgp
#' @export
#' @param model An object of class \linkS4class{lgpmodel}.
#' @param ... Optional arguments passed to \code{rstan::sampling}, for example
#' \code{iter}, \code{chains} or \code{control}. See
#' \code{\link[rstan]{sampling}} for the possible arguments.
lgp_sampling <- function(model, ...) {
  rstan::sampling(object = stanmodels[["lgp"]], data = model@stan_input, ...)
}

#' @rdname lgp
#' @export
#' @param model An object of class \linkS4class{lgpmodel}.
#' @param ... Optional arguments passed to \code{rstan::optimizing}, such as
#' \code{algorithm}. See \code{\link[rstan]{optimizing}} for the possible
#' arguments.
lgp_optimizing <- function(model, ...) {
  rstan::optimizing(object = stanmodels[["lgp"]], data = model@stan_input, ...)
}
