#' Main function of the `lgpr` package
#'
#' @export
#' @description
#' \itemize{
#'   \item The \code{lgp} function is a wrapper for both \code{lgp_model}
#'   and \code{lgp_fit}.
#'   \item The \code{lgp_model} function creates an object of class
#'   \linkS4class{lgpmodel}.
#'   \item The \code{lgp_fit} function takes an \linkS4class{lgpmodel} object
#'   and fits the model using \code{rstan::sampling}.
#' }
#' @inheritParams lgp_model
#' @inheritParams lgp_fit
#' @return a fit object
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
  fit <- lgp_fit(model = model, ...)
  list(model = model, fit = fit)
}
