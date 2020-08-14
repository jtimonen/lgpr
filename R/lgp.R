#' The main function of the `lgpr` package
#'
#' @export
#' @description This is a wrapper for both \code{\link{lgp_model}}
#'  and \code{\link{lgp_fit}}. It first creates an
#'  \linkS4class{lgpmodel} object and then fits the model,
#'  finally returning an \linkS4class{lgpfit} object.
#' @inheritParams lgp_model
#' @inheritParams lgp_fit
#' @return An object of class \code{lgpfit}.
lgp <- function(formula,
                data,
                likelihood = "gaussian",
                prior = prior_default(),
                c_hat = NULL,
                num_trials = NULL,
                id_variable = "id",
                options = NULL,
                disease_options = NULL,

                parallel = FALSE,
                threshold = 0.95,
                relevance_method = "f_mean",
                ...) {

  # Create the model
  model <- lgp_model(
    formula, data, likelihood, prior, c_hat, num_trials,
    id_variable, options, disease_options
  )

  # Fit the model
  # fit <- lgp_fit(
  #  model = model,
  #  threshold = threshold,
  #  parallel = parallel,
  #  relevance_method = relevance_method,
  #  ...
  # )
  fit <- model

  # Return the lgpfit object
  return(fit)
}
