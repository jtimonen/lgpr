#' The main function of the `lgpr` package
#'
#' @export
#' @description This is a wrapper for both \code{\link{lgp_model}}
#'  and \code{\link{lgp_fit}}. It first creates an
#'  \code{lgpmodel} object and then fits the model,
#'  finally returning an \code{lgpfit} object. Note that the
#'  covariate types are automatically inferred from the given
#'  \code{data}. If you wish to change these, see the
#'  arguments
#'  \itemize{
#'     \item \code{id_variable}
#'     \item \code{time_variable}
#'     \item \code{disAge_variable}
#'     \item \code{continuous_vars} and
#'     \item \code{categorical_vars}.
#'  }
#' @inheritParams lgp_model
#' @inheritParams lgp_fit
#' @return An object of class \code{lgpfit}.
lgp <- function(formula,
                data,
                likelihood = "Gaussian",
                prior = prior_default(),
                uncertain_effect_time = FALSE,
                equal_effect = TRUE,
                id_variable = "id",
                time_variable = "age",
                disAge_variable = NULL,
                continuous_vars = NULL,
                categorical_vars = NULL,
                offset_vars = NULL,
                C_hat = NULL,
                DELTA = 1e-8,
                sample_F = NULL,
                parallel = FALSE,
                skip_postproc = FALSE,
                threshold = 0.95,
                variance_mask = TRUE,
                N_trials = NULL,
                relevance_method = "f_mean",
                verbose = FALSE,
                ...) {
  # Create the model
  model <- lgp_model(
    formula = formula,
    data = data,
    likelihood = likelihood,
    prior = prior,
    uncertain_effect_time = uncertain_effect_time,
    equal_effect = equal_effect,
    C_hat = C_hat,
    DELTA = DELTA,
    sample_F = sample_F,
    id_variable = id_variable,
    time_variable = time_variable,
    disAge_variable = disAge_variable,
    continuous_vars = continuous_vars,
    categorical_vars = categorical_vars,
    offset_vars = offset_vars,
    variance_mask = variance_mask,
    N_trials = N_trials,
    skip_gen_quant = skip_postproc,
    verbose = verbose
  )

  if (verbose) {
    show(model)
  }

  # Fit the model
  fit <- lgp_fit(
    model = model,
    threshold = threshold,
    parallel = parallel,
    skip_postproc = skip_postproc,
    relevance_method = relevance_method,
    verbose = verbose,
    ...
  )

  # Return the lgpfit object
  return(fit)
}
