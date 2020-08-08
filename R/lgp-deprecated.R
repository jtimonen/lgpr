#' Use lgp() with the new syntax instead.
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
lgp_deprecated <- function(formula,
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
  stop("Use lgp() with the new syntax instead.")
}
