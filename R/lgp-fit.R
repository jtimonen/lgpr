
#' Fit an lgp model
#'
#' @export
#' @description Samples the posterior of an additive Gaussian process regression
#' model using \code{\link{rstan}}.
#' @param model An object of class \code{lgpmodel}.
#' @param ... Optional arguments passed to \code{rstan::sampling}, for example
#' \code{iter}, \code{chains} or \code{control}. See
#' \code{\link[rstan]{sampling}} for the possible arguments.
#' @param parallel This argument has been deprecated and throws and error if
#' \code{TRUE}. Use the \code{cores} argument of \code{rstan::sampling} instead.
#' @param skip_postproc In this mode the postprocessing after running Stan is
#' skipped.
#' @param threshold Component selection threshold for relevance sum.
#' @param relevance_method Component relevance determination method.
#' Must be either \code{"f_mean"} or \code{"alpha"}.
#' @param verbose should some output be printed?
#' @param ... Additional arguments that are passed to \code{rstan::sampling}.
#' @return An object of class \code{lgpfit}.
#' @seealso For all possible additional arguments (\code{...}), see
#' \code{?rstan::sampling}. For creating the \code{lgpmodel} input, see
#' \code{\link{lgp_model}}.
lgp_fit <- function(model,
                    threshold = 0.95,
                    parallel = FALSE,
                    skip_postproc = FALSE,
                    relevance_method = "f_mean",
                    verbose = FALSE,
                    ...) {

  # Set Stan auto_write to true
  rstan::rstan_options(auto_write = TRUE)

  # Check if parallelization should be used
  if (parallel) {
    stop(
      "The 'parallel' argument has been deprecated. Use the 'cores' argument ",
      "instead. See ?rstan::sampling for help."
    )
  }

  # Run stan
  stan_dat <- model@stan_dat
  stan_fit <- rstan::sampling(
    object = stanmodels[["lgp"]],
    data = stan_dat, ...
  )

  # Initialize the 'lgpfit' object
  ver <- "NA"
  tryCatch(
    {
      ver <- get_pkg_description()$Version
    },
    error = function(e) {
      # Do nothing
    }
  )
  fit <- new("lgpfit",
    stan_fit = stan_fit, model = model,
    pkg_version = ver
  )

  if (verbose) {
    cat("* Begin postprocessing. \n")
  }
  fit@diagnostics <- assess_convergence(fit, skip_F_gen = TRUE)

  # Finalize the 'lgpfit' object
  tryCatch(
    {
      if (!skip_postproc) {
        fit <- postproc(fit,
          threshold = threshold,
          relevance_method = relevance_method,
          verbose = verbose
        )
        if (verbose) {
          show(fit)
        }
      } else {
        if (verbose) {
          cat("* Skipped postprocessing.\n")
        }
      }
    },
    error = function(e) {
      warning(e)
      cat("* Postprocessing failed, the error is printed as a warning.\n")
    }
  )

  # Return
  return(fit)
}
