#' The 'lgpr' package.
#'
#' @description Longitudinal Gaussian Process regression.
#' The package features
#' \itemize{
#'   \item Additive Gaussian process modeling of longitudinal data
#'   \item Posterior inference of the model (hyper)parameters using Stan
#'   \item Computation of covariate relevances
#'   \item Specialized modeling of a non-stationary disease effect
#'   \item Functions for visualizing longitudinal data, posterior samples and
#'   model predictions
#'   \item Gaussian, Poisson, binomial or negative binomial observation models
#' }
#' @author Juho Timonen (first.last at iki.fi)
#' @keywords Gaussian processes, longitudinal data, Stan, covariate selection,
#' interpretable models
#'
#' @section Basic usage:
#' \itemize{
#' \item See the main function documentation: \code{\link{lgp}}.
#'  \item See tutorials at
#'  \url{https://jtimonen.github.io/lgpr-usage/index.html}
#' }
#' @section Citation:
#'
#' \emph{An interpretable probabilistic machine learning method for
#' heterogeneous longitudinal studies}. Juho Timonen, Henrik Mannerstrom,
#' Aki Vehtari and Harri Lahdesmaki, 2019.
#' \url{https://arxiv.org/abs/1912.03549}
#'
#' @docType package
#' @name lgpr-package
#' @aliases lgpr
#' @useDynLib lgpr, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan get_rng get_stream
#'
#' @references
#' \enumerate{
#'   \item Carpenter, B. et al. (2017).
#'   \emph{Stan: A probabilistic programming language}. Journal of Statistical
#'    Software 76(1).
#'   \item Jonah Gabry, Ben Goodrich and Martin Lysy (2019).
#'   \emph{rstantools: Tools for Developing R Packages Interfacing with 'Stan'}.
#'   R package version 2.0.0.
#'   \item Gabry, J. and Mahr, T. (2019).
#'   \emph{bayesplot: Plotting for Bayesian Models}. R package version 1.7.0,
#'    http://mc-stan.org/bayesplot.
#'   \item Stan Development Team (2019). \emph{RStan: the R interface to Stan.}
#'   R package version 2.19.2. http://mc-stan.org/.
#' }
#'
NULL

#' Main functions of the package
#'
#' @export
#' @description
#' \itemize{
#'   \item The \code{lgp} function is a wrapper for both \code{create_model}
#'   and \code{sample_model}, and its returns an object of class
#'   \linkS4class{lgpfit}.
#'   \item The \code{\link{create_model}} function creates an object of class
#'   \linkS4class{lgpmodel}.
#'   \item The \code{sample_model} function takes an \linkS4class{lgpmodel}
#'   object, fits the model using \code{rstan::sampling}, and returns an object
#'   of class \linkS4class{lgpfit}.
#'   \item The \code{optimize_model} function takes an \linkS4class{lgpmodel}
#'   object and fits the model using \code{rstan::optimizing}.
#' }
#' @inheritParams create_model
#' @inheritParams sample_model
#' @inheritParams optimize_model
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
                prior_only = FALSE,
                verbose = FALSE,
                sample_f = !(likelihood == "gaussian"),
                id_variable = "id",
                time_variable = "age",
                ...) {

  # Create and fit the model
  model <- create_model(
    formula, data, likelihood, prior, c_hat, num_trials, options,
    prior_only, verbose, sample_f, id_variable, time_variable
  )
  sample_model(model = model, ...)
}

#' @rdname lgp
#' @export
#' @param model An object of class \linkS4class{lgpmodel}.
#' @param ... Optional arguments passed to \code{rstan::sampling}, for example
#' \code{iter}, \code{chains} or \code{control}. See
#' \code{\link[rstan]{sampling}} for the possible arguments.
sample_model <- function(model, ...) {
  object <- stanmodels[[model@stan_model_name]]
  data <- model@stan_input
  sfit <- rstan::sampling(object = object, data = data, check_data = TRUE, ...)
  new("lgpfit", model = model, stan_fit = sfit)
}

#' @rdname lgp
#' @export
#' @param model An object of class \linkS4class{lgpmodel}.
#' @param ... Optional arguments passed to \code{rstan::optimizing}, such as
#' \code{algorithm}. See \code{\link[rstan]{optimizing}} for the possible
#' arguments.
optimize_model <- function(model, ...) {
  object <- stanmodels[[model@stan_model_name]]
  data <- model@stan_input
  rstan::optimizing(object = object, data = data, check_data = TRUE, ...)
}
