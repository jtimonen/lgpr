#' The 'lgpr' package.
#'
#' @description Longitudinal Gaussian Process regression.
#' The package features additive Gaussian process modeling of longitudinal
#' data with with interpretable covariate effects and covariate relevance
#' assesment. Models can include non-stationary, heterogeneous and temporally
#' uncertain effects. Bayesian inference of the model (hyper)parameters using
#' \code{\link[rstan]{rstan}}. Functions for visualizing longitudinal data,
#' posterior draws and model predictions are also provided.
#'
#' @author Juho Timonen (first.last at iki.fi)
#' @keywords Gaussian processes, longitudinal data, Stan, covariate relevances,
#' interpretable models
#'
#' @section Basic usage:
#' \itemize{
#' \item See documentation of the main function \code{\link{lgp}}.
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

#' Main function of the lgpr package
#'
#' @export
#' @description
#' Creates a model by calling \code{\link{create_model}} and fits by calling
#' \code{\link{sample_model}}.
#' @inheritParams create_model
#' @inheritParams sample_model
#' @inheritParams optimize_model
#'
#' @section Model formula syntax:
#' There are two ways to define the model formula:
#' \enumerate{
#'   \item Using a common \code{\link[stats]{formula}}-like syntax, like in
#'   \code{y ~ age +} \code{age|id} \code{ + sex}. Terms can consits of a
#'   single variable, such as \code{age}, or an interaction of two variables,
#'   such as \code{age|id}. In single-variable terms, the variable can be either
#'   continuous (numeric) or categorical (factor), whereas in interaction terms
#'   the variable on the left-hand side of the vertical bar (\code{|}) has to
#'   be continuous and the one on the right-hand side has to be categorical.
#'   Formulae specified using this syntax are translated to the advanced format
#'   so that
#'   \itemize{
#'     \item single-variable terms become \code{gp(x)} if
#'     variable \code{x} is numeric and \code{zs(x)} if \code{x} is a factor
#'     \item interaction terms \code{x|z} become \code{gp(x)*zs(z)}
#'   }
#'   \item Using the advanced syntax, like in \code{y ~ gp(age) +}
#'   \code{gp(age)*zs(id) +} \code{het(id)*gp_vm(disAge)}.
#'   This creates \linkS4class{lgprhs} objects, which consist of
#'  \linkS4class{lgpterm}s, which consist of \linkS4class{lgpexpr}s.
#'  This approach must be used if creating nonstationary, heterogeneous or
#'  temporally uncertain components.
#' }
#' Either one of the approaches should be used and they should not be mixed.
#'
#' @section Defining priors:
#' The \code{priors} argument must be a named list.
#' \enumerate{
#'   \item Ok.
#' }
#'
#' @name lgp
#' @family main functions
#' @return Returns an object of the S4 class \linkS4class{lgpfit}.
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
                ...) {

  # Create and fit the model
  model <- create_model(
    formula, data, likelihood, prior, c_hat, num_trials, options,
    prior_only, verbose, sample_f
  )
  sample_model(model = model, ...)
}

#' Fitting a model
#'
#' @description
#' \itemize{
#'   \item The \code{sample_model} function takes an \linkS4class{lgpmodel}
#'   object, fits it using  \code{\link[rstan]{sampling}}, and
#'   returns an object of class \linkS4class{lgpfit}.
#'   \item The \code{optimize_model} function takes an \linkS4class{lgpmodel}
#'   object and fits it using \code{\link[rstan]{optimizing}}.
#' }
#' @name sample_model
#' @family main functions
NULL

#' @rdname sample_model
#' @export
#' @param model An object of class \linkS4class{lgpmodel}.
#' @param ... Optional arguments passed to
#' \code{\link[rstan]{sampling}} or \code{\link[rstan]{optimizing}}.
sample_model <- function(model, ...) {
  object <- stanmodels[[model@stan_model_name]]
  data <- model@stan_input
  sfit <- rstan::sampling(object = object, data = data, check_data = TRUE, ...)
  new("lgpfit", model = model, stan_fit = sfit)
}

#' @rdname sample_model
#' @export
optimize_model <- function(model, ...) {
  object <- stanmodels[[model@stan_model_name]]
  data <- model@stan_input
  rstan::optimizing(object = object, data = data, check_data = TRUE, ...)
}

#' Set the GP mean vector, taking TMM or other normalization
#' into account
#'
#' @export
#' @description Creates the \code{c_hat} input for \code{lgp},
#' so that it accounts for normalization between data points in the
#' \code{"poisson"} or \code{"nb"} observation model
#' @param y response variable, vector of length \code{n}
#' @param norm_factors normalization factors, vector of length \code{n}
#' @return a vector of length \code{n}, which can be used as
#' the \code{c_hat} input to the \code{lgp} function
adjusted_c_hat <- function(y, norm_factors) {
  L1 <- length(norm_factors)
  L2 <- length(y)
  if (L1 != L2) stop("inputs must have same length!")
  if (sum(y < 0) > 0) stop("y cannot have negative values!")
  if (sum(round(y) != y) > 0) stop("y must have only integer values!")
  if (sum(norm_factors <= 0) > 0) stop("norm_factors must be all positive!")
  c_hat <- log(mean(y))
  c_hat <- rep(c_hat, length(y))
  c_hat <- c_hat + log(norm_factors)
  return(c_hat)
}
