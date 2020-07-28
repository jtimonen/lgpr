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
#' @author Juho Timonen (first.last at aalto.fi)
#' @keywords Gaussian processes, longitudinal data, Stan, covariate selection,
#' interpretable models
#'
#' @section Basic usage:
#' \itemize{
#' \item See the main function \code{\link{lgp}} for creating and fitting
#' additive longitudinal GP models.
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
#' @importFrom rstan sampling
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
