#' The 'lgpr' package.
#'
#' @description Bayesian nonparametric modeling and covariate selection for longitudinal data. 
#' The package enables 
#' \itemize{
#'   \item Additive Gaussian process modeling of longitudinal data
#'   \item Full posterior inference for the above models using No-U-Turn Hamiltonian Monte Carlo sampling
#'   \item Selecting the relevant covariates that explain the target variable
#'   \item Specialized modeling of a non-stationary disease effect
#'   \item Functions for visualizing longitudinal data, posterior samples and model predictions
#' }
#' @author Juho Timonen (first.last at aalto.fi)
#' @keywords Gaussian processes, Nonparametric modeling, Covariate selection, Longitudinal data, Stan
#'
#' @section Basic usage:
#' \itemize{
#' \item See the main function \code{\link{lgp}} for creating and fitting additive 
#'  longitudinal GP models.
#' \item Predictions outside the data can be computed using the function \code{\link[lgpr]{lgp_predict}}.
#' \item See documentation of the function \code{\link[lgpr]{simulate_data}} for generating
#'  artificial data.
#' \item For visualizing the data and results, see for example the functions 
#'   \itemize{
#'     \item \code{\link[lgpr]{plot_data}}
#'     \item \code{\link[lgpr]{plot_samples}}
#'     \item \code{\link[lgpr]{plot_components}}
#'     \item \code{\link[lgpr]{plot_predictions}}
#'     \item \code{\link[lgpr]{plot_simdata}}
#'   }
#' }
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
#'   \item Carpenter, B. et al. (2017). \emph{Stan: A probabilistic programming language}. Journal of Statistical Software 76(1).
#'   \item Gabry, J. and Goodrich, B. (2018). \emph{rstantools: Tools for Developing R Packages
#'   Interfacing with 'Stan'}. R package version 1.5.1.
#'   \item Gabry, J. and Mahr, T. (2018). \emph{bayesplot: Plotting for Bayesian Models}. R package version 1.6.0.
#'   \item Stan Development Team (2018). \emph{RStan: the R interface to Stan}. R package version 2.17.4. http://mc-stan.org
#' }
#'
NULL
