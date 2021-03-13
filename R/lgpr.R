#' The 'lgpr' package.
#'
#' @description
#' \itemize{
#'   \item A package for Bayesian additive Gaussian process (GP) modeling with
#' interpretable kernels designed for longitudinal data.
#'   \item Features inference of covariate effects
#'   and covariate relevance assesment.
#'   \item Models can include non-stationary, heterogeneous and temporally
#' uncertain effects.
#'   \item Bayesian inference of the model (hyper)parameters using
#' \code{\link[rstan]{rstan}}.
#'   \item  Functions for visualizing longitudinal data,
#' posterior draws, model predictions and inferred covariate effects are
#' also provided.
#' }
#'
#' @author Juho Timonen (first.last at iki.fi)
#' @keywords Gaussian processes, longitudinal data, Stan, covariate relevances,
#' interpretable models
#'
#' @section Overview:
#' See the documentation of
#' \itemize{
#'  \item \code{\link{lgp}} for info on how to specify and fit models
#'  \item \code{\link{ppc}} for prior and posterior predictive checks
#'  \item \code{\link{relevances}} for component/covariate relevance assessment
#'  \item \code{\link{select}} for component/covariate selection
#'  \item \code{\link{pred}} for computing model predictions and
#'  inferred components
#'  \item \code{\link{plot_pred}} and \code{\link{plot_f}} for visualizing
#'  predictions and inferred components
#' }
#'
#' @section Tutorials:
#' See tutorials at \url{https://jtimonen.github.io/lgpr-usage/index.html}.
#'
#' @section Citation:
#' Run \code{citation("lgpr")} to get citation information.
#'
#' @section Feedback:
#' Bug reports, PRs, enchancement ideas or user experiences in general are
#' welcome and appreciated. Create an issue in Github or email the author.
#'
#' @docType package
#' @name lgpr-package
#' @aliases lgpr
#' @useDynLib lgpr, .registration = TRUE
#' @import methods rstantools
#' @importClassesFrom rstan stanfit
#'
#' @references
#' \enumerate{
#'   \item Carpenter, B. et al. (2017).
#'   \emph{Stan: A probabilistic programming language}. Journal of Statistical
#'    Software 76(1).
#'
#' }
#'
NULL
