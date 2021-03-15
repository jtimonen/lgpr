#' The 'lgpr' package.
#'
#' @description Interpretable nonparametric modeling of longitudinal data
#' using additive Gaussian process regression. Contains functionality
#' for inferring covariate effects and assessing covariate relevances.
#' Models are specified using a convenient formula syntax, and can include
#' shared, group-specific, non-stationary, heterogeneous and temporally
#' uncertain effects. Bayesian inference for model parameters is performed
#' using Stan (\code{\link[rstan]{rstan}}).
#'
#' @author Juho Timonen (first.last at iki.fi)
#' @keywords Gaussian processes, longitudinal data, Stan, covariate relevances,
#' interpretable models
#'
#'
#' @section Specifying and fitting models:
#' Core functionality of the package consists of creating and fitting an
#' additive GP model. Main functions are:
#' \itemize{
#'  \item \code{\link{lgp}}: Specify and fit a model with one command.
#'  \item \code{\link{create_model}}: Define an \linkS4class{lgpmodel} object.
#'  \item \code{\link{sample_model}}: Sample model parameters and create an
#'  an \linkS4class{lgpfit} object.
#'  \item \code{\link{optimize_model}} (experimental): Optimize model
#'  parameters.
#' }
#'
#' @section Data:
#' The data that you wish to analyze with \code{lgpr} should be in an R
#' \code{data.frame} where columns correspond to measured variables and rows
#' corresponds to obervations. Some functions that can help working with such
#' data frames are:
#' \itemize{
#'  \item \code{\link{plot_data}}: Visualizing data.
#'  \item \code{\link{new_x}}: Creating new test points where the posterior
#'  distribution of any function component or sum of all components, or the
#'  posterior predictive distribution can be computed after model fitting.
#'  \item Other functions: \code{\link{add_factor}},
#'  \code{\link{add_factor_crossing}}, \code{\link{add_dis_age}},
#'  \code{\link{adjusted_c_hat}}.
#' }
#'
#' @section Tutorials and case studies:
#' See \url{https://jtimonen.github.io/lgpr-usage/index.html}.
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
