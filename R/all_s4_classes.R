#' @include validate.R
NULL

#' An S4 class to represent an lgp expression
#'
#' @slot covariate name of a covariate
#' @slot fun function name
lgpexpr <- setClass("lgpexpr",
  representation = representation(
    covariate = "character",
    fun = "character"
  ),
  prototype(covariate = "", fun = ""),
  validity = check_lgpexpr
)

#' An S4 class to represent one formula term
#'
#' @slot factors a list of at most two \linkS4class{lgpexpr}s
lgpterm <- setClass("lgpterm", slots = c(factors = "list"))

#' An S4 class to represent the right-hand side of an lgp formula
#'
#' @slot summands a list of one or more \linkS4class{lgpterm}s
lgprhs <- setClass("lgprhs", slots = c(summands = "list"))

#' An S4 class to represent an lgp formula
#'
#' @slot terms an object of class \linkS4class{lgpsum}
#' @slot response name of the response variable
#' @slot call original formula call
lgpformula <- setClass("lgpformula",
  representation = representation(
    call = "character",
    response = "character",
    terms = "lgprhs"
  ),
  validity = check_lgpformula
)

#' An S4 class to represent variable scaling and its inverse
#'
#' @slot fun function that performs normalization
#' @slot fun_inv inverse function of \code{fun}
lgpscaling <- setClass("lgpscaling",
  representation = representation(
    fun = "function",
    fun_inv = "function"
  ),
  prototype = prototype(
    fun = function(x) {
      x
    },
    fun_inv = function(x) {
      x
    }
  ),
  validity = check_lgpscaling
)

#' An S4 class to represent an lgp model
#'
#' @slot formula An object of class \linkS4class{lgpformula}
#' @slot stan_dat The data to be given as input to \code{rstan::sampling}
#' @slot scalings Variable scaling functions and their inverse operations.
#' Must be a named list with each element is an \linkS4class{lgpscaling}
#' object.
#' @slot info Model info.
lgpmodel <- setClass("lgpmodel",
  representation = representation(
    model_formula = "lgpformula",
    stan_input = "list",
    scalings = "list",
    info = "list"
  )
)

#' An S4 class to represent the output of the \code{lgp_fit} function
#'
#' @slot stan_fit The \code{stanfit} object returned by \code{rstan::sampling}.
#' @slot model The \code{lgpmodel} object returned by \code{lgp_model}.
#' @slot relevances Inferred component relevances.
#' @slot selection Component selection info.
#' @slot pkg_version Package version number.
#' @slot diagnostics  A data frame with columns
#' \code{c("Rhat", "Bulk_ESS", "Tail_ESS")}.
lgpfit <- setClass("lgpfit",
  slots = c(
    stan_fit = "stanfit",
    model = "lgpmodel",
    relevances = "list",
    selection = "list",
    pkg_version = "character",
    diagnostics = "data.frame"
  )
)
