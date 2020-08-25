#' Printing S4 object info using the show generic
#'
#' @param object an object of some S4 class
#' @return the object invisibly
#' @name show

#' @rdname show
setMethod(
  f = "show", signature = "lgpformula",
  definition = function(object) {
    cat(as.character(object))
    invisible(object)
  }
)

#' @rdname show
setMethod(
  f = "show",
  signature = "lgpfit",
  definition = function(object) {
    call <- object@model@model_formula@call
    info <- paste0("Model formula: ", call, "\n\n")
    cat(info)
    print(object@stan_fit, pars = "f_post", include = FALSE)
  }
)

#' @rdname show
setMethod(
  f = "show", signature = "lgpmodel",
  definition = function(object) {
    cat(as.character(object))
    invisible(object)
  }
)

#' @rdname show
setMethod(
  f = "show", signature = "lgpsim",
  definition = function(object) {
    cat(as.character(object))
    invisible(object)
  }
)
