#' Character representations of different S4 objects
#'
#' @param object an object of some S4 class
#' @return a character representation of the object
#' @name as_character
NULL

#' @rdname as_character
setMethod(
  f = "as.character", signature = "lgpexpr",
  definition = function(x) {
    paste0(x@fun, "(", x@covariate, ")")
  }
)

#' @rdname as_character
setMethod(
  f = "as.character", signature = "lgpterm",
  definition = function(x) {
    facs <- x@factors
    L <- length(facs)
    c1 <- as.character(facs[[1]])
    if (L == 1) {
      desc <- paste0("(1st order):   ", c1)
    } else {
      c2 <- as.character(facs[[2]])
      desc <- paste0("(interaction): ", c1, "*", c2)
    }
    return(desc)
  }
)

#' @rdname as_character
setMethod(
  f = "as.character", signature = "lgprhs",
  definition = function(x) {
    s <- x@summands
    L <- length(s)
    desc <- ""
    for (j in seq_len(L)) {
      desc <- paste0(desc, "Term ", j, " ", as.character(s[[j]]), "\n")
    }
    return(desc)
  }
)

#' @rdname as_character
setMethod(
  f = "as.character", signature = "lgpformula",
  definition = function(x) {
    desc <- paste0("Formula:              ", x@call, "\n")
    desc <- paste0(desc, "Response variable:    ", x@y_name, "\n")
    desc <- paste0(desc, as.character(x@terms))
    return(desc)
  }
)

#' @rdname as_character
setMethod(
  f = "as.character", signature = "lgpmodel",
  definition = function(x) {
    str1 <- covariate_names(x, type = "continuous")
    str2 <- covariate_names(x, type = "categorical")
    desc <- "An object of class lgpmodel.\n"
    desc <- paste0(desc, "\nFormula: ", x@model_formula@call)
    desc <- paste0(desc, "\nContinuous: {", str1, "}")
    desc <- paste0(desc, "\nCategorical: {", str2, "}")
    return(desc)
  }
)
