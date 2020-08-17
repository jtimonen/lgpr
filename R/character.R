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

#' Character format of an lgpterm
#'
#' @param x an object of class \linkS4class{lgpterm}
#' @param verbose should the format be more verbose
#' @return a string
term_as_character <- function(x, verbose = TRUE) {
  facs <- x@factors
  desc <- sapply(facs, as.character)
  desc <- paste(desc, collapse = ", ")
  return(desc)
}

#' @rdname as_character
setMethod(
  f = "as.character", signature = "lgpterm",
  definition = function(x) {
    term_as_character(x, verbose = TRUE)
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
      desc <- paste0(
        desc, "Term ", j, " expressions:   ",
        as.character(s[[j]]), "\n"
      )
    }
    return(desc)
  }
)


#' @rdname as_character
setMethod(
  f = "as.character", signature = "lgpformula",
  definition = function(x) {
    desc <- "An object of class 'lgpformula'.\n\n"
    desc <- paste0(desc, "Call:                 ", x@call, "\n")
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
