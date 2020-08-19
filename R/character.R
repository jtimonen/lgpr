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


#' Is an object printable
#'
#' @param object an object
is_printable <- function(object) {
  d <- dim(object)
  if (is.null(d)) {
    return(TRUE)
  } else if (prod(d) == 0) {
    return(FALSE)
  }
  TRUE
}

#' Print a list in a more compact format
#'
#' @param input a named list
print_list <- function(input) {
  nam <- names(input)
  printed <- c()
  skipped <- c()
  for (name in nam) {
    f <- input[[name]]
    if (is_printable(f)) {
      printed <- c(printed, name)
    } else {
      skipped <- c(skipped, name)
    }
  }

  print(input[printed])
  str <- paste(skipped, collapse = ", ")
  msg <- paste0(
    "Did not print fields with at least one zero dimension:\n    ",
    str, "\n"
  )
  cat(msg)
  invisible(input)
}

#' Print the Stan input of an lgpmodel
#'
#' @param model an object of class \linkS4class{lgpmodel}
print_stan_input <- function(model) {
  print_list(model@stan_input)
  invisible(model)
}
