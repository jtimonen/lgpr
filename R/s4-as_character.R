#' Character representations of different S4 objects
#'
#' @param x an object of some S4 class
#' @return a character representation of the object
#' @name as_character
NULL

#' @rdname as_character
setMethod("as.character", "lgpexpr", function(x) {
  paste0(x@fun, "(", x@covariate, ")")
})

#' @rdname as_character
setMethod("as.character", "lgpterm", function(x) {
  facs <- x@factors
  desc <- sapply(facs, as.character)
  desc <- paste(desc, collapse = ", ")
  return(desc)
})

#' @rdname as_character
setMethod("as.character", "lgprhs", function(x) {
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
})

#' @rdname as_character
setMethod("as.character", "lgpformula", function(x) {
  return(x@call)
})

#' @rdname as_character
setMethod("as.character", "lgpsim", function(x) {
  n <- dim(x@data)[1]
  desc <- paste0("A simulated dataset with ", n, " data points.\n")
  return(desc)
})

#' @rdname as_character
setMethod("as.character", "lgpmodel", function(x) {
  str0 <- as.character(x@model_formula)
  str1 <- get_covariate_names(x, type = "continuous")
  str2 <- get_covariate_names(x, type = "categorical")
  desc <- paste0("Formula: ", str0)
  desc <- paste0(desc, "\nContinuous: {", str1, "}")
  desc <- paste0(desc, "\nCategorical: {", str2, "}")
  return(desc)
})
