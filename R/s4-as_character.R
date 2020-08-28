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
  desc <- paste(desc, collapse = "*")
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
