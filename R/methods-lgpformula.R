#' Character representations of different formula objects
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
setMethod("as.character", "lgpformula", function(x) {
  return(x@call)
})

#' Operations on formula terms and expressions
#'
#' @param e1 The first sum, term or expression
#' @param e2 The second sum, term or expression
#' @name operations
#' @return
#' The behaviour and return type depend on the types of \code{e1} and \code{e2}.
#' You can
#' \itemize{
#'   \item sum (\code{+}) two \linkS4class{lgprhs}'s to yield an
#'   \linkS4class{lgprhs}
#'   \item sum (\code{+}) two \linkS4class{lgpterm}'s to yield an
#'   \linkS4class{lgprhs}
#'   \item sum (\code{+}) an \linkS4class{lgprhs} and an \linkS4class{lgpterm}
#'   to yield an \linkS4class{lgprhs}
#'   \item multiply (\code{*}) two \linkS4class{lgpterm}'s to yield
#'   an \linkS4class{lgpterm}
#' }
NULL

#' @rdname operations
setMethod(
  "+", signature(e1 = "lgprhs", e2 = "lgprhs"),
  function(e1, e2) {
    new("lgprhs", summands = c(e1@summands, e2@summands))
  }
)

#' @rdname operations
setMethod(
  "+", signature(e1 = "lgpterm", e2 = "lgpterm"),
  function(e1, e2) {
    new("lgprhs", summands = list(e1, e2))
  }
)

#' @rdname operations
setMethod(
  "+", signature(e1 = "lgprhs", e2 = "lgpterm"),
  function(e1, e2) {
    e1 + new("lgprhs", summands = list(e2))
  }
)

#' @rdname operations
setMethod(
  "*", signature(e1 = "lgpterm", e2 = "lgpterm"),
  function(e1, e2) {
    new("lgpterm", factors = c(e1@factors, e2@factors))
  }
)

#' Printing formula object info using the show generic
#'
#' @name show
#' @param object an object of some S4 class
#' @return the object invisibly
NULL

#' @rdname show
setMethod("show", "lgpformula", function(object) {
  cat(as.character(object))
  invisible(object)
})

#' @rdname show
setMethod("show", "lgprhs", function(object) {
  s <- object@summands
  L <- length(s)
  cat("An object of class lgprhs.\n\n")
  for (j in seq_len(L)) {
    t <- as.character(s[[j]])
    r <- paste0("Term ", j, ": ", t, "\n")
    cat(r)
  }
  invisible(object)
})

#' @rdname show
setMethod("show", "lgpterm", function(object) {
  cat("lgpterm: ")
  cat(as.character(object))
  cat("\n")
  invisible(object)
})
