#' Operations on formula terms and expressions
#'
#' @description
#' \itemize{
#'   \item Sum two \linkS4class{lgprhs}'s
#'   \item Sum two \linkS4class{lgpterm}'s
#'   \item Sum an \linkS4class{lgprhs} and an \linkS4class{lgpterm}
#'   \item Multiply two \linkS4class{lgpexpr}'s
#' }
#' @param e1 The first sum, term or expression
#' @param e2 The second sum, term or expression
#' @name operations
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
  "*", signature(e1 = "lgpexpr", e2 = "lgpexpr"),
  function(e1, e2) {
    new("lgpterm", factors = list(e1, e2))
  }
)
