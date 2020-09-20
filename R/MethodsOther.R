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
setMethod("as.character", "lgpformula", function(x) {
  return(x@call)
})

#' Operations on formula terms and expressions
#'
#' @description
#' \itemize{
#'   \item Sum two \linkS4class{lgprhs}'s to yield an \linkS4class{lgprhs}
#'   \item Sum two \linkS4class{lgpterm}'s to yield an \linkS4class{lgprhs}
#'   \item Sum an \linkS4class{lgprhs} and an \linkS4class{lgpterm}
#'   to yield an \linkS4class{lgprhs}
#'   \item Multiply two \linkS4class{lgpterm}'s to yield
#'   an \linkS4class{lgpterm}
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
  "*", signature(e1 = "lgpterm", e2 = "lgpterm"),
  function(e1, e2) {
    new("lgpterm", factors = c(e1@factors, e2@factors))
  }
)

#' Printing S4 object info using the show generic
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
  invisible(object)
})

#' @rdname show
setMethod("show", "lgpsim", function(object) {
  msg <- class_info("lgpsim")
  cat(msg)
  invisible(object)
})

#' @rdname show
setMethod("show", "lgpmodel", function(object) {
  msg <- class_info("lgpmodel")
  cat(msg)
  cat("\n")
  model_summary(object)
})

#' @rdname show
setMethod("show", "lgpfit", function(object) {
  msg <- class_info("lgpfit")
  cat(msg)
  cat("\n")
  fit_summary(object)
})

#' Class info for show methods
#'
#' @param class_name class name
#' @return a string
class_info <- function(class_name) {
  str <- paste0(
    "An object of class ", class_name, ". See ?",
    class_name, " for more info."
  )
  return(str)
}


#' Get simulated components from an lgpsim object
#'
#' @export
#' @param simdata an object of class \linkS4class{lgpsim}
#' @return a data frame
get_sim_components <- function(simdata) {
  df <- simdata@components
  nams <- colnames(df)
  idx <- which(nams == "f")
  df <- df[, 1:(idx - 1)]
  nams <- colnames(df)
  prettify <- function(str, i) {
    str <- gsub(str, pattern = ".", replacement = ", ", fixed = TRUE)
    str <- paste0("f[", i, "](", str, ")")
    return(str)
  }
  J <- length(nams)
  for (j in seq_len(J)) {
    nams[j] <- prettify(nams[j], j)
  }
  colnames(df) <- nams
  return(df)
}
