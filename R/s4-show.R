#' Printing S4 object info using the show generic
#'
#' @description All of these call the corresponding \code{as.character}
#' method and then \code{\link{show_default}} which prints the description
#' and returns the original object invisibly.
#'
#' @name show
#' @param object an object of some S4 class
#' @return the object invisibly
#' @seealso See also \code{\link{show_default}} and \code{\link{as_character}}.
NULL

#' @rdname show
setMethod("show", "lgpformula", function(object) show_default(object))

#' @rdname show
setMethod("show", "lgpmodel", function(object) show_default(object))

#' @rdname show
setMethod("show", "lgpfit", function(object) show_default(object))

#' @rdname show
setMethod("show", "lgpsim", function(object) show_default(object))

#' Print the character representation of an S4 object
#'
#' @param object an object of some S4 class
#' @return the object invisibly
show_default <- function(object) {
  cat(as.character(object))
  invisible(object)
}
