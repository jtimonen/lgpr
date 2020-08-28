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
