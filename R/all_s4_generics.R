#' Summary of model prior
#'
#' @param object an object of class \linkS4class{lgpmodel} or
#' \linkS4class{lgpfit}
#' @param digits number of digits to show for floating point numbers
#' @param ... other arguments
#' @return a data frame
#' @name prior_summary
NULL

#' @rdname prior_summary
setGeneric("prior_summary", function(object, ...) {
  standardGeneric("prior_summary")
})

#' Prior and posterior predictive checks
#'
#' @param object an object of class \linkS4class{lgpmodel}
#' @param ... other arguments
#' @name ppc
#' @aliases prior_predict, posterior_predict
NULL

#' @rdname ppc
setGeneric("posterior_predict", function(object, ...) {
  standardGeneric("posterior_predict")
})

#' @rdname ppc
setGeneric("prior_predict", function(object, ...) {
  standardGeneric("prior_predict")
})
