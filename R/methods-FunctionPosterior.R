#' @describeIn FunctionPosterior Print a summary about the object.
setMethod("show", "FunctionPosterior", function(object) {
  comps <- levels(dollar(object@components, "component"))
  paramsets <- levels(dollar(object@total, "paramset"))
  eval_points <- levels(dollar(object@total, "eval_point"))
  num_comps <- length(comps)
  DIMS <- c(length(paramsets), length(eval_points))
  desc <- class_info_fp("FunctionPosterior", num_comps, DIMS)
  cat(desc)
})

#' @describeIn FunctionDraws Print a summary about the object.
setMethod("show", "FunctionDraws", function(object) {
  df <- object@components
  comps <- levels(dollar(df, "component"))
  paramsets <- levels(dollar(df, "paramset"))
  eval_points <- levels(dollar(df, "eval_point"))
  num_comps <- length(comps)
  DIMS <- c(length(paramsets), length(eval_points))
  desc <- class_info_fp("FunctionDraws", num_comps, DIMS)
  cat(desc)
})
