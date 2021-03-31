#' @describeIn FunctionPosterior Print a summary about the object.
setMethod("show", "FunctionPosterior", function(object) {
  comp_names <- component_names(object@model)
  n_eval_points <- nrow(object@x)
  DIMS <- c(object@num_paramsets, n_eval_points)
  desc <- class_info_fp("FunctionPosterior", comp_names, DIMS)
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
