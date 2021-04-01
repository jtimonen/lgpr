#' @describeIn FunctionPosteriors Print a summary about the object.
setMethod("show", "FunctionPosteriors", function(object) {
  comp_names <- component_names(object@model)
  n_eval_points <- nrow(object@x)
  DIMS <- c(object@num_paramsets, n_eval_points)
  desc <- class_info_fp("FunctionPosteriors", comp_names, DIMS)
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

# Class info for show methods of function posterior objects
class_info_fp <- function(class_name, comp_names, D) {
  num_comps <- length(comp_names)
  comp_str <- paste(comp_names, collapse = ", ")
  desc <- class_info(class_name)
  desc <- paste0(desc, "\n - ", num_comps, " components: ", comp_str)
  desc <- paste0(desc, "\n - ", D[1], " parameter set(s)")
  desc <- paste0(desc, "\n - ", D[2], " evaluation points")
  return(desc)
}
