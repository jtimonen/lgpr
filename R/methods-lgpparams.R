#' @describeIn lgpparams Get number of parameter sets.
setMethod("num_paramsets", "lgpparams", function(object) {
  dollar(object@input, "num_paramsets")
})
