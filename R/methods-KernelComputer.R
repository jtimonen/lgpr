#' @describeIn KernelComputer Print a summary about the object.
setMethod("show", "KernelComputer", function(object) {
  desc <- class_info("KernelComputer")
  add_str <- paste0(
    " Three same matrices: ",
    three_matrices_are_same(object), "."
  )
  desc <- paste0(desc, add_str)
  cat(desc)
})


#' @describeIn KernelComputer Get number of components.
setMethod("num_components", "KernelComputer", function(object) {
  K_const <- dollar(object@K_input, "K_const")
  length(K_const)
})

#' @describeIn KernelComputer Get number of evaluation points.
setMethod("num_evalpoints", "KernelComputer", function(object) {
  dollar(object@Ks_input, "n1")
})

# Get number of observations from KernelComputer object
get_num_obs_kc <- function(object) {
  dollar(object@Ks_input, "n2")
}

#' @describeIn KernelComputer Get number of parameter sets.
setMethod("num_paramsets", "KernelComputer", function(object) {
  dollar(object@input, "num_paramsets")
})

#' @describeIn KernelComputer Get component names.
setMethod("component_names", "KernelComputer", function(object) {
  object@comp_names
})

# Determine if K = Ks = Kss for a KernelComputer
no_separate_output_points <- function(kc) {
  L <- length(kc@Ks_input)
  return(L == 0)
}
