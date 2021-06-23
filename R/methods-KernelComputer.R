#' @describeIn KernelComputer Print a summary about the object.
setMethod("show", "KernelComputer", function(object) {
  S <- num_paramsets(object)
  J <- num_components(object)
  P <- num_evalpoints(object)
  desc <- class_info("KernelComputer")
  add_str <- paste0(
    "\n - ", S, " parameter sets",
    "\n - ", J, " components",
    "\n - ", P, " evaluation points "
  )
  add_str <- paste0(add_str)
  add_str <- paste0(
    add_str, "\n - full_covariance = ", object@full_covariance, "\n"
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
  dollar(object@K_input, "n1")
}

#' @describeIn KernelComputer Get number of parameter sets.
setMethod("num_paramsets", "KernelComputer", function(object) {
  dollar(object@input, "num_paramsets")
})

#' @describeIn KernelComputer Get component names.
setMethod("component_names", "KernelComputer", function(object) {
  object@comp_names
})
