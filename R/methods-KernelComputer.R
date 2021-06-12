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

# Determine if K = Ks = Kss for a KernelComputer
three_matrices_are_same <- function(kc) {
  is.null(dollar(kc@init, "Ks_init"))
}
