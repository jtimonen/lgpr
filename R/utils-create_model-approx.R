#' Parse the given approximation options
#'
#' @param approx A named list with the following possible fields:
#' \itemize{
#'   \item \code{num_bf} Number of basis functions (0 = no approximation).
#'   \item \code{scale_bf} Scale of the domain to be used in basis
#'   function approximation. Has no effect if \code{num_bf = 0}.
#' }
#' If \code{approx} is \code{NULL}, default options are used. The defaults
#' are equivalent to \code{options = list(num_bf = 0, scale_bf = 1.5)}.
#' @return a named list of parsed options
create_model.approx_options <- function(approx) {

  # Default options
  input <- approx
  opts <- list(num_bf = 0, scale_bf = 1.5)

  # Replace defaults if found from input
  for (opt_name in names(opts)) {
    if (opt_name %in% names(input)) {
      opts[[opt_name]] <- input[[opt_name]]
    }
  }

  # Validate
  num_bf <- dollar(opts, "num_bf")
  scale_bf <- dollar(opts, "scale_bf")
  check_non_negative_all(num_bf)
  check_non_negative_all(scale_bf)
  return(opts)
}
