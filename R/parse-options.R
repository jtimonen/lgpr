#' Parse the given modeling options
#'
#' @param options A named list with the following possible fields:
#' \itemize{
#'   \item \code{sample_f} Determines if the function values are be sampled
#'   (must be \code{TRUE} if likelihood is not \code{"gaussian"}).
#'   \item \code{skip_generated} If this is true, the generated quantities
#'   block of Stan is skipped.
#'   \item \code{delta} Amount of added jitter to ensure positive definite
#'   covariance matrices.
#'   \item \code{verbose} Should more verbose output be printed?
#' }
#' @return a named list of parsed options
parse_options <- function(options = NULL) {
  input <- options
  
  # Set defaults
  opts <- list(
    verbose = FALSE,
    skip_generated = FALSE,
    sample_f = FALSE,
    delta = 1e-8
  )
  
  # Replace defaults if found from input
  for (opt_name in names(opts)) {
    if (opt_name %in% names(input)) {
      opts[[opt_name]] <- input[[opt_name]]
    }
  }
  
  # Format for Stan input
  list(
    verbose = as.numeric(opts$verbose),
    skip_generated = as.numeric(opts$skip_generated),
    sample_f = as.numeric(opts$sample_f),
    delta = opts$delta
  )
}

#' Parse the given modeling options for the disease component
#'
#' @param disease_options A named list with the following possible fields:
#' \itemize{
#'   \item \code{heterogeneous} Is the disease effect magnitude assumed to be
#'   heterogeneous over the diseased individuals?
#'   \item \code{uncertain} Do we wish to model uncertainty in the disease
#' effect time?
#'   \item \code{vm_params} Parameters of the variance masking function. Use
#'   \code{NA} to avoid using variance masking at all.
#' }
#' @return a named list of parsed options
parse_disease_options <- function(disease_options = NULL) {
  input <- disease_options
  
  # Set defaults
  opts <- list(
    uncertain = FALSE,
    heterogeneous = FALSE,
    vm_params = c(0.95, 1.0)
  )
  
  # Replece defaults if found from input
  for (opt_name in names(opts)) {
    if (opt_name %in% names(input)) {
      opts[[opt_name]] <- input[[opt_name]]
    }
  }
  
  # Format for Stan input
  if (is.na(opts$vm_params) || is.null(opts$vm_params)) {
    vm_params <- array(0, dim = c(0, 2))
  } else {
    vm_params <- array(opts$vm_params, dim = c(1, 2))
  }
  list(
    uncertain = as.numeric(opts$uncertain),
    heterogeneous = as.numeric(opts$heterogeneous),
    vm_params = vm_params
  )
}
