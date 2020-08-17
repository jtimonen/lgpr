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
    is_verbose = as.numeric(opts$verbose),
    is_generated_skipped = as.numeric(opts$skip_generated),
    is_f_sampled = as.numeric(opts$sample_f),
    delta = opts$delta
  )
}
