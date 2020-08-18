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


#' Check that data contains a variable with a certain name
#'
#' @param var_name the variable to be searched for
#' @param data an object of class \code{data.frame}
#' @return \code{TRUE} if the variable is found
check_in_data <- function(var_name, data) {
  d_names <- colnames(data)
  ok <- (var_name %in% d_names)
  if (!ok) {
    str <- paste(d_names, collapse = ", ")
    msg <- paste0(
      "The variable '", var_name, "' not found in <data>! ",
      " Found data columns = {", str, "}."
    )
    stop(msg)
  }
  return(TRUE)
}

#' Create the function that does a standardizing transform and its inverse
#'
#' @param y variable measurements
#' @param name variable name
#' @return an object of class \linkS4class{lgpscaling}
create_scaling <- function(y, name) {
  if (length(y) < 2) {
    stop("length of <y> must be at least 2!")
  }
  m <- mean(y)
  std <- stats::sd(y)
  if (std == 0) {
    stop("the varible measurements have zero variance!")
  }
  fun <- function(x) {
    (x - m) / std
  }
  fun_inv <- function(x) {
    x * std + m
  }
  new("lgpscaling", fun = fun, fun_inv = fun_inv, var_name = name)
}
