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

#' Parse the covariates and model components from given data and formula
#'
#' @inheritParams parse_response
#' @return parsed input to stan and covariate scaling
parse_data <- function(data, model_formula) {

  # Check that all covariates exist in data
  x_names <- rhs_variables(model_formula@terms)
  x_names <- unique(x_names)
  for (name in x_names) {
    check_in_data(name, data)
  }

  # Create the inputs to Stan
  covariates <- stan_data_covariates(data, x_names)
  components <- stan_data_components(model_formula, covariates)
  to_stan <- c(covariates$to_stan, components$to_stan)

  # Return
  list(
    to_stan = to_stan,
    x_cont_scalings = covariates$x_cont_scalings,
    x_cat_levels = covariates$x_cat_levels
  )
}


#' Create covariate data for Stan input
#'
#' @description Creates the following Stan data input list fields:
#' \itemize{
#'   \item \code{num_cov_cat}
#'   \item \code{num_cov_cont}
#'   \item \code{x_cat}
#'   \item \code{x_cat_num_levels}
#'   \item \code{x_cont}
#'   \item \code{x_cont_mask}
#' }
#' @param data a data frame
#' @param x_names unique covariate names
#' @return a named list with fields
#' \itemize{
#'   \item \code{to_stan}: a list of stan data
#'   \item \code{x_cont_scaling}: normalization function and inverse for each
#'   continuous covariate
#'   \item \code{x_cat_levels}: names of the levels of each categorical
#'   covariate before conversion from factor to numeric
#' }
stan_data_covariates <- function(data, x_names) {
  num_obs <- dim(data)[1]

  x_cont <- list()
  x_cont_mask <- list()
  x_cont_scalings <- list()
  x_cont_names <- c()

  x_cat <- list()
  x_cat_levels <- list()
  x_cat_num_levels <- 0
  x_cat_names <- c()

  num_cat <- 0
  num_cont <- 0

  for (name in x_names) {
    X_RAW <- data[[name]]
    c_x <- class(X_RAW)
    if (c_x == "factor") {

      # A categorical covariate
      num_cat <- num_cat + 1
      n_na <- sum(is.na(X_RAW))
      if (n_na > 0) {
        msg <- paste0(n_na, " missing values for factor '", name, "'!")
        stop(msg)
      }
      x_cat[[num_cat]] <- as.numeric(X_RAW)
      x_cat_num_levels[num_cat] <- length(levels(X_RAW))
      x_cat_levels[[num_cat]] <- levels(X_RAW)
      x_cat_names[num_cat] <- name
    } else if (c_x == "numeric") {

      # A continuous covariate
      num_cont <- num_cont + 1
      is_na <- is.na(X_RAW)
      x_cont_mask[[num_cont]] <- as.numeric(is_na)
      X_NONAN <- X_RAW
      X_NONAN[is_na] <- 0
      normalizer <- create_scaling(X_NONAN, name)
      x_cont_scalings[[num_cont]] <- normalizer
      x_cont[[num_cont]] <- normalizer@fun(X_NONAN)
      x_cont_names[num_cont] <- name
    } else {
      msg <- paste0(
        "Covariate '", name, "' has invalid type '", c_x,
        "'! Must be one of {'factor', 'numeric'}"
      )
      stop(msg)
    }
  }

  # Convert lists to matrices
  x_cat <- list_to_matrix(x_cat, num_obs)
  x_cont <- list_to_matrix(x_cont, num_obs)
  x_cont_mask <- list_to_matrix(x_cont_mask, num_obs)

  # Name lists and matrix rows
  names(x_cont_scalings) <- x_cont_names
  names(x_cat_levels) <- x_cat_names
  rownames(x_cat) <- x_cat_names
  rownames(x_cont) <- x_cont_names
  rownames(x_cont_mask) <- x_cont_names

  # Create Stan data
  to_stan <- list(
    num_cov_cont = num_cont,
    num_cov_cat = num_cat,
    x_cat = x_cat,
    x_cat_num_levels = x_cat_num_levels,
    x_cont = x_cont,
    x_cont_mask = x_cont_mask
  )

  # Return
  list(
    to_stan = to_stan,
    x_cont_scalings = x_cont_scalings,
    x_cat_levels = x_cat_levels
  )
}

#' An lgpterm to numeric representation for Stan
#'
#' @param term an object of class \linkS4class{lgpterm}
#' @param covariates a list returned by \code{\link{stan_data_covariates}}
#' @return a vector of 5 integers
term_to_numeric <- function(term, covariates) {
  facs <- term@factors
  if (length(facs) == 1) {
    out <- rep(1, 9)
  } else {
    out <- rep(2, 9)
  }
  return(out)
}

#' Map a list of terms to their "names"
#'
#' @param rhs an object of class \linkS4class{lgprhs}
#' @return a character vector
term_names <- function(rhs) {
  terms <- rhs@summands
  J <- length(terms)
  names <- c()
  for (j in seq_len(J)) {
    term <- terms[[j]]
    names <- c(names, term_as_character(term, verbose = FALSE))
  }
  return(names)
}

#' Create model component data for Stan input
#'
#' @param model_formula an object of class \linkS4class{lgpformula}
#' @param covariates a list returned by \code{\link{stan_data_covariates}}
#' @return a named list with fields
#' \itemize{
#'   \item \code{to_stan}: a list of stan data
#'   \item TODO: something?
#' }
stan_data_components <- function(model_formula, covariates) {
  terms <- model_formula@terms@summands
  print(covariates) # PTRIN
  J <- length(terms)
  comps <- array(0, dim = c(J, 9))
  for (j in seq_len(J)) {
    comps[j, ] <- term_to_numeric(terms[[j]], covariates)
  }
  colnames(comps) <- c(
    "comp", "ker", " ",
    "heter", "ns", "vm",
    "uncrt", "i_cat", "i_cont"
  )
  rownames(comps) <- term_names(model_formula@terms)
  to_stan <- list(components = as.matrix(comps))

  # Return
  list(
    to_stan = to_stan
  )
}
