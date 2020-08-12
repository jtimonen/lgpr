#' Parse the covariates and model components from given data and formula
#'
#' @inheritParams parse_response
#' @param id_variable Name of the unique subject identifier variable
#' (default = \code{"id"}).
#' @return parsed input to stan and covariate scaling
parse_data <- function(data, model_formula, id_variable) {

  # Check that all covariates exist in data
  x_names <- rhs_variables(model_formula@terms)
  x_names <- unique(x_names)
  for (name in x_names) {
    check_in_data(name, data)
  }

  # Check that the id variable is in data
  check_in_data(id_variable, data)

  # Create x_cat, x_cont, and x_mask
  parsed <- stan_data_covariates(data, x_names)
  scaling <- parsed$x_cont_scaling

  # Create the list that will go as input to stan
  to_stan <- list()
  to_stan <- c(to_stan, parsed$to_stan)

  # Return
  list(to_stan = to_stan, scaling = scaling)
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
#'   \item \code{x_cont_normalized}
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
  x_cont_scaling <- list()
  x_cont_normalized <- list()

  x_cat <- list()
  x_cat_levels <- list()
  x_cat_num_levels <- 0

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
    } else if (c_x == "numeric") {

      # A continuous covariate
      num_cont <- num_cont + 1
      is_na <- is.na(X_RAW)
      x_cont_mask[[num_cont]] <- as.numeric(is_na)
      X_NONAN <- X_RAW
      X_NONAN[is_na] <- 0
      x_cont[[num_cont]] <- X_NONAN
      normalizer <- create_scaling(X_NONAN, name)
      x_cont_scaling[[num_cont]] <- normalizer
      x_cont_normalized[[num_cont]] <- normalizer@fun(X_NONAN)
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
  x_cont_normalized <- list_to_matrix(x_cont_normalized, num_obs)

  # Return
  to_stan <- list(
    num_cov_cont = num_cont,
    num_cov_cat = num_cat,
    x_cat = x_cat,
    x_cat_num_levels = x_cat_num_levels,
    x_cont = x_cont,
    x_cont_mask = x_cont_mask,
    x_cont_normalized = x_cont_normalized
  )

  list(
    to_stan = to_stan,
    x_cont_scaling = x_cont_scaling,
    x_cat_levels = x_cat_levels
  )
}
