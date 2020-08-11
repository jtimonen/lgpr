#' Count numbers of different categories for each categorical variable
#'
#' @param X the design matrix
#' @param D a vector of length 6
#' @return a numeric vector
set_num_levels <- function(X, D) {
  N_cat <- rep(0, 1 + D[5] + D[6])
  j0 <- 2 + D[3] + D[4]
  x1 <- X[, 1]
  N_cat[1] <- length(unique(x1)) # should equal to number of individuals
  if (D[5] > 0) {
    for (j in 1:D[5]) {
      xj <- X[, j0 + j]
      N_cat[1 + j] <- length(unique(xj))
    }
  }
  if (D[6] > 0) {
    for (j in 1:D[6]) {
      xj <- X[, j0 + D[5] + j]
      N_cat[1 + D[5] + j] <- length(unique(xj))
    }
  }
  return(as.array(N_cat))
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
#' @param y response variable measurements
#' @return an object of class \linkS4class{lgpscaling}
create_scaling <- function(y) {
  if (length(y) < 2) {
    stop("length of <y> must be at least 2!")
  }
  y_m <- mean(y)
  y_std <- stats::sd(y)
  if (y_std == 0) {
    stop("<y> has zero variance!")
  }
  fun <- function(y) {
    (y - y_m) / y_std
  }
  fun_inv <- function(y) {
    y * y_std + y_m
  }
  new("lgpscaling", fun = fun, fun_inv = fun_inv)
}

#' Check that the response is numeric and compatibile with observation model
#'
#' @param y the response variable measurements
#' @param obs_model observation model (integer encoding)
check_response <- function(y, obs_model) {

  # Check that y is numeric
  c_y <- class(y)
  if (c_y != "numeric") {
    stop("the response variable must be numeric! found = ", c_y)
  }

  # Check compatibility with observation model
  if (obs_model != 1) {
    diff <- max(abs(y - round(y)))
    if (diff > 0) {
      stop(
        "the response variable should contain only integers ",
        "with this observation model"
      )
    }
    if (any(y < 0)) {
      stop(
        "the response variable measurements cannot be negative ",
        "with this observation model!"
      )
    }
  }
  return(TRUE)
}

#' Parse the response variable from given data and formula
#'
#' @param data A \code{data.frame} where each column corresponds to one
#' variable, and each row is one observation. Continuous covariates and the
#' response variable must have type \code{"numeric"} and categorical covariates
#' must have type \code{"factor"}. Missing values should be indicated with
#' \code{NaN} or \code{NA}. The response variable cannot contain missing
#' values.
#' @param likelihood Determines the observation model. Must be either
#' \code{"gaussian"} (default), \code{"poisson"}, \code{"nb"} (negative
#' binomial) or \code{"binomial"}. To use Bernoulli likelihood, use
#' \code{likelihood="binomial"} and set \code{num_trials} as a vector of ones.
#' @param model_formula An object of class \code{lgpformula}.
#' @return a named list of parsed options
parse_response <- function(data, likelihood, model_formula) {
  obs_model <- likelihood_as_int(likelihood)

  # Check that data is a data.frame and contains the response variable
  y_name <- model_formula@response
  c_data <- class(data)
  if (c_data != "data.frame") {
    stop("<data> must be a data.frame! found = ", c_data)
  }
  check_in_data(y_name, data)
  y <- data[[y_name]]

  # Check that the response is numeric and compatible with observation model
  check_response(y, obs_model)

  # Create y_cont and y_disc inputs for Stan
  num_obs <- length(y)
  if (obs_model != 1) {
    y_disc <- array(y, dim = c(1, num_obs))
    y_cont <- array(y, dim = c(0, num_obs))
    y_scaling <- NA
  } else {
    y_scaling <- create_scaling(y) # create scaling and inverse
    y <- y_scaling@fun(y) # standardize the response
    y_disc <- array(y, dim = c(0, num_obs))
    y_cont <- array(y, dim = c(1, num_obs))
  }
  to_stan <- list(num_obs = num_obs, y_disc = y_disc, y_cont = y_cont)

  # Return also the scaling and its inverse
  list(
    y_to_stan = to_stan,
    y_scaling = y_scaling
  )
}


#' Parse the covariates and model components from given data and formula
#'
#' @inheritParams parse_response
#' @param id_variable Name of the unique subject identifier variable
#' (default = \code{"id"}).
#' @return parsed input to stan and covariate scaling
parse_covariates <- function(data, model_formula, id_variable) {

  # Check that all covariates exist in data
  x_names <- rhs_variables(model_formula@terms)
  for (name in x_names) {
    check_in_data(name, data)
  }

  # Check that the id variable is in data
  check_in_data(id_variable, data)

  # TODO: get covariate types, mask NaNs, check NaNs,
  types <- c()
  for (name in x_names) {
    x <- data[[name]]
    c_x <- class(x)
    c_type <- if (c_x == "factor") "discrete" else "continuous"
    types <- c(types, c_type)
  }

  # Create the list that will go as input to stan
  to_stan <- list(
    num_subjects = 1,
    num_cases = 0,
    num_cov_cont = sum(types == "continuous"),
    num_cov_disc = sum(types == "discrete"),
    num_levels = 1,
    x_disc = 1,
    x_cont = 1
  )

  # Return
  x_scaling <- NA
  list(x_to_stan = to_stan, x_scaling = x_scaling)
}
