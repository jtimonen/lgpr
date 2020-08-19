#' Easily add a categorical covariate to a data frame
#'
#' @export
#' @param data the original data frame
#' @param x A named vector containing the category for each individual.
#' The names should specify the individual id.
#' @param id_var name of the id variable in \code{data}
#' @return A data frame with one column added. The new column will
#' have same name as the variable passed as input \code{x}.
add_categorical_covariate <- function(data, x, id_var = "id") {
  name <- deparse(substitute(x))
  if (name %in% colnames(data)) {
    stop("The data frame already contains a variable called '", name, "'!")
  }
  x_id <- as.numeric(names(x))
  data_id <- data[[id_var]]
  uid <- unique(x_id)
  xx <- rep(0, length(data_id))
  for (id in uid) {
    i_data <- which(data_id == id)
    i_new <- which(x_id == id)
    xx[i_data] <- x[i_new]
  }
  data[[name]] <- xx
  return(data)
}


#' Create the disease-related age covariate vector based on the
#' disease initiation times and add it to the data frame
#'
#' @export
#' @param data the original data frame
#' @param t_init A named vector containing the observed initiation or onset
#' time for each individual. The names, i.e. \code{names(t_init)}, should
#' specify the individual id.
#' @param id_var name of the id variable in \code{data}
#' @param time_var name of the time variable in \code{data}
#' @return A data frame with one column added. The new column will
#' be called \code{'diseaseAge'}. For controls, the value of diseaseAge
#' will be set to NaN.
add_disease_ages <- function(data, t_init, id_var = "id", time_var = "age") {
  if ("diseaseAge" %in% colnames(data)) {
    stop("The data frame already contains a variable called 'diseaseAge'!")
  }
  x_id <- as.numeric(names(t_init))
  data_id <- data[[id_var]]
  data_age <- data[[time_var]]
  uid <- unique(x_id)
  dage <- rep(NaN, length(data_id))
  for (id in uid) {
    i_data <- which(data_id == id)
    i_new <- which(x_id == id)
    dage[i_data] <- data_age[i_data] - t_init[i_new]
  }
  data$diseaseAge <- dage
  return(data)
}


#' Create the GP mean input for \code{lgp}, so that it accounts
#' for normalization between data points in the Poisson or NB
#' observation model
#'
#' @export
#' @param y response variable, vector of length \code{n}
#' @param norm_factors normalization factors, vector of length \code{n}
#' @return a vector of length \code{n}, which can be used as
#' the \code{c_hat} input to the \code{lgp} function
adjusted_c_hat <- function(y, norm_factors) {
  if (length(norm_factors) != length(y)) {
    stop("inputs must have same length!")
  }
  if (sum(y < 0) > 0) {
    stop("y cannot have negative values!")
  }
  if (sum(round(y) != y) > 0) {
    stop("y must have only integer values!")
  }
  if (sum(norm_factors <= 0) > 0) {
    stop("norm_factors must be all positive!")
  }

  c_hat <- log(mean(y))
  c_hat <- rep(c_hat, length(y))
  c_hat <- c_hat + log(norm_factors)
  return(c_hat)
}
