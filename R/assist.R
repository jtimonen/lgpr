#' Validate data array
#'
#' @param data a data frame
#' @param id_name name of the id variable
#' @return a data frame
#' @name valid_data
NULL

#' @export
#' @rdname valid_data
validate_data <- function(data, id_name = "id") {
  id <- dollar(data, id_name)
  check_type(id, "factor")
  D <- dim(data)[2]
  for (j in seq_len(D)) {
    x <- data[, j]
    if (is.factor(x)) {
      nam <- names(data)[j]
      validate_factor(x, id, nam)
    }
  }
  TRUE
}

#' @rdname valid_data
#' @param x a factor
#' @param id the id factor
#' @param name factor name
validate_factor <- function(x, id, name) {
  check_type(x, "factor")
  check_type(id, "factor")
  check_lengths(x, id)
  for (lev in levels(id)) {
    inds <- which(id == lev)
    nu <- length(unique(x[inds]))
    msg <- paste0(
      "measurements corresponding to <id> level ", lev, " do not",
      " all have the same level for factor '", name, "'!"
    )
    if (nu > 1) stop(msg)
  }
  TRUE
}


#' Easily add a categorical covariate to a data frame
#'
#' @export
#' @param data the original data frame
#' @param x A named vector containing the category for each individual.
#' The names should specify the individual id.
#' @param id_var name of the id variable in \code{data}
#' @return A data frame with one column added. The new column will
#' have same name as the variable passed as input \code{x}.
#' @family user assist functions
add_factor <- function(data, x, id_var = "id") {
  check_type(data, "data.frame")
  name <- deparse(substitute(x))
  bad <- name %in% colnames(data)
  if (bad) stop("<data> already contains a variable called '", name, "'!")
  x_id <- as.numeric(names(x))
  L <- length(x_id)
  if (L < 1) stop("<x> must be a named vector!")
  data_id <- dollar(data, id_var)
  uid <- unique(x_id)
  new_factor <- rep(0, length(data_id))
  for (id in uid) {
    i_data <- which(data_id == id)
    i_new <- which(x_id == id)
    new_factor[i_data] <- x[i_new]
  }
  data[[name]] <- as.factor(new_factor)
  return(data)
}

#' Easily add the disease-related age variable to a data frame
#'
#' @export
#' @description Creates the disease-related age covariate vector based on the
#' disease initiation times and adds it to the data frame
#' @param data the original data frame
#' @param t_init A named vector containing the observed initiation or onset
#' time for each individual. The names, i.e. \code{names(t_init)}, should
#' specify the individual id.
#' @param id_var name of the id variable in \code{data}
#' @param time_var name of the time variable in \code{data}
#' @return A data frame with one column added. The new column will
#' be called \code{dis_age}. For controls, its value will be \code{NaN}.
#' @family user assist functions
add_dis_age <- function(data, t_init, id_var = "id", time_var = "age") {
  check_type(data, "data.frame")
  bad <- "dis_age" %in% colnames(data)
  if (bad) stop("<data> already contains a variable called 'dis_age'!")
  x_id <- as.numeric(names(t_init))
  L <- length(x_id)
  if (L < 1) stop("<t_init> must be a named vector!")
  data_id <- dollar(data, id_var)
  data_age <- dollar(data, time_var)
  uid <- unique(x_id)
  dage <- rep(NaN, length(data_id))
  for (id in uid) {
    i_data <- which(data_id == id)
    i_new <- which(x_id == id)
    dage[i_data] <- data_age[i_data] - t_init[i_new]
  }
  data$dis_age <- dage
  return(data)
}

#' Create the GP mean vector
#'
#' @export
#' @description Creates the \code{c_hat} input for \code{lgp},
#' so that it accounts for normalization between data points in the
#' Poisson or NB observation model
#' @param y response variable, vector of length \code{n}
#' @param norm_factors normalization factors, vector of length \code{n}
#' @return a vector of length \code{n}, which can be used as
#' the \code{c_hat} input to the \code{lgp} function
#' @family user assist functions
adjusted_c_hat <- function(y, norm_factors) {
  L1 <- length(norm_factors)
  L2 <- length(y)
  if (L1 != L2) stop("inputs must have same length!")
  if (sum(y < 0) > 0) stop("y cannot have negative values!")
  if (sum(round(y) != y) > 0) stop("y must have only integer values!")
  if (sum(norm_factors <= 0) > 0) stop("norm_factors must be all positive!")
  c_hat <- log(mean(y))
  c_hat <- rep(c_hat, length(y))
  c_hat <- c_hat + log(norm_factors)
  return(c_hat)
}


#' Split data into training and test sets
#'
#' @name split
#' @param data a data frame
#' @param var_name name of a factor in the data
#' @description
#' \itemize{
#'   \item \code{split_by_factor} splits according to given factor
#'   \item \code{split_within_factor} splits according to given
#'   data point indices within the same level of a factor
#'   \item \code{split_within_factor_random} selects k points
#'   from each level of a factor uniformly at random as test data
#'   \item \code{split_random} splits uniformly at random
#'   \item \code{split_data} splits according to given data rows
#' }
#' @return a named list with names \code{train}, \code{test}, \code{i_train}
#' and \code{i_test}
#' @family user assist functions
NULL

#' @export
#' @rdname split
#' @param test the levels of the factor that will be used as test data
split_by_factor <- function(data, test, var_name = "id") {
  check_type(data, "data.frame")
  fac <- dollar(data, var_name)
  check_type(fac, "factor")
  i_test <- which(fac %in% test)
  split_data(data, i_test)
}

#' @export
#' @rdname split
#' @param idx_test indices point indices with the factor
split_within_factor <- function(data, idx_test, var_name = "id") {
  check_type(data, "data.frame")
  id <- dollar(data, var_name)
  check_type(id, "factor")
  uid <- unique(id)
  i_test <- c()
  for (i in uid) {
    inds <- which(id == i)
    i_test <- c(i_test, inds[idx_test])
  }
  split_data(data, i_test)
}

#' @export
#' @rdname split
#' @param k_test desired number of test data points per each level of the factor
split_within_factor_random <- function(data, k_test = 1, var_name = "id") {
  check_type(data, "data.frame")
  id <- dollar(data, var_name)
  check_type(id, "factor")
  uid <- unique(id)
  i_test <- c()
  for (i in uid) {
    inds <- which(id == i)
    L <- length(inds)
    i_sel <- sample.int(L, k_test)
    i_test <- c(i_test, inds[i_sel])
  }
  split_data(data, i_test)
}

#' @export
#' @rdname split
#' @param p_test desired proportion of test data
#' @param n_test desired number of test data points (if NULL, \code{p_test}
#' is used to compute this)
split_random <- function(data, p_test = 0.2, n_test = NULL) {
  check_type(data, "data.frame")
  n_total <- dim(data)[1]
  if (is.null(n_test)) {
    n_test <- round(p_test * n_total)
  }
  i_test <- sample.int(n_total, size = n_test)
  split_data(data, i_test)
}

#' @export
#' @rdname split
#' @param i_test test data row indices
#' @param sort_ids should the test indices be sorted into increasing order
split_data <- function(data, i_test, sort_ids = TRUE) {
  check_type(data, "data.frame")
  i_test <- if (sort_ids) sort(i_test, decreasing = FALSE) else sort_ids
  n_total <- dim(data)[1]
  i_train <- setdiff(seq_len(n_total), i_test)

  # Return
  list(
    train = data[i_train, ],
    test = data[i_test, ],
    i_train = i_train,
    i_test = i_test
  )
}
