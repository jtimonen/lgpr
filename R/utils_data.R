#' Easily add a categorical covariate to a data frame
#'
#' @export
#' @param data the original data frame
#' @param x A named vector containing the category for each individual.
#' The names should specify the individual id.
#' @param id_var name of the id variable in \code{data}
#' @return A data frame with one column added. The new column will
#' have same name as the variable passed as input \code{x}.
#' @family data utilities
add_factor <- function(data, x, id_var = "id") {
  check_type(data, "data.frame")
  name <- deparse(substitute(x))
  bad <- name %in% colnames(data)
  if (bad) stop("<data> already contains a variable called '", name, "'!")
  check_named(x)
  x_id <- as.numeric(names(x))
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
#' @family data utilities
add_dis_age <- function(data, t_init, id_var = "id", time_var = "age") {
  check_type(data, "data.frame")
  bad <- "dis_age" %in% colnames(data)
  if (bad) stop("<data> already contains a variable called 'dis_age'!")
  check_named(t_init)
  x_id <- as.numeric(names(t_init))
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
#' @family data utilities
NULL

#' @rdname split
#' @param test the levels of the factor that will be used as test data
split_by_factor <- function(data, test, var_name = "id") {
  check_type(data, "data.frame")
  fac <- dollar(data, var_name)
  check_type(fac, "factor")
  i_test <- which(fac %in% test)
  split_data(data, i_test)
}

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

#' @rdname split
#' @param k_test desired number of test data points per each level of the
#' factor
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

#' Create prediction points
#'
#' @description Replaces a continuous variable \code{x} in the data frame, and
#' possibly another continuous variable \code{x_ns} derived from it, with new
#' values, for each level of a grouping factor (usually id)
#' @export
#' @param data a data frame
#' @param group_by name of the grouping varible, must be a factor
#' in \code{data}
#' @param x of the variable along which to extend,
#' must be a numeric in \code{data}
#' @param x_ns of a nonstationary variable derived from \code{x},
#' must be a numeric in \code{data}
#' @param x_values the values of \code{x} to set for each individual
#' @return a data frame containing the following columns
#' \itemize{
#'  \item all factors in the original \code{data}
#'  \item \code{x}
#'  \item \code{x_ns} (unless it is NULL)
#' }
#'
#' @family data utilities
new_x <- function(data, x_values, group_by = "id", x = "age", x_ns = NULL) {
  check_type(data, "data.frame")
  check_not_null(group_by)
  check_not_null(x)
  check_not_null(x_values)
  check_in_data(group_by, data)
  check_in_data(x, data)
  df <- pick_one_row_each(data, group_by)
  k <- length(x_values)
  col_names <- if (is.null(x_ns)) x else c(x, x_ns)
  df <- select_factors_and(df, col_names)
  N <- nrow(df)
  inds <- rep(seq_len(N), each = k)
  df <- df[inds, ]
  x_val_rep <- rep(x_values, times = N)
  df[[x]] <- x_val_rep
  if (!is.null(x_ns)) {
    t0 <- get_teff_obs(data, group_by, x, x_ns)
    t0 <- as.numeric(t0)
    df[[x_ns]] <- x_val_rep - rep(t0, each = k)
  }
  rownames(df) <- NULL
  return(df)
}

#' Get observed effect times from a data frame
#'
#' @inheritParams new_x
#' @return a named vector, where the names are the levels of \code{group_by}
#' @family data utilities
get_teff_obs <- function(data, group_by = "id", x = "age",
                         x_ns = "diseaseAge") {
  check_type(data, "data.frame")
  df <- pick_one_row_each(data, "id")
  times <- dollar(df, x) - dollar(df, x_ns)
  names(times) <- dollar(df, group_by)
  return(times)
}

#' For each unique value of a factor, pick one row from data
#'
#' @param data a data frame
#' @param fac name of a factor in \code{data}
#' @return a data frame
#' @family data utilities
pick_one_row_each <- function(data, fac) {
  check_type(data, "data.frame")
  z <- dollar(data, fac)
  check_type(z, "factor")
  rows <- c()
  for (lev in levels(z)) {
    inds <- which(z == lev)
    rows <- c(rows, inds[1])
  }
  data[rows, ]
}

#' Select data columns which are factors or have a certain name
#'
#' @param data a data frame
#' @param valid names of variables that do not need to be factors
#' @return a data frame
#' @family data utilities
select_factors_and <- function(data, valid) {
  check_type(data, "data.frame")
  col_inds <- c()
  D <- ncol(data)
  nams <- names(data)
  for (j in seq_len(D)) {
    a <- is.factor(data[, j])
    b <- nams[j] %in% valid
    if (a || b) {
      col_inds <- c(col_inds, j)
    }
  }
  data[, col_inds]
}


#' Get type of each data column
#'
#' @param data a data frame
#' @return a named list
#' @family data utilities
data_types <- function(data) {
  check_type(data, "data.frame")
  allowed <- c("factor", "numeric")
  D <- ncol(data)
  nams <- names(data)
  types <- list()
  for (j in seq_len(D)) {
    check_allowed(class(data[, j])[1], allowed)
    types[[j]] <- class(data[, j])[1]
  }
  names(types) <- nams
  return(types)
}

#' Add a crossing of two factors to a data frame
#'
#' @param df a data frame
#' @param fac1 name of first factor, must be found in \code{df}
#' @param fac2 name of second factor, must be found in \code{df}
#' @param new_name name of the new factor
#' @return a data frame
#' @family data utilities
add_factor_crossing <- function(df, fac1, fac2, new_name) {
  a <- dollar(df, fac1)
  b <- dollar(df, fac2)
  check_not_null(new_name)
  check_type(a, "factor")
  check_type(b, "factor")
  df[[new_name]] <- interaction(a, b, sep = "*")
  return(df)
}

#' Data frame and additional information to long format
#'
#' @param df a data frame in wide format
#' @return a data frame where first column is a factor that determines the
#' column variable name in the original data and second column contains
#' the actual values
#' @family data utilities
to_long_format <- function(df) {
  nam <- colnames(df)
  x <- c()
  J <- ncol(df)
  N <- nrow(df)
  for (j in seq_len(J)) {
    x <- c(x, as.vector(df[, j]))
  }
  fac <- as.factor(rep(nam, each = N))
  out <- data.frame(fac, x)
  return(out)
}

#' Repeat a data frame vertically
#'
#' @param df a data frame
#' @param times number of times to repeat
#' @return a data frame
#' @family data utilities
rep_df <- function(df, times) {
  out <- c()
  for (j in seq_len(times)) out <- rbind(out, df)
  return(out)
}
