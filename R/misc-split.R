
#' Split data into training and test data according to given individuals
#'
#' @export
#' @param data a data frame
#' @param test_ids test data individual identifiers
#' @param id_variable name of id variable
#' @return a \code{list(train, test)}
split_data_by_id <- function(data, test_ids, id_variable = "id") {
  id <- data[[id_variable]]
  i_test <- which(id %in% test_ids)
  ret <- split_data(data, i_test)
  return(ret)
}


#' Split data into training and test data according to time point indices
#'
#' @export
#' @param data a data frame
#' @param test_idx indices of test time points
#' @param id_variable name of id variable
#' @return a \code{list(train, test)}
split_data_by_timepoint <- function(data, test_idx, id_variable = "id") {
  id <- data[[id_variable]]
  uid <- unique(id)
  i_test <- c()
  for (i in uid) {
    inds <- which(id == i)
    i_test <- c(i_test, inds[test_idx])
  }
  ret <- split_data(data, i_test)
  return(ret)
}

#' Split data into training and test data randomly
#'
#' @export
#' @param data a data frame
#' @param p_test desired proportion of test data
#' @param n_test desired number of test data points (if NULL, p_test is used
#' to compute this)
#' @return a \code{list(train, test)}
split_data_random <- function(data, p_test = 0.1, n_test = NULL) {
  n_total <- dim(data)[1]
  if (is.null(n_test)) {
    n_test <- round(p_test * n_total)
  }
  i_test <- sample.int(n_total, size = n_test)
  ret <- split_data(data, i_test)
  return(ret)
}


#' Split data into training and test data by selecting randomly k points
#' from each individual
#'
#' @export
#' @param data a data frame
#' @param n_test desired number of test data points per individual
#' @param id_variable name of id variable
#' @return a \code{list(train, test)}
split_data_random_each <- function(data, n_test = 1, id_variable = "id") {
  id <- data[[id_variable]]
  uid <- unique(id)
  i_test <- c()
  for (i in uid) {
    inds <- which(id == i)
    L <- length(inds)
    i_sel <- sample.int(L, n_test)
    i_test <- c(i_test, inds[i_sel])
  }
  ret <- split_data(data, i_test)
  return(ret)
}

#' Split data into training and test data according to given row indices
#'
#' @param data a data frame
#' @param i_test test data row indices
#' @param sort_ids should the test indices be sorted into increasing order
#' @return a \code{list(train, test)}
split_data <- function(data, i_test, sort_ids = TRUE) {
  if (sort_ids) {
    i_test <- sort(i_test, decreasing = FALSE)
  }
  n_total <- dim(data)[1]
  i_train <- setdiff(1:n_total, i_test)
  data_train <- data[i_train, ]
  data_test <- data[i_test, ]
  ret <- list(
    train = data_train, test = data_test,
    i_train = i_train, i_test = i_test
  )
  return(ret)
}
