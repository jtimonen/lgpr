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


#' Parse the given data
#'
#' @param data A \code{data.frame} where each column corresponds to one
#' variable, and each row is one observation. Continuous covariates must have
#' type \code{"numeric"} and categorical ones must have type \code{"factor"}.
#' Missing values should be indicated with \code{NaN} or \code{NA}.
#' @return a named list of parsed options
parse_data <- function(data, formula){
  return(list(0))
}
