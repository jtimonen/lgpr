
#' Get some variables related to diseased individuals
#'
#' @param D an integer vector of length 6
#' @param X the design matrix
#' @param X_notnan a binary vector of length n
#' @param uncertain_effect_time Boolean value
#' @param equal_effect Boolean value
#' @param TSCL time scaling function and its inverse
#' @return a list
get_diseased_info <- function(D, X, X_notnan, uncertain_effect_time,
                              equal_effect, TSCL) {
  X_id <- X[, 1]
  test_id_variable(X_id)
  M_max <- max(table(X_id))
  MAPS <- get_case_row_mappings(X_notnan, X_id)
  ONS <- get_onset_info(D, X, MAPS, TSCL)

  # Return
  ret <- list(
    HMGNS = as.numeric(equal_effect),
    UNCRT = as.numeric(uncertain_effect_time),
    N_cases = length(MAPS$caseID_nrows),
    caseID_to_rows = MAPS$caseID_to_rows,
    caseID_nrows = MAPS$caseID_nrows,
    row_to_caseID = MAPS$row_to_caseID,
    M_max = M_max,
    T_observed = ONS$t_obs,
    T_last = ONS$t_last
  )
  return(ret)
}

#' Test that the id variable is correctly specified to allow
#' disease effect modeling
#'
#' @param X_id the id covariate in X
#' @return a list
test_id_variable <- function(X_id) {
  if (!is.numeric(X_id)) {
    stop("id_variable must be numeric (integers)!")
  }
  if (is.unsorted(X_id)) {
    stop("Values of id_variable must be increasing integers!")
  }
  m <- min(X_id)
  if (m != 1) {
    stop("smallest value of id_variable must be 1!")
  }
  M <- max(X_id)
  uid <- unique(X_id)
  N <- length(uid)
  if (M != N) {
    stop(
      paste0(
        "id_variable should take values from the set ",
        "{1, 2, ..., N}, where N is the number of ",
        "unique individuals (", N, ")"
      )
    )
  }
}

#' Create case ID to rows and back mappings
#'
#' @description Create mappings
#' \itemize{
#'   \item from case ID to data rows (caseID_to_rows, caseID_nrows)
#'   \item  from row number to case ID (row_to_caseID)
#' }
#' @param X_notnan binary vector indicating if diseaseAge is available
#' for that measurement
#' @param X_id the id covariate in X
#' @param only_R2C should this return only the rows-to-caseID mapping
#' @return a list
get_case_row_mappings <- function(X_notnan, X_id, only_R2C = FALSE) {
  UNI <- unique(cbind(X_id, X_notnan))
  uid <- UNI[, 1]
  is_case <- UNI[, 2]
  N_cases <- sum(is_case)
  inn <- which(X_notnan == 1)
  tab <- table(X_id[inn])
  nrows <- as.numeric(tab)
  M_max <- max(table(X_id))

  # Create mappings
  # - from case ID to data rows (C2R, nrows)
  # - from row number to case ID (R2C)
  C2R <- matrix(0, N_cases, M_max)
  R2C <- rep(0, length(X_id))
  i_case <- 0
  j <- 0
  for (id in uid) {
    j <- j + 1
    if (is_case[j] == 1) {
      i_case <- i_case + 1
      i_rows <- which(X_id == id)
      if (!only_R2C) {
        C2R[i_case, 1:nrows[i_case]] <- i_rows
      }
      R2C[i_rows] <- i_case
    }
  }

  # Return
  if (only_R2C) {
    return(R2C)
  } else {
    ret <- list(
      caseID_to_rows = C2R,
      caseID_nrows = nrows,
      row_to_caseID = R2C
    )
    return(ret)
  }
}


#' Get disease onset info
#'
#' @description This returns
#' \itemize{
#'   \item a vector of observed onsets
#'   \item mapping from case ID to average sampling interval before
#'   the observed disease onset
#' }
#' @param D an integer vector of length 6
#' @param X the design matrix
#' @param MAPS mappings created by \code{get_case_row_mappings}
#' @param TSCL time scaling function and its inverse
#' @return two vectors of length \code{N_cases}
get_onset_info <- function(D, X, MAPS, TSCL) {
  M <- MAPS$caseID_nrows
  C2R <- MAPS$caseID_to_rows
  if (D[3] == 0) {
    A <- rep(0, 0)
    L <- rep(0, 0)
  } else {
    X_age <- as.numeric(X[, 2])
    X_age <- TSCL$fun_inv(X_age) # age back to original scale
    if (any(X_age <= 0)) {
      stop("There seems to be a measurement at age <=0 !")
    }
    X_disAge <- as.numeric(X[, 3])
    N_cases <- length(M)
    A <- rep(0, N_cases)
    L <- rep(0, N_cases)
    for (k in 1:N_cases) {
      i_row <- C2R[k, 1:M[k]]
      a <- X_age[i_row]
      d <- X_disAge[i_row]
      A[k] <- a[1] - d[1]
      L[k] <- max(a)
    }
  }
  return(list(t_obs = A, t_last = L))
}
