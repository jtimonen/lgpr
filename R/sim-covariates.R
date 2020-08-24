#' Simulate an input data frame X
#'
#' @param names Covariate names.
#' @inheritParams sim_create_x_D
#' @inheritParams sim_create_x_check
#' @inheritParams sim_create_x_dis_age
#' @inheritParams sim_create_x_other
#' @param t_jitter Standard deviation of the jitter added to the given
#' measurement times.
#' @return a list
sim_create_x <- function(N,
                         covariates,
                         names,
                         n_categs,
                         t_data,
                         t_jitter,
                         t_effect_range,
                         continuous_info) {

  # Validate input
  D <- sim_create_x_D(covariates)
  checked <- sim_create_x_check(n_categs, D, t_data, t_effect_range)

  # Initialize data
  k <- length(t_data)
  id <- rep(1:N, each = k)
  id <- as.factor(id)
  age <- sim_draw_measurement_times(N, t_data, t_jitter)
  X <- data.frame(id, age)

  # Effect times and disease ages
  parsed_dis <- sim_create_x_dis_age(X, k, N, D, age, checked$t_effect_range)
  dis_age <- parsed_dis$dis_age
  X <- parsed_dis$X

  # Other covariates
  cinfo <- continuous_info
  parsed_x <- sim_create_x_other(X, k, N, D, n_categs, dis_age, cinfo)
  X <- parsed_x$X
  colnames(X) <- names

  # Return
  list(
    info = checked$info,
    onsets = parsed_dis$teff,
    N_cases = parsed_dis$N_cases,
    par_cont = parsed_x$par_cont,
    X = X
  )
}

#' Helper function for sim_create_x
#'
#' @param covariates Integer vector that defines the types of covariates
#' (other than id and age). If not given, only the id and age
#' covariates are created. Different integers correspond to the following
#' covariate types:
#' \itemize{
#'   \item 0 = disease-related age
#'   \item 1 = other continuous covariate
#'   \item 2 = a categorical covariate that interacts with age
#'   \item 3 = a categorical covariate that acts as a group offset
#'   \item 4 = a categorical covariate that that acts as a group offset AND
#'   is restricted to have value 0 for controls and 1 for cases
#' }
#' @return an array of five integers
sim_create_x_D <- function(covariates) {
  D <- rep(0, 5)
  D[1] <- sum(covariates == 0)
  D[2] <- sum(covariates == 1)
  D[3] <- sum(covariates == 2)
  D[4] <- sum(covariates == 3)
  D[5] <- sum(covariates == 4)
  return(D)
}

#' Input check helper function for sim_create_x
#'
#' @param n_categs An integer vector defining the number of categories
#' for each categorical covariate, so that \code{length(n_categs)} equals to
#' the number of 2's and 3's in the \code{covariates} vector.
#' @param D covariate type array
#' @param t_data Measurement times.
#' @param t_effect_range Time interval from which the disease effect times are
#' sampled uniformly. Alternatively, This can any function that returns the
#' (possibly randomly generated) real disease effect time for one individual.
#' @returns a list
sim_create_x_check <- function(n_categs, D, t_data, t_effect_range) {

  # Check length of n_categs
  L1 <- length(n_categs)
  L2 <- sum(D[3:4])
  if (L1 != L2) {
    msg <- paste0(
      "Length of <n_categs> must be same as the number of 2's",
      " and 3's in the <covariates> vector!",
      " Found = ", L1, ", should be = ", L2, "."
    )
    stop(msg)
  }

  # Check values of n_categs
  n_invalid <- sum(n_categs <= 1)
  if (n_invalid > 0) {
    stop("<n_categs> must only contain integers larger than 1!")
  }

  # Parse effect time range
  info <- ""
  if (is.character(t_effect_range)) {
    if (t_effect_range == "auto") {
      ran <- range(t_data)
      mrn <- mean(ran)
      t_effect_range <- c(mrn, mrn)
      if (D[1] == 1) {
        info <- paste0(
          "Disease effect time range set to [",
          t_effect_range[1], ", ", t_effect_range[2], "].\n"
        )
      }
    } else {
      stop("invalid t_effect_range!")
    }
  }

  # Return
  list(
    t_effect_range = t_effect_range,
    info = info
  )
}

#' Helper function for sim_create_x
#'
#' @param X current covariate matrix
#' @param k number of time points per individual
#' @param N Number of individuals.
#' @param D covariate type array
#' @param age the age covariate
#' @param t_effect_range Parsed version of \code{t_effect_range}
#' @param return list
sim_create_x_dis_age <- function(X, k, N, D, age, t_effect_range) {
  if (D[1] > 0) {
    N_cases <- round(N / 2)
    if (is.function(t_effect_range)) {
      teff <- rep(0, N_cases)
      for (j in 1:N_cases) {
        teff[j] <- t_effect_range()
      }
    } else {
      teff <- stats::runif(N_cases,
        min = t_effect_range[1],
        max = t_effect_range[2]
      )
    }
    teff <- c(teff, rep(NaN, N - N_cases))
    dis_age <- sim_onsets_to_dis_age(teff, age, k)
    X <- cbind(X, dis_age)
  } else {
    dis_age <- NULL
    teff <- numeric()
    N_cases <- NaN
  }

  # Return
  list(
    X = X,
    teff = teff,
    N_cases = N_cases,
    dis_age = dis_age
  )
}



#' Helper function for sim_create_x
#'
#' @inheritParams sim_create_x_dis_age
#' @inheritParams sim_create_x_check
#' @param continuous_info Info for generating continuous covariates. Must be a
#' list containing fields \code{lambda} and \code{mu}, which have length 3.
#' The continuous covariates are generated so that \code{x <- sin(a*t + b) + c},
#' where
#' \itemize{
#'   \item \code{t <- seq(0, 2*pi, length.out = k)}
#'   \item \code{a <- mu[1] + lambda[1]*stats::runif(1)}
#'   \item \code{b <- mu[2] + lambda[2]*stats::runif(1)}
#'   \item \code{c <- mu[3] + lambda[3]*stats::runif(1)}
#' }
#' @param dis_age a disease-age vector
#' @returns a list
sim_create_x_other <- function(X, k, N, D, n_categs, dis_age, continuous_info) {
  if (D[2] > 0) {
    mu <- continuous_info$mu
    lambda <- continuous_info$lambda
    CONT <- sim_draw_continuous(N, k, D[2], mu, lambda)
    cont <- CONT$C
    par_cont <- list(a = CONT$A, b = CONT$B, offset = CONT$OFS)
    X <- cbind(X, cont)
  } else {
    par_cont <- list(a = NULL, b = NULL, offset = NULL)
  }

  if ((D[3] + D[4]) > 0) {
    categ <- sim_draw_categorical(N, k, n_categs)
    X <- cbind(X, categ)
  }

  if (D[5] > 0) {
    if (D[1] == 0) {
      stop("cannot include a 4 in covariates if 0 is not included!")
    }
    if (D[5] == 1) {
      group <- as.factor(as.numeric(!is.nan(dis_age)))
      X <- cbind(X, group)
    } else {
      stop("only one 4 can be included in the covariates vector")
    }
  }

  # Return
  list(X = X, par_cont = par_cont)
}


#'  Draw the age covariate
#'
#' @param N number of individuals
#' @param t_data a vector of length \code{k}
#' @param t_jitter Standard deviation of the jitter added to the given
#' measurement times.
#' @return a vector of length \code{N*k}
sim_draw_measurement_times <- function(N, t_data, t_jitter) {
  k <- length(t_data)
  age <- rep(0, N * k)
  if (t_jitter < 0) {
    stop("t_jitter must be positive!")
  }
  for (i in 1:N) {
    idx <- (1 + (i - 1) * k):(i * k)
    t <- t_data + stats::rnorm(k, mean = 0, sd = t_jitter)
    age[idx] <- sort(t)
  }
  return(age)
}


#' Compute the disease-related ages
#' @param age the age covariate, a vector of length \code{N*k}
#' @param onsets true disease effect times, a vector of length \code{N}
#' @param k number of measurements per individual
#' @return the diseaseAge covariate, a vector of length \code{N*k}
sim_onsets_to_dis_age <- function(onsets, age, k) {
  N <- length(onsets)
  n <- N * k
  diseaseAge <- rep(0, n)
  for (i in 1:N) {
    inds <- ((i - 1) * k + 1):(i * k)
    if (!is.nan(onsets[i])) {
      diseaseAge[inds] <- age[inds] - onsets[i]
    } else {
      diseaseAge[inds] <- NaN
    }
  }
  return(diseaseAge)
}


#' Indepedently draw continuous variables for each individual
#'
#' @param N number of individuals
#' @param k number of timepoints
#' @param D number of variables
#' @param mu a vector of length 3
#' @param lambda a vector of length 3
#' @return a matrix of size \code{N} x \code{D}
sim_draw_continuous <- function(N, k, D, mu, lambda) {
  C <- matrix(0, N * k, D)
  A <- matrix(0, N, D)
  B <- matrix(0, N, D)
  OFS <- matrix(0, N, D)
  for (j in 1:D) {
    for (i in c(1:N)) {
      inds <- ((i - 1) * k + 1):(i * k)
      ttt <- seq(0, 2 * pi, length.out = k)
      a <- mu[1] + lambda[1] * stats::runif(1)
      b <- mu[2] + lambda[2] * stats::runif(1)
      ofs <- mu[3] + lambda[3] * stats::runif(1)
      A[i, j] <- a
      B[i, j] <- b
      OFS[i, j] <- ofs
      C[inds, j] <- sin(a * ttt + b) + ofs
    }
  }
  return(list(A = A, B = B, OFS = OFS, C = C))
}


#' Independently draw categorical variables for each individual
#'
#' @param N number of individuals
#' @param k number of timepoints
#' @param v vector of numbers of different categories
#' @return a data frame with \code{N} rows and \code{D} = \code{length(v)}
#' columns
sim_draw_categorical <- function(N, k, v) {
  D <- length(v)
  C <- matrix(0, N, D)
  for (i in seq_len(D)) {
    C[, i] <- sample.int(v[i], N, replace = TRUE)
  }
  C_seq <- seq_len(nrow(C))
  rows <- rep(C_seq, each = k)
  C <- data.frame(C[rows, ])
  for (i in seq_len(D)) {
    C[, i] <- as.factor(C[, i])
  }
  return(C)
}
