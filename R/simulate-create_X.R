#' Simulate an input data frame X
#'
#' @param N Number of individuals.
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
#' @param names Covariate names.
#' @param n_categs An integer vector defining the number of categories
#' for each categorical covariate, so that \code{length(n_categs)} equals to
#' the number of 2's and 3's in the \code{covariates} vector.
#' @param t_data Measurement times.
#' @param t_jitter Standard deviation of the jitter added to the given
#' measurement times.
#' @param t_effect_range Time interval from which the disease effect times are
#' sampled uniformly. Alternatively, This can any function that returns the
#' (possibly randomly generated) real disease effect time for one individual.
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
#' @param verbose verbosity mode
#' @return \code{list(X, onsets, par_cont)}
create_X <- function(N,
                     covariates,
                     names,
                     n_categs,
                     t_data,
                     t_jitter,
                     t_effect_range,
                     continuous_info,
                     verbose) {
  onset_range <- t_effect_range
  D <- rep(0, 5)
  D[1] <- sum(covariates == 0)
  D[2] <- sum(covariates == 1)
  D[3] <- sum(covariates == 2)
  D[4] <- sum(covariates == 3)
  D[5] <- sum(covariates == 4)
  if (length(n_categs) != sum(D[3:4])) {
    stop(paste("Length of n_categs must be same as the number of 2's",
      " and 3's in the covariates vector!",
      " Found = ", length(n_categs), " should be = ",
      sum(D[3:4]), ".",
      sep = ""
    ))
  }

  n_invalid <- sum(n_categs <= 1)
  if (n_invalid > 0) stop("n_categs must only contain integers larger than 1!")
  k <- length(t_data)
  id <- rep(1:N, each = k)
  age <- drawMeasurementTimes(N, t_data, t_jitter)

  # Parse onset range
  if (is.character(onset_range)) {
    if (onset_range == "auto") {
      ran <- range(t_data)
      mrn <- mean(ran)
      onset_range <- c(mrn, mrn)
      if (D[1] == 1) {
        if (verbose) {
          cat("Disease effect time range set to [", onset_range[1], ", ",
            onset_range[2], "]. \n",
            sep = ""
          )
        }
      }
    } else {
      stop("invalid t_effect_range!")
    }
  }

  # Id and onsets
  if (D[1] > 0) {
    N_cases <- round(N / 2)
    if (is.function(onset_range)) {
      onsets <- rep(0, N_cases)
      for (j in 1:N_cases) {
        onsets[j] <- onset_range()
      }
    } else {
      onsets <- stats::runif(N_cases, min = onset_range[1],
                             max = onset_range[2])
    }
    onsets <- c(onsets, rep(NaN, N - N_cases))
    disAge <- onsetsToDiseaseAge(onsets, age, k)
  } else {
    onsets <- NULL
    disAge <- c()
    N_cases <- NaN
  }

  # Other
  if (D[2] > 0) {
    CONT <- drawContinuous(N, k, D[2],
                           continuous_info$mu, continuous_info$lambda)
    cont <- CONT$C
    par_cont <- list(A = CONT$A, B = CONT$B, OFS = CONT$OFS)
  } else {
    cont <- c()
    par_cont <- NULL
  }
  if ((D[3] + D[4]) > 0) {
    categ <- drawCategorical(N, k, n_categs)
  } else {
    categ <- c()
  }
  X <- cbind(id, age, disAge, cont, categ)
  if (D[5] > 0) {
    if (D[1] == 0) {
      stop("cannot include a 4 in covariates if 0 is not included!")
    }
    if (D[5] == 1) {
      group <- as.numeric(!is.nan(disAge))
      X <- cbind(X, group)
    } else {
      stop("only one 4 can be included in the covariates vector")
    }
  }
  colnames(X) <- names
  out <- list(
    X = data.frame(X),
    onsets = onsets,
    par_cont = par_cont,
    N_cases = N_cases
  )
  return(out)
}


#'  Draw the age covariate
#'
#' @param N number of individuals
#' @param t_data a vector of length \code{k}
#' @param t_jitter Standard deviation of the jitter added to the given
#' measurement times.
#' @return a vector of length \code{N*k}
drawMeasurementTimes <- function(N, t_data, t_jitter) {
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
onsetsToDiseaseAge <- function(onsets, age, k) {
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
drawContinuous <- function(N, k, D, mu, lambda) {
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
#' @return a matrix of size \code{N} x \code{D}, where \code{D <- length(v)}
drawCategorical <- function(N, k, v) {
  D <- length(v)
  C <- matrix(0, N, D)
  for (i in c(1:D)) {
    C[, i] <- sample.int(v[i], N, replace = T)
  }
  C_seq <- seq_len(nrow(C))
  rows <- rep(C_seq, each = k)
  C <- C[rows, ]
  return(C)
}
