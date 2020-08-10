#' Create the covariate matrix that is given to stan
#'
#' @param data the data frame that was passed to \code{lgp}
#' @param varInfo original variable type info
#' @param verbose can this print some info?
#' @param formula the model formula
#' @param types the types returned by \code{\link{check_data}}
#' @return a list
create_covariates_stan <- function(data, varInfo, types, formula, verbose) {

  # Create the design matrix X
  XD <- stan_input_X_and_D(data, varInfo, types, formula, verbose)

  # Do this because as.matrix() sometimes creates character matrices
  X <- data.matrix(XD$X)
  D <- XD$D

  # Get location of NaNs because Stan does not accept NaNs
  if (D[3] > 0) {
    X_notnan <- 1 - is.nan(X[, 3])
    X[which(X_notnan == 0), 3] <- 0
  } else {
    X_notnan <- rep(0, length(X[, 1]))
  }
  X_notnan <- as.vector(X_notnan)

  # If there are still NaNs, throw an error
  n_nans <- sum(is.nan(X))
  if (n_nans > 0) {
    stop("Only the diseaseAge column of the data can contain NaNs!")
  }

  # Return
  covariates <- list(
    X = X,
    X_notnan = X_notnan,
    D = D
  )

  return(covariates)
}


#' Validate the `data` input to \code{lgp} and resolve covariate types
#'
#' @param data the data frame that was passed to \code{lgp}
#' @param varInfo variable type info
#' @param verbose can this print some info?
#' @return a list
check_data_ <- function(data, varInfo, verbose) {
  cols <- colnames(data)
  idvar <- varInfo$id_variable
  timevar <- varInfo$time_variable
  davar <- varInfo$disAge_variable
  convars <- varInfo$continuous_vars
  catvars <- varInfo$categorical_vars
  offvars <- varInfo$offset_vars
  d <- length(cols)
  types <- rep(-1, d)

  if (is.null(idvar)) {
    stop("ID variable cannot be NULL!")
  }
  if (is.null(timevar)) {
    stop("Time variable cannot be NULL!")
  }

  # Find response variable
  respvar <- varInfo$response_variable
  idx0 <- which(cols == respvar)
  types[idx0] <- 0

  # Data must have columns corresponding to 'idvar' and 'timevar'
  idx1 <- which(cols == idvar)
  idx2 <- which(cols == timevar)
  if (!(idvar %in% cols)) {
    stop("The data frame must contain a column called ", idvar, "!")
  }
  if (!(timevar %in% cols)) {
    stop("The data frame must contain a column called ", timevar, "!")
  }
  types[idx1] <- 1
  types[idx2] <- 2

  # Disease-related age covariate
  if (is.null(davar)) {
    if ("diseaseAge" %in% cols) {
      davar <- "diseaseAge"
      idx3 <- which(cols == davar)
      if (verbose) {
        msg <- paste0(
          "* Interpreting 'diseaseAge' as the ",
          "disease-related age variable.\n"
        )
        cat(msg)
      }
    } else {
      idx3 <- c()
    }
  } else {
    idx3 <- which(cols == davar)
    if (length(idx3) == 0) {
      stop(
        "The given disease-related age variable ", davar,
        " not found in the data!"
      )
    }
  }
  types[idx3] <- 3

  # Get types of remaining covariates
  remain <- setdiff(1:d, c(idx0, idx1, idx2, idx3))
  for (i in remain) {
    x <- data[, i]
    B <- all.equal(as.integer(x), x)
    if (cols[i] %in% offvars) {
      types[i] <- 6
    } else if (cols[i] %in% catvars) {
      types[i] <- 5
    } else if (cols[i] %in% convars) {
      types[i] <- 4
    } else {
      if (is.logical(B)) {
        typestr <- "categorical"
        types[i] <- 5
        catvars <- c(catvars, cols[i])
      } else {
        typestr <- "continuous"
        types[i] <- 4
        convars <- c(convars, cols[i])
      }
      if (verbose) {
        cat("* Covariate '", cols[i], "' resolved to type '",
          typestr, "'.\n",
          sep = ""
        )
      }
    }
  }

  # Update varInfo
  varInfo$disAge_variable <- davar
  varInfo$continuous_vars <- convars
  varInfo$categorical_vars <- catvars

  return(list(types = types, varInfo = varInfo))
}



#' Predictor covariates and types to Stan input
#'
#' @description Reorders covariates and takes only those that are needed
#' @param data a data frame containing the covariates
#' @param formula model formula
#' @param types types of the covariates
#' @param varInfo original variable type info
#' @param verbose can this print some info?
#' @return X and needed types and updated varInfo
stan_input_X_and_D <- function(data, varInfo, types, formula, verbose) {

  # Predictors
  trm <- stats::terms(formula)
  predictors <- attr(trm, "term.labels")

  # Types
  i1 <- which(types == 1)
  i2 <- which(types == 2)
  i3 <- which(types == 3)
  i4 <- which(types == 4)
  i5 <- which(types == 5)
  i6 <- which(types == 6)

  # Order
  order <- c(i1, i2, i3, i4, i5, i6)
  X <- data[, order]
  types <- types[order]

  # Names
  idvar <- varInfo$id_variable
  timevar <- varInfo$time_variable
  cn <- colnames(X)
  used_names <- c(idvar, timevar, predictors)

  # Take only the needed covariates
  i_use <- which(cn %in% used_names)
  X <- X[, i_use]
  types <- types[i_use]

  # Check unused
  unused <- cn[!(cn %in% used_names)]
  if (verbose) {
    if (length(unused) > 0) {
      cat("* The following data columns will not be used: {")
      cat(paste(unused, collapse = ", "))
      cat("}\n")
    }
  }

  # Create D
  D <- rep(0, 6)
  D[1] <- as.numeric(idvar %in% predictors)
  D[2] <- as.numeric(timevar %in% predictors)
  D[3] <- sum(types == 3)
  D[4] <- sum(types == 4)
  D[5] <- sum(types == 5)
  D[6] <- sum(types == 6)

  if (D[2] == 0) {
    if (D[1] != 0) {
      stop(
        "cannot model id effect as time-dependent ",
        "if the time variable is not in model!"
      )
    }
    if (D[5] != 0) {
      stop(
        "cannot model categorical covariate effect as time-dependent",
        " if the time variable is not in model!"
      )
    }
  }

  return(list(X = X, D = D))
}


#' Standardize continuous input variables in X
#' @param X the design matrix
#' @param D the covariate types, a vector of length 6
#' @return updated X and info about scaling
standardize_inputs <- function(X, D) {

  # Standardize age
  x_age <- X[, 2]
  t_m <- mean(x_age)
  t_std <- stats::sd(x_age)

  # Create the function that does the age transform, also store its inverse map
  sclfun_t <- function(t) {
    (t - t_m) / t_std
  }
  sclfun_t_inv <- function(t) {
    t * t_std + t_m
  }
  TSCL <- list(fun = sclfun_t, fun_inv = sclfun_t_inv)
  X[, 2] <- sclfun_t(x_age)

  # Standardize other continuous covariates
  idx <- 2 + D[3]
  inds <- (idx + 1):(idx + D[4])
  CSCL <- c()
  cntr <- 0
  if (D[4] > 0) {
    for (j in inds) {
      cntr <- cntr + 1
      xj <- X[, j]
      m_j <- mean(xj)
      sd_j <- stats::sd(xj)
      if (sd_j == 0) {
        stop("a continuous covariate has zero variance!")
      }
      sclfun <- function(x) {
        (x - m_j) / sd_j
      }
      sclfun_inv <- function(x) {
        x * sd_j + m_j
      }
      cscl <- list(fun = sclfun, fun_inv = sclfun_inv)
      CSCL[[cntr]] <- cscl
      X[, j] <- sclfun(xj)
    }
  }
  return(list(X = X, TSCL = TSCL, CSCL = CSCL))
}
