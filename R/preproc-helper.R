#' Check that given data is consistent with the given formula
#'
#' @param formula an object of class \code{\link{lgpformula}}
#' @param data an object of class \code{data.frame}
check_data <- function(formula, data) {
  trm <- stats::terms(formula)
  vars <- rownames(attr(trm, "factors"))
  resp <- attr(trm, "response")
  if (!resp) stop("The formula does not contain a response variable")
  tord <- attr(trm, "order")
  if (sum(tord > 1) > 0) {
    stop("Only first-order terms are allowed in the model formula!")
  }
  yName <- vars[resp]
  if (!(yName %in% colnames(data))) {
    stop(paste("The data frame does not contain the response variable", yName))
  }
  vars_seq <- seq_len(length(vars))
  for (i in vars_seq) {
    if (!(vars[i] %in% colnames(data))) {
      stop(paste(
        "Variable", vars[i], "not found in the data frame!",
        "Type ?lgp for help."
      ))
    }
  }
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


#' Get the (scaled) response variable
#'
#' @description Gets and possibly scales the response variable.
#' @param data the data frame given as input to \code{lgp}
#' @param varInfo variable type info
#' @param standardize should the response be standardized to
#' unit variance and zero mean
#' @param LH likelihood as integer
#' @return a list with the (scaled) response variable
#'
get_response <- function(data, varInfo, standardize, LH) {
  yName <- varInfo$response_variable
  response <- data[yName]
  response <- unlist(response)
  y_m <- mean(response)
  y_std <- stats::sd(response)
  if (y_std == 0) {
    stop("The response has zero variance!")
  }

  # Do some checks and update info
  lh_not_01 <- !(LH %in% c(0, 1))
  if (standardize) {
    if (lh_not_01) {
      stop(
        "Standardization of response is only possible if ",
        "likelihood is 'Gaussian' or 'none'!"
      )
    }
  }

  # Create the function that does the transform, also store its inverse map
  if (standardize) {
    sclfun <- function(y) {
      (y - y_m) / y_std
    }
    sclfun_inv <- function(y) {
      y * y_std + y_m
    }
  } else {
    sclfun <- function(y) {
      y
    }
    sclfun_inv <- function(y) {
      y
    }
  }

  # Apply the response scaling
  response <- sclfun(response)

  # Check the response for negative values or non-integer values
  if (lh_not_01) {
    if (sum(response < 0) > 0) {
      msg <- paste0(
        "The response variable contains negative values. ",
        "Only the likelihoods 'Gaussian' and 'none' are allowed ",
        "in such case!\n"
      )
      stop(msg)
    }
    notint <- sum(response - round(response))
    if (notint > 0) {
      msg <- paste0(
        "The response variable contains non-integer values. ",
        "Only the likelihoods 'Gaussian' and 'none' are allowed ",
        "in such case!\n"
      )
      stop(msg)
    }
  }

  # Return
  ret <- list(
    response = response,
    SCL = list(fun = sclfun, fun_inv = sclfun_inv)
  )
  return(ret)
}


#' Get some dimension variables that the Stan model needs as input
#'
#' @param X the design matrix
#' @param D a vector of length 6
#' @return a list
get_model_dims <- function(X, D) {
  N_tot <- length(unique(X[, 1]))
  n <- dim(X)[1]
  d <- dim(X)[2]
  ret <- list(n = n, d = d, D = D, N_tot = N_tot)
  return(ret)
}


#' Check that variable types make sense
#'
#' @param varInfo a named list
#' @return nothing
check_varInfo <- function(varInfo) {

  # Check that id variable is not in offsets
  if (varInfo$id_variable %in% varInfo$offset_vars) {
    stop(
      "The subject identifier variable cannot currently be included in",
      " 'offset_vars'. If you wish to model the effect of 'id_variable'",
      " as a constant offset,",
      " you can create another covariate with the same values and",
      " use it in your 'formula' and 'offset_vars' instead."
    )
  }

  # TODO: more checks
}


#' Count numbers of different categories for each categorical variable
#'
#' @param X the design matrix
#' @param D a vector of length 6
#' @return a numeric vector
set_N_cat <- function(X, D) {
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



#' Set C_hat (Non-gaussian observation models)
#'
#' @param C_hat the \code{C_hat} argument given as input to \code{lgp_model}
#' @param response response variable
#' @param LH likelihood as int
#' @param N_trials the N_trials data (binomial likelihood)
#' @return a real number
set_C_hat <- function(C_hat, response, LH, N_trials) {
  nb_or_pois <- LH %in% c(2, 3)
  binomial <- LH == 4
  if (is.null(C_hat)) {
    if (nb_or_pois) {
      C_hat <- log(mean(response))
    } else if (binomial) {
      p <- mean(response / N_trials)
      C_hat <- log(p / (1 - p))
    } else {
      C_hat <- 0
    }
  } else {
    if (LH == 1) {
      stop(
        "Only give the C_hat argument if observation model is not Gaussian!",
        " With Gaussian likelihood, you should use",
        " C_hat = NULL, in which case the GP mean will be set to zero",
        " and the response variable is standardized to have mean zero."
      )
    }
  }
  n <- length(response)
  L <- length(C_hat)
  if (L != 1 && L != n) {
    stop("C_hat must have length 1 or equal to the number of data points!")
  }
  if (L == 1) {
    C_hat <- rep(C_hat, n)
  }
  return(C_hat)
}


#' Set N_trials (binomial and Bernoulli observation models)
#'
#' @param N_trials the \code{N_trials} argument given as input to
#' \code{lgp_model}
#' @param response response variable
#' @param LH likelihood as int
#' @return a numeric vector
set_N_trials <- function(N_trials, response, LH) {
  if (is.null(N_trials)) {
    N_trials <- rep(1, length(response))
  } else {
    if (LH != 4) {
      stop("Only give the N_trials argument if likelihood is binomial!")
    }
    if (length(N_trials) == 1) {
      N_trials <- rep(N_trials, length(response))
    }
    if (length(N_trials) != length(response)) {
      stop("Invalid length of N_trials!")
    }
  }
  return(N_trials)
}
