#' Validate the formula of \code{lgp}
#'
#' @description Checks if the input `formula` to \code{lgp_model} 
#' are valid with the given data
#' @param formula the formula that was passed to \code{lgp_model}
#' @param data the data frame that was passed to \code{lgp_model}
#' @return nothing
check_formula <- function(formula, data){
  
  trm  <- stats::terms(formula)
  vars <- rownames(attr(trm, "factors"))
  resp <- attr(trm, "response")
  if(!resp) stop("The formula does not contain a response variable")
  comp  <- attr(trm, "term.labels")
  tord  <- attr(trm, "order")
  if(sum(tord > 1) > 0) stop("Only first-order terms are allowed in the model formula!")
  yName <- vars[resp]
  if(!(yName  %in% colnames(data))) stop(paste('The data frame does not contain the response variable', yName))
  for(i in 1:length(vars)){
    if(!(vars[i] %in% colnames(data))) stop(paste("Variable", vars[i], "not found in the data frame!"))
  }
}


#' Validate the `data` input to \code{lgp} and resolve covariate types
#'
#' @param data the data frame that was passed to \code{lgp}
#' @param varInfo variable type info 
#' @param verbose can this print some info?
#' @return a list 
check_data <- function(data, varInfo, verbose){
  
  cols    <- colnames(data)
  idvar   <- varInfo$id_variable
  timevar <- varInfo$time_variable
  davar   <- varInfo$disAge_variable
  convars <- varInfo$continuous_vars
  catvars <- varInfo$categorical_vars
  offvars <- varInfo$offset_vars
  d       <- length(cols)
  types   <- rep(-1, d)
  
  if(is.null(idvar)){stop("ID variable cannot be NULL!")}
  if(is.null(timevar)){stop("Time variable cannot be NULL!")}
  
  # Find response variable
  respvar     <- varInfo$response_variable
  idx0        <- which(cols==respvar)
  types[idx0] <- 0
  
  # Data must have columns corresponding to'idvar' and 'timevar'
  idx1    <- which(cols==idvar)
  idx2    <- which(cols==timevar)
  if(!(idvar  %in% cols)) stop('The data frame must contain a column called ', idvar, '!')
  if(!(timevar %in% cols)) stop('The data frame must contain a column called ', timevar, '!')
  types[idx1] <- 1
  types[idx2] <- 2
  
  # Disease-related age covariate
  if(is.null(davar)){
    if("diseaseAge" %in% cols){
      davar <- "diseaseAge"
      idx3  <- which(cols==davar)
      if(verbose){
        cat("* Interpreting 'diseaseAge' as the disease-related age variable.\n")
      }
    }else{
      idx3 <- c()
    }
  }else{
    idx3 <- which(cols==davar)
    if(length(idx3)==0){
      stop("The given disease-related age variable ", davar, " not found in the data!")
    }
  }
  types[idx3] <- 3
  
  # Get types of remaining covariates
  remain <- setdiff(1:d, c(idx0,idx1,idx2,idx3))
  for(i in remain){
    x <- data[,i]
    B <- all.equal(as.integer(x),x)
    if(cols[i] %in% offvars){
      types[i] <- 6
    }else if(cols[i] %in% catvars){
      types[i] <- 5
    }else if(cols[i] %in% convars){
      types[i] <- 4
    }else{
      if(is.logical(B)){
        typestr  <- "categorical"
        types[i] <- 5
        catvars  <- c(catvars, cols[i])
      }else{
        typestr  <- "continuous"
        types[i] <- 4
        convars  <- c(convars, cols[i])
      }
      if(verbose){
        cat("* Covariate '", cols[i], "' resolved to type '",
            typestr,"'.\n", sep="")
      }
    }
  }
  
  # Update varInfo
  varInfo$disAge_variable  <- davar
  varInfo$continuous_vars  <- convars
  varInfo$categorical_vars <- catvars
  
  return(list(types=types, varInfo=varInfo))
  
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
stan_input_X_and_D <- function(data, varInfo, types, formula, verbose){
  
  # Predictors
  trm <- stats::terms(formula)
  predictors <- attr(trm, "term.labels")
  
  # Types
  i0 <- which(types==0)
  i1 <- which(types==1)
  i2 <- which(types==2)
  i3 <- which(types==3)
  i4 <- which(types==4)
  i5 <- which(types==5)
  i6 <- which(types==6)
  
  # Order
  order <- c(i1,i2,i3,i4,i5,i6)
  X     <- data[,order]
  types <- types[order]
  
  # Names
  idvar      <- varInfo$id_variable
  timevar    <- varInfo$time_variable
  cn         <- colnames(X)
  used_names <- c(idvar, timevar, predictors)
  
  # Take only the needed covariates
  i_use <- which(cn %in% used_names)
  used  <- cn[i_use]
  X     <- X[,i_use]
  types <- types[i_use]
  
  # Check unused
  unused <- cn[!(cn %in% used_names)]
  if(verbose){
    if(length(unused)>0){
      cat("* The following data columns will not be used: {")
      cat(paste(unused, collapse=", "))
      cat("}\n")
    }
  }
  
  # Create D
  D    <- rep(0, 6)
  D[1] <- as.numeric(idvar %in% predictors)
  D[2] <- as.numeric(timevar %in% predictors)
  D[3] <- sum(types==3)
  D[4] <- sum(types==4)
  D[5] <- sum(types==5)
  D[6] <- sum(types==6)
  
  if(D[2]==0){
    if(D[1]!=0){
      stop("cannot model id effect as time-dependent if the time variable is not in model!")
    }
    if(D[5]!=0){
      stop("cannot model categorical covariate effect as time-dependent",
           " if the time variable is not in model!")
    }
  }
  
  return(list(X = X, D = D))
}


#' Standardize continuous input variables in X
#' @param X the design matrix
#' @param D the covariate types, a vector of length 6
#' @return updated X and info about scaling
standardize_inputs <- function(X, D){
  
  # Standardize age
  x_age     <- X[,2]
  t_m       <- mean(x_age)
  t_std     <- stats::sd(x_age)
  
  # Create the function that does the age transform, also store its inverse map
  sclfun_t     <- function(t){(t - t_m)/t_std}
  sclfun_t_inv <- function(t){t*t_std + t_m}
  TSCL         <- list(fun = sclfun_t, fun_inv = sclfun_t_inv)
  X[,2]        <- sclfun_t(x_age)
  
  # Standardize other continuous covariates
  idx   <- 2 + D[3]
  inds  <- (idx + 1):(idx + D[4])
  CSCL  <- c()
  cntr  <- 0
  if(D[4] > 0){
    for(j in inds){
      cntr         <- cntr + 1
      xj           <- X[,j]
      m_j          <- mean(xj)
      sd_j         <- stats::sd(xj)
      if(sd_j == 0){
        stop("a continuous covariate has zero variance!")
      }
      sclfun       <- function(x){(x - m_j)/sd_j}
      sclfun_inv   <- function(x){x*sd_j + m_j}
      cscl         <- list(fun = sclfun, fun_inv = sclfun_inv)
      CSCL[[cntr]] <- cscl
      X[,j]        <- sclfun(xj)
    }
  }
  return(list(X = X, TSCL = TSCL, CSCL = CSCL))
}


#' Get the (scaled) response variable
#'
#' @description Gets and possibly scales the response variable. 
#' @param data the data frame given as input to \code{lgp}
#' @param varInfo variable type info
#' @param standardize should the response be standardized to unit variance and zero mean
#' @param likelihood the likelihood
#' @return a list with the (scaled) response variable
#'
get_response <- function(data, varInfo, standardize, likelihood){
  
  yName      <- varInfo$response_variable
  response   <- data[yName]
  response   <- unlist(response)
  y_m        <- mean(response)
  y_std      <- stats::sd(response)
  if(y_std==0){
    stop("The response has zero variance!")
  }
  
  # Do some checks and update info
  if(standardize){
    if(!(likelihood %in% c("none", "Gaussian"))){
      stop("Standardization of response is only possible if likelihood is 'Gaussian' or 'none'!")
    }
  }
  
  # Create the function that does the transform, also store its inverse map
  if(standardize){
    sclfun     <- function(y){(y - y_m)/y_std}
    sclfun_inv <- function(y){y*y_std + y_m}
  }else{
    sclfun     <- function(y){y}
    sclfun_inv <- function(y){y}
  }
  
  # Apply the response scaling
  response <- sclfun(response)
  
  # Check the response for negative values or non-integer values
  if(!(likelihood %in% c("none", "Gaussian"))){
    if(sum(response<0) > 0){
      msg <- paste("The response variable contains negative values.",
                   "Only the likelihoods 'Gaussian' and 'none' are allowed in such case!\n", sep="")
      stop(msg)
    }
    notint <- sum(response - round(response))
    if(notint > 0){
      msg <- paste("The response variable contains non-integer values.",
                   " Only the likelihoods 'Gaussian' and 'none' are allowed in such case!\n", sep="")
      stop(msg)
    }
  }
  
  # Return
  ret <- list(response = response,
              SCL      = list(fun = sclfun, fun_inv = sclfun_inv) )
  return(ret)
}


#' Set a lot of generic variables that the Stan model needs as input
#'
#' @param X the design matrix
#' @param D a vector of length 6
#' @param likelihood the `likelihood` input to \code{lgp}
#' @return a list
get_model_dims <- function(X, D, likelihood){
  
  # Dimensions
  N_tot <- length(unique(X[,1]))
  n     <- dim(X)[1]
  d     <- dim(X)[2]
  
  # Get likelihood
  if(likelihood=="none"){
    LH <- 0
  }else if(likelihood=="Gaussian"){
    LH <- 1
  }else if(likelihood=="Poisson"){
    LH <- 2
  }else if(likelihood=="NB"){
    LH <- 3
  }else{
    stop("The likelihood must be either 'none', 'Gaussian', 'Poisson', or 'NB'!")
  }
  
  # Return
  ret <- list(
    n     = n,
    d     = d,
    D     = D,
    N_tot = N_tot,
    LH    = LH
  )
  
  return(ret)
}



