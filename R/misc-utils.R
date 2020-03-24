
#' Create the GP mean input for \code{lgp}, so that it accounts
#' for normalization between data points in the Poisson or NB
#' observation model
#'
#' @export
#' @param y response variable, vector of length \code{n}
#' @param norm_factors normalization factors, vector of length \code{n}
#' @return a vector of length \code{n}, which can be used as
#' the \code{C_hat} input to the \code{lgp} function
adjusted_Chat <- function(y, norm_factors){
  
  if(length(norm_factors)!=length(y)){stop("inputs must have same length!")}
  if(sum(y<0) > 0){stop("y cannot have negative values!")}
  if(sum(round(y)!=y) > 0){stop("y must have only integer values!")}
  if(sum(norm_factors<=0) > 0){stop("norm_factors must be all positive!")}
  
  C_hat <- log(mean(y))
  C_hat <- rep(C_hat, length(y))
  C_hat <- C_hat + log(norm_factors)
  return(C_hat)
}


#' Repeat a vector as a rows of an array
#'
#' @param v a vector of length \code{m}
#' @param n number of times to repeat
#' @return returns an array of size \code{n} x \code{m}
repvec <- function(v, n){
  m <- length(v)
  A <- matrix(rep(v,n), n, m, byrow = TRUE)
  return(as.array(A))
}


#' Convert the Stan likelihood encoding to a string
#'
#' @param LH an integer
#' @return a string
likelihood_as_str <- function(LH){
  if(LH==1 || LH==-1){
    str <- "Gaussian"
  }else if(LH==0){
    str <- "none"
  }else if(LH==2){
    str <- "Poisson"
  }else if(LH==3){
    str <- "NB"
  }else if(LH==4){
    str <- "binomial"
  }else{
    str <- "Unknown likelihood"
  }
  return(str)
}


#' Convert likelihood string to Stan encoding
#'
#' @param likelihood a string
#' @return an integer
likelihood_as_int <- function(likelihood){
  likelihood <- tolower(likelihood)
  if(likelihood=="none"){
    LH <- 0
  }else if(likelihood=="gaussian"){
    LH <- 1
  }else if(likelihood=="poisson"){
    LH <- 2
  }else if(likelihood=="nb"){
    LH <- 3
  }else if(likelihood=="binomial"){
    LH <- 4
  }else{
    stop("likelihood must be either 'none', 'Gaussian', 'Poisson', 'binomial', or 'NB'!")
  }
  return(LH)
}


#' Get names of model components
#'
#' @param stan_dat The data that was passed to \code{rstan::sampling}
#' @return names of model components
lgp_component_names <- function(stan_dat){
  D        <- stan_dat$D
  compnam  <- rownames(stan_dat$X)
  idvar    <- compnam[1]
  timevar  <- compnam[2]
  
  compnam[1] <- paste(idvar,'*', timevar, sep="")
  idx      <- rep(1, length(compnam))
  idx[1:2] <- c(D[1],D[2])
  types    <- c(1,2, rep(3,D[3]), rep(4,D[4]), rep(5, D[5]), rep(6, D[6]))
  i_inter <- which(types==5)
  for(j in 1:length(i_inter)){
    ind <- i_inter[j]
    compnam[ind] <- paste(compnam[ind],"*",timevar,sep="")
  }
  cn  <- compnam[which(idx==1)]
  cn  <- gsub("\\*", ", ", cn)
  cn  <- paste("f[(",1:length(cn),")](",cn,")", sep="")
  return(cn)
}


#' Get names of model covariates
#'
#' @param stan_dat The data that was passed to \code{rstan::sampling}
#' @return names of model components
lgp_covariate_names <- function(stan_dat){
  D        <- stan_dat$D
  covnam   <- rownames(stan_dat$X)
  idx      <- rep(1, length(covnam))
  idx[1:2] <- c(D[1],D[2])
  covnam  <- covnam[which(idx==1)]
  return(covnam)
}


#' Matrix to data frame without editing column names
#'
#' @param M a matrix
#' @return a data frame
matrix_to_df <- function(M){
  cn  <- colnames(M)
  df  <- data.frame(M)
  colnames(df) <- cn
  return(df)
}


#' Get model info
#'
#' @param object an object of class \code{lgpmodel} or \code{lgpfit}
#' @param print should this print the info?
#' @return the info as a string
model_info <- function(object, print = TRUE){
  if(class(object)=="lgpfit"){
    object <- object@model
  }
  info  <- object@info
  
  # Model formula and likelihood
  str    <- paste('  Model:\n    f = ', paste(info$component_names, collapse=" + "), sep = "")
  LH     <- object@stan_dat$LH
  LH_str <- likelihood_as_str(LH)
  yvar   <- info$varInfo$response_variable
  str    <- paste(str, "\n    Response variable: ", yvar, sep="")
  str    <- paste(str, "\n    Observation model: ", LH_str, sep="")
  str    <- paste(str, "\n")
  
  # Covariate types
  t1 <- info$varInfo$id_variable
  t2 <- info$varInfo$time_variable
  t3 <- info$varInfo$disAge_variable
  t4 <- info$varInfo$continuous_vars
  t5 <- info$varInfo$categorical_vars
  t6 <- info$varInfo$offset_vars
  
  str <- paste(str, "  Variable types:\n", sep="")
  str <- paste(str, "    - Identifier variable: ", t1, "\n", sep = "")
  str <- paste(str, "    - Time variable: ", t2, "\n", sep ="")
  if(!is.null(t3)){
    str <- paste(str, "    - Disease-related age variable: ", t3, "\n", sep="")
  }
  if(!is.null(t4)){
    str <- paste(str, "    - Other continuous variables: ", paste(t4, collapse = ", "), "\n", sep="")
  }
  if(!is.null(t5)){
    str <- paste(str, "    - Other categorical variables: ", paste(t5, collapse = ", "), "\n", sep="")
  }
  if(!is.null(t6)){
    str <- paste(str, "    - Group offset variables: ", paste(t6, collapse = ", "), "\n", sep="")
  }
  
  # Print info
  if(print){
    cat(str)
  }
  return(invisible(str))
}


#' Extract observed disease onset times from diseaseAge covariate vector
#'
#' @param id the id covariate, vector of length \code{n}
#' @param age the age covariate, vector of length \code{n}
#' @param disAge the observed disease-related age covariate, vector of length \code{n}
#' @return vector of observed onset times
get_obs_onset_times <- function(id, age, disAge){
  uid          <- unique(id)
  t_ons        <- rep(0, length(uid))
  formatter    <- function(x){formatC(x, width = 2, format = "d", flag = "0")}
  uid_str      <- formatter(uid) # zero padded format for the IDs
  names(t_ons) <- uid_str
  
  j <- 0
  for(i in uid){
    j    <- j + 1
    inds <- which(id==i)
    i1   <- inds[1]
    if(is.nan(disAge[i1])){
      t_ons[j] <- NaN
    }else{
      t_ons[j] <- age[i1] - disAge[i1]
    }
    
  }
  return(t_ons)
}


#' Extract samples of T_effect
#'
#' @param fit an object of class \code{lgpfit}
#' @return a matrix
extract_t_effect_samples <- function(fit){
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  if(fit@model@stan_dat$UNCRT==0){
    stop("uncertainty was not modeled!")
  }else{
    SMP <- rstan::extract(fit@stan_fit, pars = "T_effect")$T_effect[,1,]
    cid <- get_case_ids(fit)
    formatter <- function(x){formatC(x, width = 2, format = "d", flag = "0")}
    cid_str <- formatter(cid) # zero padded format for the IDs
    colnames(SMP) <- cid_str
    return(SMP)  
  }
}


#' Get case ids in original data
#'
#' @param fit an object of class \code{lgpfit}
#' @return a character vector
get_case_ids <- function(fit){
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  sd  <- fit@model@stan_dat
  id  <- as.numeric(sd$X[1,])
  i1  <- sd$caseID_to_rows[,1]
  cid <- id[i1]
  return(cid)
}


#' Get lgpr version description
#'
#' @export
#' @return package description
get_pkg_description <- function(){
  lgprLib <- dirname(system.file(package = "lgpr"))
  descr   <- suppressWarnings(utils::packageDescription("lgpr",
                                                        lib.loc = lgprLib))
  return(descr)
}


#' Get main stan model of the package
#'
#' @export
#' @return an object of class stanmodel
get_stan_model <- function(){
  return(stanmodels[["lgp"]])
}


#' Get formula of a full model with all covariates included
#' @export
#' @param data a data frame, where the response variable is the last column
#' @return a formula
full_model_formula <- function(data){
  cn   <- colnames(data)
  ddd  <- length(cn)
  rhs  <- paste(cn[1:(ddd-1)], collapse = " + ")
  form <- stats::formula(paste(cn[ddd], "~", rhs))
  return(form)
}

#' Create a full model with all covariates included
#' @export
#' @param data a data frame
#' @param ... additional parameters to \code{\link{lgp_model}}
#' @return a ggplot object
full_model <- function(data, ...){
  form  <- full_model_formula(data)
  model <- lgp_model(form, data, ...)
  return(model)
}

#' PRED object to arrays
#' @param PRED an object returned by \code{\link{lgp_predict}}
#' @return a list containing two arrays
PRED_to_arrays <- function(PRED){
  L <- PRED$LIST
  S <- length(L)
  FFF_1 <- L[[1]]$mu_cmp
  nnn <- dim(FFF_1)[1]
  ddd <- dim(FFF_1)[2]
  MMM <- array(0, c(S,nnn,ddd+1))
  SSS <- array(0, c(S,nnn,ddd+1))
  for(s in 1:S){
    MMM_s <- L[[s]]$mu_cmp
    SSS_s <- sqrt(L[[s]]$s2_cmp)
    mu    <- L[[s]]$mu_f
    std   <- sqrt(L[[s]]$s2_f)
    MMM[s,,] <- as.matrix(cbind(MMM_s, mu))
    SSS[s,,] <- as.matrix(cbind(SSS_s, std))
  }
  return(list(MMM=MMM,SSS=SSS))
}

#' Create an example fit object
#' 
#' @export
#' @param N number of individuals
#' @param t time points
#' @param iter number of iterations
#' @param chains number of chains
#' @param ... other arguments to the \code{\link{lgp}} call
#' @return an object of class \code{\link{lgpfit}}
example_fit <- function(N = 4, 
                        t = 10*c(1,2,3,4,5),
                        iter = 100,
                        chains = 1,
                        ...){
  fit <- lgp(formula = y ~ id + age, 
             data    = simulate_data(N = 4, t)$data, 
             iter    = iter,
             chains  = chains, 
             refresh = 0,
             ...)
  return(fit)
}

#' Returns a valid example call of the \code{lgp} function
#' witi valid data input
#' 
#' @export
#' @return a string
example_call <- function(){
  code <- "
  fit <- lgp(formula = y ~ id + age, 
  data    = simulate_data(N = 4, t_data=10*c(1,2,3,4,5))$data, 
  iter    = 100,
  chains  = 1, 
  refresh = 0,
  verbose = FALSE)
  "
  return(code)
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
add_categorical_covariate <- function(data, x, id_var="id"){
  name     <- deparse(substitute(x))
  if(name %in% colnames(data)){
    stop("The data frame already contains a variable called '", name, "'!")
  }
  x_id     <- as.numeric(names(x))
  data_id  <- data[[id_var]]
  uid      <- unique(x_id)
  xx       <- rep(0, length(data_id))
  for(id in uid){
    i_data <- which(data_id == id)
    i_new <- which(x_id == id)
    xx[i_data] <- x[i_new]
  }
  data[[name]] <- xx
  return(data)
}


#' Create the disease-related age covariate vector based on the
#' disease initiation times and add it to the data frame
#' 
#' @export
#' @param data the original data frame
#' @param t_init A named vector containing the observed initiation or onset
#' time for each individual. The names, i.e. \code{names(t_init)}, should 
#' specify the individual id.
#' @param id_var name of the id variable in \code{data}
#' @param time_var name of the time variable in \code{data}
#' @return A data frame with one column added. The new column will
#' be called \code{'diseaseAge'}. For controls, the value of diseaseAge
#' will be set to NaN.
add_diseaseAges <- function(data, t_init, id_var="id", time_var="age"){
  if("diseaseAge" %in% colnames(data)){
    stop("The data frame already contains a variable called 'diseaseAge'!")
  }
  x_id     <- as.numeric(names(t_init))
  data_id  <- data[[id_var]]
  data_age <- data[[time_var]]
  uid      <- unique(x_id)
  dage     <- rep(NaN, length(data_id))
  for(id in uid){
    i_data <- which(data_id == id)
    i_new <- which(x_id == id)
    dage[i_data] <- data_age[i_data] - t_init[i_new]
  }
  data$diseaseAge <-dage
  return(data)
}
