

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
  cn  <- paste("f_",1:length(cn),"(",cn,")", sep="")
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
get_onset_times <- function(id, age, disAge){
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


#' Extract samples of T_onset
#'
#' @param fit an object of class \code{lgpfit}
#' @return a matrix
extract_t_onset_samples <- function(fit){
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  if(fit@model@stan_dat$UNCRT==0){
    stop("uncertainty was not modeled!")
  }else{
    SMP <- rstan::extract(fit@stan_fit, pars = "T_onset")$T_onset[,1,]
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
  id  <- as.numeric(sd$X_id)
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

