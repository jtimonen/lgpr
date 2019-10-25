#' Get priors as a format that can be input to Stan
#'
#' @param D an integer vector of length 6
#' @param prior The \code{prior} argument supplied to \code{lgp()}.
#' @param HMGNS Is diseaseAge assumed to have a homogenous effect (1) or not (0)?
#' @param UNCRT Boolean value, is uncertainty of disease onset modeled?
#' @param N_cases number of case individuals
#' @param T_observed observed disease onsets
#' @param T_last last time point for each diseased individual
#' @return a list with all things related to priors that Stan needs
prior_to_stan <- function(D, prior, HMGNS, UNCRT, N_cases, T_observed, T_last){
  
  if(length(D)!=6){
    stop("D must be an integer vector with length 6!")
  }
  
  # Prior to Stan-readable format
  ta1 <- parse_prior_distribution(prior$alpha$idAge)$typ
  pa1 <- parse_prior_distribution(prior$alpha$idAge)$par
  ta2 <- parse_prior_distribution(prior$alpha$sharedAge)$typ
  pa2 <- parse_prior_distribution(prior$alpha$sharedAge)$par
  ta3 <- parse_prior_distribution(prior$alpha$diseaseAge)$typ
  pa3 <- parse_prior_distribution(prior$alpha$diseaseAge)$par
  ta4 <- parse_prior_distribution(prior$alpha$continuous)$typ
  pa4 <- parse_prior_distribution(prior$alpha$continuous)$par
  ta5 <- parse_prior_distribution(prior$alpha$categorical)$typ
  pa5 <- parse_prior_distribution(prior$alpha$categorical)$par
  ta6 <- parse_prior_distribution(prior$alpha$offset)$typ
  pa6 <- parse_prior_distribution(prior$alpha$offset)$par
  
  tl1 <- parse_prior_distribution(prior$lengthscale$idAge)$typ
  pl1 <- parse_prior_distribution(prior$lengthscale$idAge)$par
  tl2 <- parse_prior_distribution(prior$lengthscale$sharedAge)$typ
  pl2 <- parse_prior_distribution(prior$lengthscale$sharedAge)$par
  tl3 <- parse_prior_distribution(prior$lengthscale$diseaseAge)$typ
  pl3 <- parse_prior_distribution(prior$lengthscale$diseaseAge)$par
  tl4 <- parse_prior_distribution(prior$lengthscale$continuous)$typ
  pl4 <- parse_prior_distribution(prior$lengthscale$continuous)$par
  tl5 <- parse_prior_distribution(prior$lengthscale$categorical)$typ
  pl5 <- parse_prior_distribution(prior$lengthscale$categorical)$par
  
  tst  <- parse_prior_distribution(prior$warp_steepness)$typ
  pst  <- parse_prior_distribution(prior$warp_steepness)$par
  tsig <- parse_prior_distribution(prior$sigma_n)$typ
  psig <- parse_prior_distribution(prior$sigma_n)$par
  tphi <- parse_prior_distribution(prior$phi)$typ
  pphi <- parse_prior_distribution(prior$phi)$par
  
  # Store integer encodings of prior types
  TYP <- list(  c(ta1,tl1),
                c(ta2,tl2),
                c(ta3,tl3,tst),
                c(ta4,tl4),
                c(ta5,tl5),
                c(ta6),
                c(tsig), # sigma_n
                c(tphi)  # phi
  )
  
  # Store prior hyperparameters
  PAR <- list(  c(pa1,pl1),
                c(pa2,pl2),
                c(pa3,pl3,pst),
                c(pa4,pl4),
                c(pa5,pl5),
                c(pa6),
                c(psig),
                c(pphi)
  )
  
  # A deprecated feature
  HS <- rep(0,6)
  
  # Kernel hyperprior types
  t_ID  <- repvec(TYP[[1]], D[1])
  t_A   <- repvec(TYP[[2]], D[2])
  t_D   <- repvec(TYP[[3]], D[3])
  t_CNT <- repvec(TYP[[4]], D[4])
  t_CAT <- repvec(TYP[[5]], D[5])
  t_OFS <- repvec(TYP[[6]], D[6])
  t_SIG <- TYP[[7]]
  t_PHI <- TYP[[8]]
  
  # Kernel hyperprior params
  p_ID  <- repvec(PAR[[1]], D[1])
  p_A   <- repvec(PAR[[2]], D[2])
  p_D   <- repvec(PAR[[3]], D[3])
  p_CNT <- repvec(PAR[[4]], D[4])
  p_CAT <- repvec(PAR[[5]], D[5])
  p_OFS <- repvec(PAR[[6]], D[6])
  p_SIG <- PAR[[7]]
  p_PHI <- PAR[[8]]
  p_BET <- c(prior$beta$shape1, prior$beta$shape2)
  
  # Parse prior of uncertain disease effect time
  ONSET <- parse_prior_onset(prior$onset, N_cases, T_observed, T_last, UNCRT)
  
  # Things returned and given to Stan
  stan_priors <- list(t_ID  = t_ID,    t_A   = t_A,   
                      t_D   = t_D,     t_CNT = t_CNT,
                      t_CAT = t_CAT,   t_OFS = t_OFS,
                      t_SIG = t_SIG,   t_PHI = t_PHI,
                      
                      p_ID  = p_ID,    p_A   = p_A,   
                      p_D   = p_D,     p_CNT = p_CNT,
                      p_CAT = p_CAT,   p_OFS = p_OFS,
                      p_SIG = p_SIG,   p_PHI = p_PHI,
                      HS    = HS,      p_BET = p_BET,
                      
                      t_ONS = ONSET$t_ONS,
                      p_ONS = ONSET$p_ONS, 
                      L_ons = ONSET$L_ons,
                      U_ons = ONSET$U_ons,
                      
                      BACKWARDS = ONSET$backwards,
                      RELATIVE  = ONSET$relative,
                      vm_params = prior$vm_params
  )
  
  return(stan_priors)
}


#' Turn a list describing an onset prior distribution into things to be given to Stan
#'
#' @param dist This is \code{prior$onset}, where \code{prior} is an argument of
#' \code{lgp_model}
#' @param N_cases number of case individuals
#' @param T_observed observed disease onsets
#' @param T_last last time point for each diseased individual
#' @param UNCRT 0 or 1
#' @return a list with things to be given to Stan
parse_prior_onset <- function(dist, N_cases, T_observed, T_last, UNCRT){
  
  # Check type
  if(is.null(dist$type)){
    stop("The onset prior must have a field 'type'!")
  }
  type <- dist$type
  
  # Parse type
  if(is.character(type)){
    parts <- strsplit(type,"_")[[1]]
    backwards <- FALSE
    relative <- FALSE
    if(length(parts)>2){
      if(parts[3]=="backwards"){
        backwards <- TRUE
      }else{
        stop("Unknown third keyword '", parts[3],"'! ")
      }
    }
    if(length(parts) > 1){
      if(parts[2]=="whole"){
        L_ons  <- rep(0, N_cases)
        U_ons  <- T_last 
      }else if(parts[2]=="before"){
        L_ons  <- rep(0, N_cases)
        U_ons  <- T_observed 
      }else if(parts[2]=="relative"){
        relative <- TRUE
        L_ons  <- rep(0, N_cases)
        U_ons  <- T_last
        if(length(parts)>2){
          stop("Do not set third keyword if using relative mode!")
        }
        if(!(parts[1] %in% c("normal", "student-t"))){
          stop("Only normal and student-t priors available in relative mode!")
        }
      }else{
        stop("Unknown second keyword '", parts[2],"'! ")
      }
    }else{
      L_ons  <- rep(0, N_cases)
      U_ons  <- T_last 
    }
    dist$type <- parts[1]
    tons   <- parse_prior_distribution(dist)$typ
    pons   <- parse_prior_distribution(dist)$par
    t_ONS  <- repvec(tons, N_cases)
    p_ONS  <- repvec(pons, N_cases)
  }else{
    stop("The field 'type' must be a string!")
  }
  
  L_ons <- repvec(L_ons, UNCRT)
  U_ons <- repvec(U_ons, UNCRT)
  
  # Return a list
  ret <- list(t_ONS     = t_ONS,
              p_ONS     = p_ONS,
              L_ons     = L_ons,
              U_ons     = U_ons,
              backwards = as.numeric(backwards),
              relative  = as.numeric(relative))
  return(ret)
}

#' Turn a list describing a prior distribution into vectors to be given to Stan
#'
#' @param dist a list with field type, and possibly others
#' @param add_correct additional correct parameter names
#' @return a list with two vectors to be given to Stan
parse_prior_distribution <- function(dist, add_correct = NULL){
  if(is.null(dist)){
    stop("dist cannot be NULL!")
  }
  types  <- c(get_prior_type(dist$type), get_transform_type(dist$transform))
  params <- c(get_prior_params(dist, add_correct), 1) # trailing one because of a deprecated feature
  if(length(params)!=3)        stop("Length of params is not 3!");
  if(length(types)!=2)         stop("Length of types is not 2!")
  if(sum(is.null(types)) > 0)  stop("types contains NULL values")
  if(sum(is.null(params)) > 0) stop("params contains NULL values")
  ret <- list(typ = types, par = params)
  return(ret)
}


#' A dictionary from distribution names to integer encoding
#'
#' @param type type of the distribution as a string
#' @return an integer
get_prior_type <- function(type){
  if(is.null(type)){
    stop("Prior type must be given!")
  }
  TYPE <- switch(
    type,
    "uniform"=1,
    "normal"=2,
    "student-t"=3,
    "gamma"=4,
    "inv-gamma"=5,
    "log-normal"=6
  )
  if(is.null(TYPE)){
    msg <- paste("Invalid prior distribution type '", type, "'.", sep="")
    msg <- paste(msg, " Valid types are 'uniform', 'normal', 'student-t', 'gamma',",
                 "'inv-gamma' and 'log-normal'.", sep="")
    stop(msg)
  }else{
    return(TYPE)
  }
}


#' A dictionary from transform names to integer encoding
#'
#' @param type Type of the transform as a string. Allowed arguments
#' are "none" or "square". If NULL, "none" is used.
#' @return an integer (0, 1 or 2)
get_transform_type <- function(type){
  if(is.null(type)){
    type <- "none"
  }
  TYPE <- switch(
    type,
    "none"=0,
    "square"=1
  )
  if(is.null(TYPE)){
    stop(paste("Invalid transform type '", type, "'.", sep=""))
  }else{
    return(TYPE)
  }
}


#' Get prior parameters
#'
#' @param dist the distribution
#' @param add_correct additional correct parameter names
#' @return a hyperparameter vector of length 2
get_prior_params <- function(dist, add_correct){
  type <- dist$type
  if(is.null(type)) stop("type must be defined!")
  if(type=="uniform"){
    l <-c()
    v <- c(0,0)
  }else if(type=="normal" || type=="log-normal"){
    l <- c("mu", "sigma")
    v <- c(dist$mu, dist$sigma)
  }else if(type=="student-t"){
    l <- c("mu", "sigma","nu")
    v <- c(dist$nu, dist$sigma)
    if(dist$mu !=0){
      stop("for the student-t prior, hyperparameter 'mu' must be zero!")
    }
  }else if(type=="gamma"){
    l <- c("shape", "rate")
    v <- c(dist$shape, dist$rate)
  }else if(type=="inv-gamma"){
    l <- c("shape", "scale")
    v <- c(dist$shape, dist$scale)
  }else{
    stop("Invalid distribution type", type, "!")
  }
  check_hyperparameter_names(dist, correct = c(l, add_correct))
  return(v)
}


#' An error message for wrong hyperparameter naming
#'
#' @param dist the distribution
#' @param correct the allowed hyperparameter names
#' @return nothing
check_hyperparameter_names <- function(dist, correct){
  type <- dist$type
  if(is.null(type)) stop("distribution type cannot be NULL!")
  nam  <- names(dist)
  nam  <- setdiff(nam, c("transform", "type"))
  
  # Check missing parameters
  mis  <- setdiff(correct, nam)
  if(length(mis) > 0){
    msg <- paste("The ", type, " distribution is missing the parameter '",
                 mis[1], "'. The needed parameters are:", sep="")
    msg <- paste(msg, paste(correct, collapse = ", "))
    stop(msg)
  }
  
  # Check extra parameters
  bad  <- setdiff(nam, correct)
  if(length(bad) > 0){
    if(type=="uniform"){
      msg <- paste("The uniform distribution does not have a parameter called '",
                   bad[1], "'. There are not any allowed parameters for the uniform distribution.")
    }else{
      msg <- paste("The ", type, " distribution does not have a parameter called '",
                   bad[1], "'. The allowed parameters are:", sep="")
      msg <- paste(msg, paste(correct, collapse = ", "))
    }
    stop(msg)
  }
}

