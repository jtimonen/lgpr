#' Create the default prior
#' @export
#' @param sigma_alpha Sigma parameter of the student-t distribution for all alpha. 
#' @return A list defining a valid \code{prior} argument for the \code{lgp} function.
prior_default <- function(sigma_alpha = 1)
{
  
  # Magnitude parameters  
  pmag <- list(type      = "student-t",
               mu        = 0,
               sigma     = sigma_alpha,
               nu        = 20)
  
  
  prior_mag <- list(idAge       = pmag,
                    sharedAge   = pmag,
                    diseaseAge  = pmag,
                    continuous  = pmag,
                    categorical = pmag,
                    offset      = pmag)
  
  # Lengthscale parameters
  plen <- list(type      = "log-normal", 
               mu        = 0, 
               sigma     = (log(1)-log(0.1))/2)
  
  prior_ls <- list(idAge       = plen,
                   sharedAge   = plen,
                   diseaseAge  = plen,
                   continuous  = plen,
                   categorical = plen)
  
  
  # Other parameters
  prior_st  <- list(type = "inv-gamma",  shape  = 9,   scale = 5)
  prior_sig <- list(type = "log-normal", mu = 0, sigma = 1, transform = "square")
  prior_phi <- list(type = "log-normal", mu = 1, sigma = 1, transform = "square")
  prior_bet <- list(shape1 = 0.2,        shape2 = 0.2)
  
  # Uncertain disease onset
  prior_ons <- list(type = "uniform_whole")
  
  prior <- list(alpha          = prior_mag,
                lengthscale    = prior_ls,
                warp_steepness = prior_st,
                sigma_n        = prior_sig,
                phi            = prior_phi,
                beta           = prior_bet,
                onset          = prior_ons,
                vm_params      = c(0.01, 1))
  
  return(prior)
}


#' Create a similar default prior as in \code{LonGP} (Cheng et. al, 2019)
#' @export
#' @description Not recommended, because a lengthscale close to 0 is possible.
#' @return A list defining a valid \code{prior} argument for the \code{lgp_model} function.
prior_LonGP <- function()
{
  
  # Start by creating the default prior
  prior <- prior_default()
  

  t4    <- list(type = "student-t", mu = 0, sigma  = 1, nu = 4)
  p_sig <- list(type = "inv-gamma", shape  = 1/2, scale = 0.01/2, transform = "square")
  
  # Edit some parts of the default prior
  prior$alpha$offset <- t4
  prior$lengthscale$idAge <- t4 # these are dangerous, lengthscale ~ 0 is possible
  prior$lengthscale$categorical <- t4 # these are dangerous, lengthscale ~ 0 is possible
  prior$sigma_n <- p_sig
  return(prior)
}


#' Validate prior by sampling the signal and noise from it
#'
#' @export
#' @param model An object of class \code{\link{lgpmodel}}.
#' @param chains how many chains are used to sample from the prior
#' @param iter for how many iterations are the chains run
#' @param parallel should the chains be run in parallel?
#' @return An object of class \code{lgpfit} and random samples of both `f` and `y`.
#'
validate_prior <- function(model,
                           chains   = 4,
                           iter     = 1000,
                           parallel = FALSE)
{
  
  stop("sorry, this function does not currenly work!")
  # Remove likelihood
  model@stan_dat$LH <- 0
  
  # Sample from the prior
  fit <- lgp_fit(model, 
                 iter     = iter, 
                 parallel = parallel, 
                 chains   = chains,
                 refresh  = 2*iter)
  
  sf    <- fit@stan_fit
  F_rng <- rstan::extract(sf, pars = "F_rng")$F_rng[,1,]
  sig   <- rstan::extract(sf, pars = "sigma_n")$sigma_n
  phi   <- rstan::extract(sf, pars = "phi")$phi
  C_hat <- model@stan_dat$C_hat
  if(model@likelihood=="Gaussian"){
    fff <- (F_rng + C_hat)
    F_mean <- colMeans(fff)
    F_std  <- apply(fff,2,stats::sd)
    mean_f <- mean(F_mean)
    std_f  <- mean(F_std)
    mean_y <- mean_f
    std_y  <- mean(F_std + sig)
  }else if(model@likelihood == "Poisson"){
    fff    <- C_hat*exp(F_rng)
    F_mean <- colMeans(fff)
    F_std  <- apply(fff,2,stats::sd)
    mean_f <- mean(F_mean)
    std_f  <- mean(F_std)
    mean_y <- mean_f
    std_y  <- mean(sqrt(F_std^2 + F_mean))
  }else if(model@likelihood == "Negative Binomial"){
    fff    <- C_hat*exp(F_rng)
    F_mean <- colMeans(fff)
    n      <- length(F_mean)
    F_std  <- apply(fff,2,stats::sd)
    mean_f <- mean(F_mean)
    std_f  <- mean(F_std)
    S      <- length(phi)
    yyy    <- matrix(0, S, n)
    for(s in 1:S){
      yyy[s,] <- MASS::rnegbin(n = n, mu = fff[s,], theta = phi[s])
    }
    Y_mean <- colMeans(yyy)
    Y_std  <- apply(yyy,2,stats::sd)
    mean_y <- mean(Y_mean)
    std_y  <- mean(Y_std)
  }
  #if(noiseType == "Gaussian"){
  #  y_n <- matrix(0, s, n)
  #  for(is in 1:s){
  #    y_n[is,] <- stats::rnorm(n = n, mean = 0, sd = sig[is])
  #  }
  #  mu_rng <- f_rng
  # y_rng <- mu_rng + y_n
  #}else if(noiseType == "Poisson"){
  #  mu_rng <- mean(data$y)*exp(f_rng)
  #  y_rng  <- matrix(0, s, n)
  ##  for(is in 1:s){
  #    y_rng[is,] <- stats::rpois(n = n, lambda = mu_rng[is,])
  #  }
  #}else if(noiseType == "NB"){
  #  mu_rng <- mean(data$y)*exp(f_rng)
  #  y_rng  <- matrix(0, s, n)
  #  for(is in 1:s){
  #    y_rng[is,] <- MASS::rnegbin(n = n, mu = mu_rng[is,], theta = phi[is])
  #  }
  #}
  
  yname <- strsplit(model@formula, "~")[[1]][1]
  yname <- substr(yname,1,nchar(yname)-1)
  y_data <- model@data[[yname]]
  cat("\n")
  cat("* Mean of ", yname, ": ", mean(y_data),  " (data)\n", sep="")
  cat("* Std of ", yname, ": ", stats::sd(y_data), " (data)\n", sep="")
  cat("* Mean of signal: ", mean_f,  " (sampled)\n", sep="")
  cat("* Std of signal: ", std_f, " (sampled)\n", sep="")
  cat("* Mean of y: ", mean_y,  " (sampled)\n", sep="")
  cat("* Std of y: ", std_y, " (sampled)\n", sep="")
  
  ret <- list(mean_f = mean_f,
              std_f  = std_f,
              mean_y = mean_y,
              std_y  = std_y,
              C_hat  = C_hat)
  cat("\n")
  return(ret)
}


#' Human-readable description of a specified prior
#' 
#' @export
#' @description Print human-readable info about the prior specification that was used or will be used
#' @param object An object of class \code{lgpfit} or a valid prior argument for the `lgp` function.
#' @return nothing
print_prior <- function(object)
{
  if(class(object)=="lgpmodel"){
    sd   <- object@stan_dat
    str  <- prior_stan_to_readable(sd)
  }else{
    D <- c(1,1,1, 1,1,1)
    pseudo_stan_data   <- prior_to_stan(D, object, 0)
    pseudo_stan_data$D <- D 
    pseudo_stan_data$LH <- 0 
    str  <- prior_stan_to_readable(pseudo_stan_data)
  }
  cat(str)
}

#' Human-readable information about the priors in the Stan data object
#'
#' @param stan_dat The list that is passed as data to \code{rstan::sampling}.
#' @return Info as a string.
prior_stan_to_readable <- function(stan_dat){
  HS   <- stan_dat$HS
  info <- " ---------- PRIOR SPECIFICATIONS ----------\n"
  
  if(sum(HS)>0){
    info <- paste(info, "* A HS prior is used.\n")
    return(info)
  }
  
  D  <- stan_dat$D
  dist <- c("Uniform","Normal","Student-t","Gamma","Inverse-Gamma","Log-Normal","Beta")
  
  info_mag <- ""
  info_ls  <- ""
  info_oth <- ""
  
  # ID
  if(D[1]>0){
    TYP <- stan_dat$t_ID
    P   <- stan_dat$p_ID
    for(j in 1:D[1]){
      parnames <- c("alpha_id","lengthscale_id")
      str1     <- prior_statement(parnames[1],TYP[j,1:2],P[j,1:3],dist)
      info_mag <- paste(info_mag, str1)
      str2     <- prior_statement(parnames[2],TYP[j,3:4],P[j,4:6],dist)
      info_ls  <- paste(info_ls, str2)
    }
  }
  
  # Age
  if(D[2]>0){
    TYP <- stan_dat$t_A
    P <- stan_dat$p_A
    for(j in 1:D[2]){
      parnames <- c("alpha_age","lengthscale_age")
      str      <- prior_statement(parnames[1],TYP[j,1:2],P[j,1:3],dist)
      info_mag <- paste(info_mag, str)
      str      <- prior_statement(parnames[2],TYP[j,3:4],P[j,4:6],dist)
      info_ls  <- paste(info_ls, str)
    }
  }
  
  # DiseaseAge
  if(D[3]>0){
    TYP <- stan_dat$t_D
    P <- stan_dat$p_D
    for(j in 1:D[3]){
      parnames <- c("alpha_diseaseAge","lengthscale_diseaseAge","warp_steepness")
      str      <- prior_statement(parnames[1],TYP[j,1:2],P[j,1:3],dist)
      info_mag <- paste(info_mag, str)
      str      <- prior_statement(parnames[2],TYP[j,3:4],P[j,4:6],dist)
      info_ls  <- paste(info_ls, str)
      str      <- prior_statement(parnames[3],TYP[j,5:6],P[j,7:9],dist)
      info_oth <- paste(info_oth, str)
    }
  }
  
  
  # Continuous
  if(D[4]>0){
    TYP <- stan_dat$t_CNT
    P <- stan_dat$p_CNT
    for(j in 1:D[4]){
      if(D[4]>1){
        parnames <- c(paste("alpha_continuous[",j,"]",sep=""),
                      paste("lengthscale_continuous[",j,"]",sep=""))
      }else{
        parnames <- c(paste("alpha_continuous",sep=""),
                      paste("lengthscale_continuous",sep=""))
      }
      
      str      <- prior_statement(parnames[1],TYP[j,1:2],P[j,1:3],dist)
      info_mag <- paste(info_mag, str)
      str      <- prior_statement(parnames[2],TYP[j,3:4],P[j,4:6],dist)
      info_ls  <- paste(info_ls, str)
    }
  }
  
  # Categorical
  if(D[5]>0){
    TYP <- stan_dat$t_CAT
    P <- stan_dat$p_CAT
    for(j in 1:D[5]){
      if(D[5]>1){
        parnames <- c(paste("alpha_categAge[",j,"]",sep=""),
                      paste("lengthscale_categAge[",j,"]",sep=""))
      }else{
        parnames <- c(paste("alpha_categAge",sep=""),
                      paste("lengthscale_categAge",sep=""))
      }
      str      <- prior_statement(parnames[1],TYP[j,1:2],P[j,1:3],dist)
      info_mag <- paste(info_mag, str)
      str      <- prior_statement(parnames[2],TYP[j,3:4],P[j,4:6],dist)
      info_ls  <- paste(info_ls, str)
    }
  }
  
  # Offset
  if(D[6]>0){
    TYP <- stan_dat$t_OFS
    P <- stan_dat$p_OFS
    for(j in 1:D[6]){
      if(D[6]>1){
        parnames <- c(paste("alpha_categOffset[",j,"]",sep="")) 
      }else{
        parnames <- c(paste("alpha_categOffset", sep=""))
      }
      str      <- prior_statement(parnames[1],TYP[j,1:2],P[j,1:3],dist)
      info_mag <- paste(info_mag, str)
    }
  }
  
  # Info about other parameters
  if(stan_dat$LH==0 || stan_dat$LH==1){
    # noise std
    TYP <- stan_dat$t_SIG
    P <- stan_dat$p_SIG
    for(j in 1:1){
      parnames <- c("sigma_n")
      str      <- prior_statement(parnames[1],TYP[1:2],P[1:3],dist)
      info_oth <- paste(info_oth, str)
    }
  }
  if(stan_dat$LH==0 || stan_dat$LH==3){
    # phi
    TYP <- stan_dat$t_PHI
    P   <- stan_dat$p_PHI
    for(j in 1:1){
      parnames <- c("phi")
      str      <- prior_statement(parnames[1],TYP[1:2],P[1:3],dist)
      info_oth <- paste(info_oth, str)
    }
  }
  if(stan_dat$HMGNS==0){
    # beta
    P        <- c(stan_dat$p_BET, 1)
    str      <- prior_statement("beta",c(7,0),P[1:3],dist)
    info_oth <- paste(info_oth, str)
  }
  if(stan_dat$UNCRT==1){
    LB <- stan_dat$L_ons
    UB <- stan_dat$U_ons
    NC <- stan_dat$N_cases
    if(NC > 0){
      info_oth <- paste(info_oth, "\n")
      for(k in 1:NC){
        TYP <- stan_dat$t_ONS[k,]
        P   <- stan_dat$p_ONS[k,]
        pname <- paste("T_onset[1,",k,"]",sep="")
        if(stan_dat$backwards==1){
          pn_base <- pname
          pname <- paste("T_obs[",k,"] - ", pname, sep="")
        }
        ps  <- prior_statement(pname, TYP[1:2], P[1:3], dist, FALSE)
        str <- paste("[lower=", LB[k], ", upper=", UB[k], "]", sep="")
        if(stan_dat$backwards==1){
          str <- paste(str, " (bound is for ", pn_base, ")", sep="")
        }
        ps <- paste(ps, str, sep = "")
        ps <- paste(ps, "\n")
        info_oth <- paste(info_oth, ps)
      }
    }
  }
  info <- cat(info, info_mag, info_ls, info_oth, sep= "\n")
  return(info)
  
}

#' Human-readable prior statement
#'
#' @param parname parameter name
#' @param TYP two integers
#' @param P three real numbers
#' @param row_change should a newline be last character?
#' @param dist list of distribution names
#' @return Sampling statement as a string.
prior_statement <- function(parname, TYP, P, dist, row_change = TRUE){
  
  P <- round(P,3)
  
  # Check if there is scaling
  if(P[3]!=1){
    parname <- paste("1/",P[3], " * ", parname, sep="")
  }
  
  # Check if there is a transform
  if(TYP[2]==1){
    parname <- paste("(",parname,")^2", sep="")
  }else if(TYP[2]==0){
    parname <- parname
  }else{
    stop("Invalid transform ", TYP[2]," for parameter '", parname,
         "'! Must be 0 (none) or 1 (squaring).")
  }
  
  # Get prior statement
  if(TYP[1] %in% c(2,4,5,6,7)){
    str <- paste(parname," ~ ", dist[TYP[1]], "(",P[1],",",P[2],")", sep="")
  }else if(TYP[1]==3){
    str <- paste(parname," ~ ", dist[TYP[1]], "(nu=",P[1],",","mu=0",",","sigma=", P[2], ")", sep="")
  }else{
    str <- paste(parname," ~ ", dist[TYP[1]], sep="")
  }
  
  # Possible newline
  if(row_change){
    str <- paste(" ", str, "\n", sep="")
  }
  
  # Return
  return(str)
}

