#' Generate an artificial longitudinal data set
#'
#' @export
#' @description Generate an artificial longitudinal data set.
#' @inheritParams create_X
#' @inheritParams create_F
#' @inheritParams create_y
#' @param N_affected Number of diseased individuals that are affected by the disease. This defaults
#' to the number of diseased individuals. This argument can only be given if \code{covariates}
#' contains a zero.
#' @param t_observed Determines how the disease onset is observed. This can be any function 
#' that takes the real disease onset as an argument and returns the (possibly randomly generated)
#' observed onset time. Alternatively, this can be a string of the form \code{"after_n"}
#' or \code{"random_p"}.
#' @param f_var variance of f
#' @param C_hat A constant added to f
#' @return A list \code{out}, where 
#' \itemize{
#'   \item \code{out$data} is a data frame containing the actual data and 
#'   \item \code{out$components} contains more points for smoother visualizations of the 
#'   generating process. 
#'   \item \code{out$onsets} contains the real disease onset times
#'   \item \code{out$p_signal} proportion of variance explained by signal
#' }
#' @examples 
#' # Generate Gaussian data
#' dat <- simulate_data(N = 4, t_data = c(6,12,24,36,48), snr = 3)
#' 
#' # Generate negative binomially distributed count data
#' dat <- simulate_data(N = 6, t_data = seq(2, 10, by = 2), noise_type = "NB", phi = 2)
simulate_data <- function(N,
                          t_data,
                          covariates      = c(),
                          names           = NULL,
                          relevances      = c(1,1, rep(1, length(covariates)) ),
                          n_categs        = rep(2, sum(covariates %in% c(2,3))),
                          t_jitter        = 0,
                          lengthscales    = rep(12, 2 + sum(covariates %in% c(0,1,2))),
                          f_var           = 1,
                          noise_type      = "Gaussian",
                          snr             = 3,
                          phi             = 1,
                          N_affected      = round(N/2),
                          onset_range     = range(t_data),
                          t_observed      = "after_0",
                          C_hat           = 0,
                          dis_fun         = NULL,
                          continuous_info = list(mu = c(pi/8,pi,-0.5), lambda = c(pi/8,pi,1))
)
{
  if(N < 2) stop("There must be at least 2 individuals!")
  if(length(t_data) < 3) stop("There must be at least 3 time points per individual!")
  
  # Input checks
  names <- sim_check_covariates(covariates, relevances, names, n_categs)
  if(N_affected > round(N/2)) stop("N_affected cannot be greater than round(N/2)!")
  if(is.unsorted(covariates)){
    stop("The covariates vector must be increasing!")
  }
  
  # Generate the covariates
  IN <- create_X(N               = N, 
                 covariates      = covariates, 
                 names           = names, 
                 n_categs        = n_categs,
                 t_data          = t_data ,
                 t_jitter        = t_jitter,
                 onset_range     = onset_range,
                 continuous_info = continuous_info)
  X <- IN$X
  
  # Compute X_affected
  k <- length(t_data)
  X_affected <- c(rep(1, N_affected*k), rep(0, (N - N_affected)*k))
  
  # Simulate the components F
  FFF <- create_F(X,
                  covariates,
                  relevances,
                  lengthscales,
                  X_affected,
                  dis_fun
  )
  
  # Total signal f is a (scaled) sum of the components plus a constant
  f <- rowSums(FFF)
  f <- sqrt(f_var)/stats::sd(f) * f
  f <- f + C_hat
  
  # Generate noise
  NOISY <- create_y(noise_type, f, snr = snr, phi = phi)
  g     <- NOISY$g
  y     <- NOISY$y
  noise <- y - g
  
  # Create the output objects
  dat     <- cbind(X, y)
  rownames(dat) <- 1:(N*k)
  comp    <- cbind(FFF, f, g, noise, y)
  
  # Real diseaseAges to observed ones
  OBSERVED <- sim_data_to_observed(dat, t_observed)
  
  # Signal and noise sums of squares
  SSR <- sum((g-mean(g))^2)
  SSE <- sum(noise^2)
  
  # Return 
  out <- list(data            = OBSERVED$dat,
              components      = comp,
              onsets          = IN$onsets,
              onsets_observed = OBSERVED$onsets_observed,
              par_ABO         = IN$par_cont,
              par_ell         = lengthscales,
              p_signal        = SSR/(SSR + SSE))
  
  return(out)
}


#' Generate noisy observations

#' @param noise_type Either "Gaussian", "Poisson" or "NB" (negative binomial).
#' @param snr The desired signal-to-noise ratio. This argument is valid
#' only with \cr
#' \code{noise_type = "Gaussian"}.
#' @param phi The dispersion parameter for negative binomial noise.
#' @param f The underlying signal.
#' @return A list \code{out}, where
#' \itemize{
#'   \item \code{out$g} is \code{f} mapped through an inverse link function and
#'   \item \code{out$y} is the noisy response variable.
#' }
create_y <- function(noise_type, f, snr, phi){
  L  <- length(f)
  y  <- rep(0, L)
  if(noise_type=="Gaussian"){
    sf      <- stats::var(f)
    sigma_n <- 1/sqrt(snr)*sqrt(sf)
    y_n     <- stats::rnorm(n = L, mean = 0, sd = 1)
    y_n     <- sigma_n*(y_n - mean(y_n))/stats::sd(y_n)
    g       <- f
    y       <- g + y_n
  }else if(noise_type=="Poisson"){
    g <- exp(f)
    for(i in 1:L){
      y[i] <- stats::rpois(n = 1, lambda = g[i])
    }
  }else if(noise_type=="NB"){
    g <- exp(f)
    for(i in 1:L){
      y[i] <- stats::rnbinom(n = 1, mu = g[i], size = phi)
    }
  }else if(noise_type=="Bernoulli"){
    g <- 1/(1 + exp(-f))
    for(i in 1:L){
      y[i] <- stats::rbinom(n = 1, size = 1, prob = g[i])
    }
  }else{
    stop(paste("Invalid input noise_type = '", noise_type, "'", sep =""))
  }
  NOISY <- list(g = g, y = y)
  return(NOISY)
}


#' Input check for the covariates-related arguments of \code{simulate_data}
#' @param covariates argument to \code{simulate_data}
#' @param relevances argument to \code{simulate_data}
#' @param names argument to \code{simulate_data}
#' @param n_cat the \code{n_categs} argument to \code{simulate_data}
#' @return the covariate names
sim_check_covariates <- function(covariates, relevances, names, n_cat){
  
  # Validity check
  d0 <- sum(covariates==0)
  if(d0 > 1) stop("There can be only one diseaseAge component!")
  L1 <- length(covariates)
  L2 <- length(relevances)
  L3 <- sum(covariates %in% c(2,3))
  if((L1+2) != L2) stop("length(relevances) must be length(covariates) + 2")
  if(L3 != length(n_cat)) stop("The argument n_cat has invalid length!")
  
  if(!is.null(names)){
    # Check input names
    L4 <- length(names)
    if(L4!=L1) stop("length(names) must be the same as length(covariates)")
    names <-c("id", "age", names)
  }else{
    # Generate default names
    names <- sim_generate_names(covariates)
  }
  
  if(sum(relevances<0)>0){
    stop("The relevances must be non-negative!")
  }
  
  return(names)
}


#' Generate names for covariates
#' @param covariates vector of covariate types
#' @return covariate names
sim_generate_names <- function(covariates){
  
  
  # Get default names
  names <- c("id", "age")
  def   <- c("x", "z", "offset")
  d0 <- sum(covariates==0)
  d1 <- sum(covariates==1)
  d2 <- sum(covariates==2)
  d3 <- sum(covariates==3)
  if(d0==1){
    names <- c(names, "diseaseAge" )
  }
  if(d1 > 0){
    if(d1 > 1) {
      names <- c(names, paste(def[1], 1:d1, sep="") )
    }else{
      names <- c(names, def[1])
    }
  }
  if(d2 > 0){
    if(d2 > 1) {
      names <- c(names, paste(def[2], 1:d2, sep="") )
    }else{
      names <- c(names, def[2])
    }
  }
  if(d3 > 0){
    if(d3 > 1) {
      names <- c(names, paste(def[3], 1:d3, sep="") )
    }else{
      names <- c(names, def[3])
    }
  }
  
  return(names)
}


#' Real generated disease ages to observed ones
#' @param dat data frame
#' @param t_observed Determines how the disease onset is observed. See documentation
#' of \code{\link{simulate_data}}.
#' @return a new data frame and observed onsets
sim_data_to_observed <- function(dat, t_observed){
  flag <- !("diseaseAge" %in% colnames(dat))
  id     <- dat$id
  uid    <- unique(id)
  N      <- length(uid)
  onsets_observed <- rep(NaN, N)
  if(flag){
    # No modifications to data frame needed
    ret <- list(dat = dat, onsets_observed = onsets_observed)
    return(ret)
  }else{
    age    <- dat$age
    disAge <- dat$diseaseAge
    j      <- 0 
    for(ID in uid){
      j     <- j + 1
      inds  <- which(id==ID)
      age_i <- age[inds]
      dag_i <- disAge[inds]
      if(is.nan(dag_i[1])){
        # not a diseased individual
      }else{
        
        # how many points are there after the real onset?
        irem <- which(dag_i > 0) 
        rem  <- age_i[irem]
        t0_real <- -dag_i[1] + age_i[1]
        print(t0_real)
        
        # sample the observed onset t0
        if(is.function(t_observed)){
          t_possible <- t_observed(t0_real)
          inds0 <- which(rem > t_possible)
          if(length(inds0)<1){
            stop("There are no data points after t = ", t_possible,"!")
          }
          idx0  <- inds0[1]
        }else{
          parsed <- sim_parse_t_obs(t_observed)
          
          if(parsed$type=="random"){
            idx0 <- rtgeom(length(irem), parsed$value)
          }else{
            idx0 <- parsed$value + 1
          }
          if(idx0 > length(rem)){
            stop("Not enough data points to go that far!")
          }
        }
        t0   <- rem[idx0]
        onsets_observed[j] <- t0
        
        # update the diseaseAge covariate
        dat$diseaseAge[inds] <- age_i - t0
      }
      
    }
    ret <- list(dat=dat, onsets_observed = onsets_observed)
    return(ret)
  }
}

#' Sample from the 'truncated geometric' distribution
#' @param p a number between 0 and 1
#' @param s an integer
#' @param n number of samples
#' @return an integer from the interval 1...n
rtgeom <- function(s, p, n = 1){
  r <- sample.int(n = s, size = n, prob = p^(0:(s-1)), replace = TRUE)
  return(r)
}


#' Parse the t_observed argument of \code{simulate_data}
#' @param t_observed a string
#' @return a list with a name and number
sim_parse_t_obs <- function(t_observed){
  
  parts <- strsplit(t_observed, "_")[[1]]
  if(parts[1]=="after"){
    type <- "after"
  }else if(parts[1]=="random"){
    type <- "random"
  }else{
    stop("unknown keyword ", parts[1])
  }
  value <- as.numeric(parts[2])
  return(list(type=type, value=value))
}
