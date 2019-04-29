
#' An S4 class to represent an lgp model
#'
#' @slot data The original unmodified data frame.
#' @slot stan_dat The data to be given as input to \code{rstan::sampling}.
#' @slot scalings Preprocessing scaling functions and their inverse operations.
#' @slot info Model info.
#'
lgpmodel <- setClass(
  # Set the name for the class
  "lgpmodel",
  
  # Define the slots
  slots = c(
    data     = "data.frame",
    stan_dat = "list",
    scalings = "list",
    info     = "list"
  )
)


#' Show a summary of an \code{lgpmodel}
#'
#' @export
#' @param object an object of class \code{lgpmodel}
#' @rdname show_lgpmodel
#' @return nothing
setMethod(f = "show",
          signature = "lgpmodel",
          definition = function(object)
          {
            print_prior <- TRUE
            str1 <- likelihood_as_str(object@stan_dat$LH)
            if(object@stan_dat$UNCRT==1){
              str2 <- "uncertain"
            }else{
              str2 <- "exact"
            }
            if(object@stan_dat$HMGNS==0){
              str3 <- "heterogeneous"
            }else{
              str3 <- "homogeneous"
            }
            D <- object@stan_dat$D
            cat("\n ---------- LGPMODEL SUMMARY ----------\n\n")
            str_dis <- paste("  Disease component \n",
                             "    - Onset type: ", str2, "\n",
                             "    - Effect type: ", str3, "\n", sep="")
            
            minfo <- model_info(object, print = FALSE)
            cat(minfo)
            cat(paste("  Likelihood: ", str1, "\n", sep=""))
            if(D[3]==1){
              cat(str_dis)
            }
            cat("\n")
            if(print_prior){
              prior_info <- prior_stan_to_readable(object@stan_dat)
              cat(prior_info)
            }
          }
)


#' An S4 class to represent the output of the \code{lgp_fit} function
#' 
#' @description All slots that are lists contain fields 'samples' and 'average'.
#' @slot stan_fit The \code{stanfit} object returned by \code{rstan::sampling}.
#' @slot model The \code{lgpmodel} object returned by \code{lgp_model}.
#' @slot components Inferred components.
#' @slot components_corrected Covariate-effect corrected components.
#' @slot component_relevances Inferred component relevances.
#' @slot covariate_relevances Inferred covariate relevances.
#' @slot covariate_selection Covariate selection info.
#' @slot signal_variance Signal variance.
#' @slot residual_variance Residual variance.
#' @slot Rhat Split Rhat statistics.
#'
lgpfit <- setClass(
  # Set the name for the class
  "lgpfit",
  
  # Define the slots
  slots = c(
    stan_fit             = "stanfit",
    model                = "lgpmodel",
    components           = "data.frame",
    components_corrected = "data.frame",
    component_relevances = "list",
    covariate_relevances = "list",
    signal_variance      = "numeric",
    residual_variance    = "numeric",
    Rhat                 = "numeric",
    covariate_selection  = "list"
  )
)


#' Show a summary of results of the \code{lgp} function
#'
#' @export
#' @param object an object of class \code{lgpfit}
#' @rdname show_lgpfit
#' @return nothing
setMethod(f = "show",
          signature = "lgpfit",
          definition = function(object)
          {
            cat("\n ---------- LGPFIT SUMMARY ----------\n\n")
            n_chains <- length(object@stan_fit@inits)
            LP       <- rstan::get_logposterior(object@stan_fit)
            runtime  <- get_runtime(object)
            n_iters  <- length(LP[[1]])
            
            cat("* Chains: ", n_chains, sep="")
            cat(", iterations: ", n_iters, "\n", sep="")
            cat("* Average runtime per chain: ",
                runtime$warmup, " s (warmup) and ",
                runtime$sampling, " s (sampling)\n", sep="")
            
            assess_convergence(object)
            sv <- object@signal_variance
            rv <- object@residual_variance
            pevf <- mean(sv/(sv+rv))
            
            cat("* Proportion of variance explained by signal = ", round(pevf, 3), "\n", sep="")
            cat("* Covariate relevances:\n")
            cat("\n")
            
            rel          <- object@covariate_relevances$average
            DF           <- round(rbind(rel),3)
            rownames(DF) <- c("PEV")
            print(DF)
            cat("\n")
          }
)


#' Visualize a fitted `lgpfit` object
#'
#' @export
#' @param fit an object of class \code{lgpfit}
#' @param color_scheme bayesplot color scheme
#' @param x does nothing
#' @param y does nothing
#' @rdname plot_lgpfit
#' @return a ggplot object
setMethod(f = "plot",
          signature = "lgpfit",
          definition = function(fit, x = 1, y = 1, color_scheme = "red")
          {
            h <- plot_relevances(fit, color_scheme)
            return(h)
          }
)
