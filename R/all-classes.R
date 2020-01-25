
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
              str2 <- "fixed"
            }
            if(object@stan_dat$HMGNS==0){
              str3 <- "heterogeneous"
            }else{
              str3 <- "homogeneous"
            }
            D <- object@stan_dat$D
            cat("\n ---------- LGPMODEL SUMMARY ----------\n\n")
            str_dis <- paste("  Disease component \n",
                             "    - Effect time: ", str2, "\n",
                             "    - Effect type: ", str3, "\n", sep="")
            
            minfo <- model_info(object, print = FALSE)
            cat(minfo)
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
#' @slot stan_fit The \code{stanfit} object returned by \code{rstan::sampling}.
#' @slot model The \code{lgpmodel} object returned by \code{lgp_model}.
#' @slot relevances Inferred component relevances.
#' @slot selection Component selection info.
#' @slot pkg_version Package version number.
#' @slot diagnostics  A data frame with columns \code{c("Rhat", "Bulk_ESS", "Tail_ESS")}.
#'
lgpfit <- setClass(
  # Set the name for the class
  "lgpfit",
  
  # Define the slots
  slots = c(
    stan_fit             = "stanfit",
    model                = "lgpmodel",
    relevances           = "list",
    selection            = "list",
    pkg_version          = "character",
    diagnostics          = "data.frame"
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
            
            # Runtime info
            tryCatch({
              n_chains <- length(object@stan_fit@inits)
              runtime  <- get_runtime(object)
              cat("* Average runtime per chain: ",
                  runtime$warmup, " s (warmup) and ",
                  runtime$sampling, " s (sampling)\n", sep="")
            }, error = function(e) {
              cat("* Unable to show runtime info. Reason:\n")
              print(e)
            })
            

            # Convergence info
            tryCatch({
              diag <- object@diagnostics
              Rhat <- diag$Rhat
              imax <- which(Rhat==max(Rhat))
              cat("* Largest R-hat value is ", round(max(Rhat), 4), 
                  " (", paste(rownames(diag)[imax], collapse = ', '), 
                  ")\n", sep="")
            }, error = function(e) {
              cat("* Unable to show convergence info. Reason:\n")
              print(e)
            })
            
            # Relevance info
            tryCatch({
              sel  <- object@selection
              rel_method <- object@relevances$method
              r3   <- sel$prob
              tr   <- sel$threshold
              cat("* Used relevance method = ", rel_method, "\n", sep ="")
              cat("* Used selection threshold = ", tr, "\n", sep ="")
              
              rel  <- object@relevances$average
              cn   <- names(rel)
              r1   <- round(rel,3)
              r2   <- cn %in% sel$selected
              DF   <- data.frame(r1,r2,r3)
              cat("\n")
              
              colnames(DF) <- c("Relevance", "Selected", "Prob." )
              print(DF)
              cat("\n")
            }, error = function(e) {
              cat("* Unable to show component relevance info. Reason:\n")
              print(e)
            })
            

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
            h <- plot_relevances(fit, color_scheme = color_scheme)
            return(h)
          }
)
