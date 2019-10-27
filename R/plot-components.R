
#' Visualize average inferred components or component samples,
#' evaluated at data points
#' @export
#' @param fit An object of class \code{lgpfit}.
#' @param time_is_xvar is age the x-axis variable?
#' @param subsamples How many samples to plot. If this is NULL,
#' average over all samples is plotted. If this is "all", all
#' samples are plotted.
#' @param ... additional arguments for \code{\link{plot_components}}
#' @return an object returned by ggpubr::ggarrange
plot_components_posterior <- function(fit, subsamples = NULL,
                                      time_is_xvar = TRUE, ...)
{
  # Check input correctness
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  model <- fit@model
  
  # Create plot
  FFF <- get_function_component_samples(fit)
  ddd <- dim(FFF)[3]
  FFF <- FFF[,,1:(ddd-1)]
  
  if(is.null(subsamples)){
    FFF_1 <- apply(FFF, c(2,3), mean)
    nnn     <- dim(FFF_1)[1]
    FFF     <- array(0, dim = c(1, nnn, ddd-1))
    FFF[1,,]<- as.matrix(FFF_1)
  }else if(is.character(subsamples)){
    if(subsamples=="all"){
      # all samples
    }else{
      stop("invalid input for subsamples!")
    }
  }else if(is.numeric(subsamples)){
    n_smp   <- dim(FFF)[1]
    i_sub   <- sample.int(n_smp, subsamples)
    FFF     <- FFF[i_sub,,] 
  }else{
    stop("invalid input for subsamples!")
  }
  
  h <- plot_components(FFF, NULL, model, time_is_xvar, ...)
  
  # Return ggplot object
  return(h)
  
}


#' Visualize the components of a simulated data set
#' @export
#' @param simData simulated data object (list)
#' @param time_is_xvar is the time variable the x-axis variable?
#' @param ... additional arguments for \code{\link{plot_components}}
#' @return an object returned by ggpubr::ggarrange
plot_components_simdata <- function(simData, time_is_xvar = TRUE, ...)
{
  data  <- simData$data
  FFF   <- simData$components
  model <- full_model(data)
  
  # Create plot
  nnn     <- dim(FFF)[1]
  ddd     <- dim(FFF)[2]
  FFF_1   <- FFF[,1:(ddd-3)]
  FFF     <- array(0, dim = c(1, nnn, ddd-3))
  FFF[1,,]<- as.matrix(FFF_1)
  
  h <- plot_components(FFF, NULL, model, time_is_xvar, ...)
  return(h)
}

#' Helper function for plotting components
#' 
#' @inheritParams plot_component 
#' @param ncol number of plot columns
#' @param nrow number of plot rows
#' @param legend legend argument for ggarrange, use "none" to remove legends
#' @param labels labels argument for ggarrange
#' @param ylim y axis limits
#' @param font_size font size for plots
#' @param theme ggplot theme
#' @param legend_dir direction of legend
#' @param xlabel x-axis label
#' @param ylabel y-axis label 
#' @return an object returned by ggpubr::ggarrange
plot_components <- function(MMM, SSS, model, time_is_xvar,
                            linealpha   = 1,
                            linetype    = 1, 
                            ncol        = NULL,
                            nrow        = NULL,
                            legend      = NULL,
                            labels      = NULL,
                            ylim        = NULL,
                            font_size   = 9,
                            theme       = ggplot2::theme_linedraw(),
                            legend_dir  = "horizontal",
                            xlabel      = NULL,
                            ylabel      = " "){
  
  GG <- list()
  sum_D <- sum(model@stan_dat$D) + 1
  for(d in 1:sum_D){
    gg <- plot_component(MMM, SSS, model, d, time_is_xvar,
                         linealpha, linetype)
    if(is.null(ylim)){
      if(!is.null(SSS)){
        ylim <- c(min(MMM - SSS), max(MMM + SSS))
      }else{
        ylim <- range(MMM)
      }
    }
    gg <- gg + theme + ggplot2::ylim(ylim)
    if(!is.null(xlabel)){gg <- gg + ggplot2::xlab(xlabel)}
    if(!is.null(ylabel)){gg <- gg + ggplot2::ylab(ylabel)}
    if(is.null(legend)){
      gg <- gg + ggplot2::theme(
        legend.justification = c(0.95,0.05), 
        legend.position = c(0.95,0.05),
        legend.key = ggplot2::element_rect(fill = "gray95"),
        legend.background = ggplot2::element_rect(fill = "gray95"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        legend.direction = legend_dir
      )
    }
    GG[[d]] <- gg + ggplot2::theme(text=ggplot2::element_text(size=font_size))
  }
  p <- ggpubr::ggarrange(plotlist = GG, ncol = ncol, nrow = nrow,
                         legend = legend,
                         labels = labels)
  return(p)
}

#' Helper function for plotting one component
#' @param MMM a n array of size n_samples x n_data x n_components
#' @param SSS a n array of size n_samples x n_data x n_components
#' @param model an object of class 'lgpmodel'
#' @param idx Index of component to be plotted.
#' @param time_is_xvar boolean value 
#' @param linealpha line alpha
#' @param linetype line type
#' @return a ggplot object
plot_component <- function(MMM, SSS, model, idx, time_is_xvar,
                           linealpha, linetype)
{
  sdat  <- model@stan_dat
  X     <- t(sdat$X)
  Xnn   <- sdat$X_notnan
  D     <- sdat$D
  SCL   <- model@scalings
  sum_D <- sum(D) + 1
  S     <- dim(MMM)[1]
  n     <- dim(MMM)[2]
  d     <- dim(MMM)[3]
  cpn   <- model@info$component_names
  cvn   <- model@info$covariate_names
  if(sum_D!=d){ stop("dim(MMM)[3] must be ", sum_D) }
  if(idx < 1 || idx > sum_D){ stop("idx must be between 1 and ", sum_D) }
  if(idx < sum_D){
    ctype <- component_index_to_type(D, idx)
    cind  <- component_index_to_covariate_index(D, idx)
  }else{
    ctype <- 7 # sum 
  }
  
  plot_eb  <- !is.null(SSS)
  MMM_i    <- MMM[,,idx]
  f        <- as.numeric(t(MMM_i))
  if(plot_eb){
    UUU_i <- MMM_i + SSS[,,idx]
    LLL_i <- MMM_i - SSS[,,idx]
    ub    <- as.numeric(t(UUU_i))
    lb    <- as.numeric(t(LLL_i))
  }
  sample <- rep(1:S, each = n)
  id     <- rep(X[,1], S)
  
  if(ctype %in% c(1,2,5,6,7)){
    time_is_xvar <- TRUE
  }
  if(time_is_xvar){
    xvar  <- rep(SCL$TSCL$fun_inv(X[,2]), S)
    xvar_name <- model@info$varInfo$time_variable
  }else{
    xvar  <- rep(X[,cind], S)
    if(ctype==4){
      xvar_name <- cvn[cind]
      warning('covariate', cvn[cind], 'is not on original scale!')
    }
  }
  
  # Create df and aes
  if(ctype %in% c(1,4,7)){
    leg      <- ggplot2::theme(legend.position = "none")
    grpvar   <- as.factor(paste(id, sample))
    df       <- data.frame(xvar, f, grpvar)
    aes      <- ggplot2::aes_string(x = 'xvar', y = 'f', group = 'grpvar')
  }else if(ctype==2){
    leg      <- ggplot2::theme(legend.position = "none")
    grpvar   <- as.factor(sample)
    df       <- data.frame(xvar, f, grpvar)
    aes      <- ggplot2::aes_string(x = 'xvar', y = 'f', group = 'grpvar')
  }else if(ctype==3){
    colorvar <- as.factor(rep(Xnn, S))
    leg      <- ggplot2::labs(color = 'Group')
    grpvar   <- as.factor(paste(id, sample))
    df       <- data.frame(xvar, f, grpvar, colorvar)
    aes      <- ggplot2::aes_string(x = 'xvar', y = 'f', 
                                    group = 'grpvar', color = 'colorvar')
  }else if(ctype == 5 || ctype == 6){
    leg      <- ggplot2::labs(color = cvn[cind])
    colorvar <- rep(X[,cind], S)
    grpvar   <- as.factor(paste(colorvar, sample))
    colorvar <- as.factor(colorvar)
    df       <- data.frame(xvar, f, grpvar, colorvar)
    aes      <- ggplot2::aes_string(x = 'xvar', y = 'f', 
                                    group = 'grpvar', color = 'colorvar')
  }else{
    stop("invalid ctype!")
  }
  
  # Create title
  if(ctype!=7){
    title <- parse(text = cpn[idx])
  }else{
    title <- parse(text = "f[sum]")
  }
  
  # Create ggplot object
  h   <- ggplot2::ggplot(df, aes) + leg
  h   <- h + ggplot2::geom_line(linetype = linetype, alpha = linealpha)
  h   <- h + ggplot2::ggtitle(title)
  h   <- h + ggplot2::xlab(xvar_name)
  # Add errorbar
  if(plot_eb){
    h <- h + ggplot2::geom_ribbon(ggplot2::aes(ymin=lb, ymax=ub))
  }
  return(h)
  
}


#' Component index to component type
#' @param D integer vector of length 
#' @param idx integer
#' @return an integer
component_index_to_type <- function(D, idx){
  all <- c(rep(1,D[1]), rep(2,D[2]), rep(3,D[3]), 
           rep(4,D[4]), rep(5,D[5]), rep(6,D[6]))
  return(all[idx])
}

#' Component index to covariate index
#' @param D integer vector of length 
#' @param idx integer
#' @return an integer
component_index_to_covariate_index <- function(D, idx){
  covr_idx <- 2 - sum(D[1:2]) + idx
  return(covr_idx)
}
