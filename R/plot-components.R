
#' Visualize inferred components
#' @export
#' @inheritParams plot_components_posterior_sub1
#' @inheritParams plot_components_posterior_sub2
#' @return an object returned by ggpubr::ggarrange
plot_components_posterior <- function(fit, subsamples = NULL,
                                      time_is_xvar = TRUE, 
                                      PRED = NULL, 
                                      marker = NULL, 
                                      sample_idx = 1,
                                      n_sd = 2,
                                      ...)
{
  if(is.null(PRED)){
    if(is.null(marker) && is.null(subsamples)){
      marker <- 16
    }
    h <- plot_components_posterior_sub1(fit, subsamples, 
                                        time_is_xvar, marker, ...)
  }else{
    h <- plot_components_posterior_sub2(fit, PRED, sample_idx,
                                        time_is_xvar, n_sd, ...)
  }
  return(h)
}



#' Helper for \code{\link{plot_components_posterior}}
#' @param fit An object of class \code{lgpfit}.
#' @param time_is_xvar is the time variable the x-axis variable
#' in all subplots?
#' @param subsamples How many samples to plot. If this is NULL,
#' average over all samples is plotted. If this is "all", all
#' samples are plotted.
#' @param marker point type
#' @param ... additional arguments for \code{\link{plot_components}}
#' @return an object returned by ggpubr::ggarrange
plot_components_posterior_sub1 <- function(fit, subsamples,
                                           time_is_xvar, marker,
                                           ...)
{
  # Check input correctness
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  model <- fit@model
  
  # Create plot
  df  <- as.data.frame(fit@stan_fit)
  FFF <- get_function_components_from_df_all(df, model)
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
  
  h <- plot_components(FFF, NULL, model, 
                       time_is_xvar, marker = marker, ...)
  
  # Return ggplot object
  return(h)
  
}

#' Helper for \code{\link{plot_components_posterior}}
#' @param fit An object of class \code{lgpfit}.
#' @param n_sd number of standard deviations (ribbon width)
#' @param PRED object returned by \code{\link{lgp_predict}}
#' @param time_is_xvar is the time variable the x-axis variable
#' in all subplots?
#' @param sample_idx Which sample to plot. 
#' @param ... additional arguments for \code{\link{plot_components}}
#' @return an object returned by ggpubr::ggarrange
plot_components_posterior_sub2 <- function(fit, PRED, sample_idx,
                                           time_is_xvar,
                                           n_sd, ...)
{
  # Check input correctness
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  model <- fit@model
  # Create plot
  AAA <- PRED_to_arrays(PRED)
  MMM <- AAA$MMM
  SSS <- n_sd * AAA$SSS
  X_test <- PRED$X_test_scaled
  
  h <- plot_components(MMM, SSS, model, time_is_xvar, X_test, ...)
  
  # Return ggplot object
  return(h)
  
}

#' Visualize the components of a simulated data set
#' @export
#' @param simData simulated data object (list)
#' @param time_is_xvar is the time variable the x-axis variable
#' in all subplots?
#' @param marker point marker
#' @param ... additional arguments for \code{\link{plot_components}}
#' @return an object returned by ggpubr::ggarrange
plot_components_simdata <- function(simData, 
                                    time_is_xvar = TRUE, 
                                    marker = 16, ...)
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
  
  h <- plot_components(FFF, NULL, model, time_is_xvar, 
                       marker = marker, ...)
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
                            X_test        = NULL,
                            sum_highlight = NULL,
                            linealpha   = 1,
                            linetype    = 1, 
                            fill_alpha  = 0.3,
                            marker      = NULL,
                            ncol        = NULL,
                            nrow        = NULL,
                            legend      = NULL,
                            labels      = NULL,
                            ylim        = NULL,
                            font_size   = 9,
                            theme       = ggplot2::theme_linedraw(),
                            legend_dir  = "horizontal",
                            xlabel      = NULL,
                            ylabel      = " ",
                            viridis_option = "viridis"){
  
  GG <- list()
  sum_D <- sum(model@stan_dat$D) + 1
  for(d in 1:sum_D){
    gg <- plot_component(MMM, SSS, model, d, time_is_xvar,
                         linealpha, linetype, fill_alpha,
                         X_test, marker, sum_highlight, viridis_option)
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
#' @param time_is_xvar is the time variable the x-axis variable 
#' @param linealpha line alpha
#' @param linetype line type
#' @param fill_alpha fill alpha for geom_ribbons
#' @param X_test optional matrix of test points
#' @param marker point type
#' @param sum_highlight name of a categorical covariate to be highlighted
#' @param viridis_option the option argument of \code{ggplot2::scale_colour_viridis_c}
#' by colour in the sum plot
#' @return a ggplot object
plot_component <- function(MMM, SSS, model, idx, time_is_xvar,
                           linealpha, linetype, fill_alpha,
                           X_test, marker, sum_highlight, viridis_option)
{
  sdat  <- model@stan_dat
  if(is.null(X_test)){
    X     <- t(sdat$X)
    Xnn   <- sdat$X_notnan
  }else{
    X     <- X_test
    Xnn   <- as.numeric(!is.nan(X_test[,3]))
  }
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
  
  plot_eb  <- !is.null(SSS) && !(ctype %in% c(1,7))
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
    xvar_name <- cvn[cind]
    if(ctype==4){
      warning('covariate', xvar_name, 'is not on original scale!')
    }
  }
  
  
  # Create df and aes
  if(ctype==1){
    leg      <- ggplot2::theme(legend.position = "none")
    grpvar   <- as.factor(paste(id, sample))
    df       <- data.frame(xvar, f, grpvar)
    aes      <- ggplot2::aes_string(x = 'xvar', y = 'f', 
                                    group = 'grpvar')
    aes_eb   <- ggplot2::aes_string(ymin='lb', ymax='ub')
  }else if(ctype==2){
    leg      <- ggplot2::theme(legend.position = "none")
    grpvar   <- as.factor(sample)
    df       <- data.frame(xvar, f, grpvar)
    aes      <- ggplot2::aes_string(x = 'xvar', y = 'f', group = 'grpvar')
    aes      <- ggplot2::aes_string(x = 'xvar', y = 'f', 
                                    group = 'grpvar')
    aes_eb   <- ggplot2::aes_string(ymin='lb', ymax='ub')
  }else if(ctype==3){
    colorvar <- as.factor(rep(Xnn, S))
    leg      <- ggplot2::labs(color = 'Group', fill = 'Group')
    grpvar   <- as.factor(paste(id, sample))
    df       <- data.frame(xvar, f, grpvar, colorvar)
    aes      <- ggplot2::aes_string(x = 'xvar', y = 'f', 
                                    group = 'grpvar', color = 'colorvar')
    aes_eb   <- ggplot2::aes_string(ymin='lb', ymax='ub', fill = 'colorvar')
  }else if(ctype==4){
    leg      <- ggplot2::theme(legend.position = "none")
    grpvar   <- as.factor(paste(id, sample))
    colorvar <- rep(X[,cind], S)
    df       <- data.frame(xvar, f, grpvar, colorvar)
    aes      <- ggplot2::aes_string(x = 'xvar', y = 'f', 
                                    group = 'grpvar',
                                    color = 'colorvar')
    aes_eb   <- ggplot2::aes_string(ymin='lb', ymax='ub')
  }else if(ctype == 5 || ctype == 6){
    leg      <- ggplot2::labs(color = cvn[cind], fill = cvn[cind])
    colorvar <- rep(X[,cind], S)
    grpvar   <- as.factor(paste(colorvar, sample))
    colorvar <- as.factor(colorvar)
    df       <- data.frame(xvar, f, grpvar, colorvar)
    aes      <- ggplot2::aes_string(x = 'xvar', y = 'f', 
                                    group = 'grpvar', color = 'colorvar')
    aes_eb   <- ggplot2::aes_string(ymin='lb', ymax='ub', fill = 'colorvar')
  }else if(ctype==7){
    if(!is.null(sum_highlight)){
      if(sum_highlight=='group'){
        colorvar <- as.factor(rep(Xnn, S))
        leg      <- ggplot2::labs(color = 'Group', fill = 'Group')
        grpvar   <- as.factor(paste(id, sample))
        df       <- data.frame(xvar, f, grpvar, colorvar)
        aes      <- ggplot2::aes_string(x = 'xvar', y = 'f', 
                                        group = 'grpvar',
                                        color = 'colorvar')
        aes_eb   <- ggplot2::aes_string(ymin='lb', ymax='ub',
                                        fill ='colorvar')
      }else {
        cind <- which(cvn==sum_highlight)
        if(length(cind)==0){stop("invalid sum_highlight (", sum_highlight,")")}
        leg      <- ggplot2::labs(color = cvn[cind], fill = cvn[cind])
        colorvar <- rep(X[,cind], S)
        grpvar   <- as.factor(paste(id,sample))
        colorvar <- as.factor(colorvar)
        df       <- data.frame(xvar, f, grpvar, colorvar)
        aes      <- ggplot2::aes_string(x = 'xvar', y = 'f', 
                                        group = 'grpvar', 
                                        color = 'colorvar')
        aes_eb   <- ggplot2::aes_string(ymin='lb', ymax='ub', 
                                        fill = 'colorvar')
      }
    }else{
      leg      <- ggplot2::theme(legend.position = "none")
      grpvar   <- as.factor(paste(id, sample))
      df       <- data.frame(xvar, f, grpvar)
      aes      <- ggplot2::aes_string(x = 'xvar', y = 'f', group = 'grpvar')
      aes_eb   <- ggplot2::aes_string(ymin='lb', ymax='ub')
    }
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
  h   <- ggplot2::ggplot(df, aes)
  if(plot_eb){
    h <- h + ggplot2::geom_ribbon(aes_eb, alpha = fill_alpha, 
                                  lty = 0)
  }
  if(ctype==4){
    h <- h + ggplot2::scale_colour_viridis_c(option = viridis_option) 
  }
  h   <- h + leg
  h   <- h + ggplot2::geom_line(linetype = linetype, alpha = linealpha)
  h   <- h + ggplot2::ggtitle(title)
  h   <- h + ggplot2::xlab(xvar_name)
  
  # Add marker
  if(!is.null(marker)){
    h <- h + ggplot2::geom_point(pch = marker)
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
