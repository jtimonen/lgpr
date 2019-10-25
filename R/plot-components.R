#' Visualize the (average) inferred components evaluated at data points
#' @export
#' @param fit An object of class \code{lgpfit}.
#' @param title optional prefix to plot title
#' @param linealpha line alpha
#' @param linetype line type
#' @param return_list should this return a list of ggplot objects?
#' @param ncol number of plot columns
#' @param nrow number of plot rows
#' @param legend legend argument for ggarrange, use "none" to remove legends
#' @param labels labels argument for ggarrange
#' @param ylim y axis limits
#' @param font_size font size for plots
#' @param theme ggplot theme
#' @param legend.direction direction of legend
#' @param xlabel x-axis label
#' @param ylabel y-axis label 
#' @param subsamples how many samples to plot
#' @return alist of ggplot objects or one combined plot
plot_components <- function(fit, 
                            title = NULL,
                            return_list = FALSE,
                            ncol = NULL,
                            nrow = NULL,
                            legend = NULL,
                            labels = NULL,
                            ylim = NULL,
                            font_size = 9,
                            theme = ggplot2::theme_linedraw(),
                            legend.direction = "horizontal",
                            xlabel = NULL,
                            ylabel = NULL,
                            linealpha = 0.2,
                            linetype = 1,
                            subsamples = 100){
  
  GG <- list()
  D <- fit@model@stan_dat$D
  sum_D <- sum(D)
  for(d in 1:sum_D){
    gg <- plot_component(fit, d, linealpha, linetype, subsamples)
    if(is.null(ylim)){
      ylim <- range(fit@components)
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
        legend.direction = legend.direction
      )
    }
    GG[[d]] <- gg + ggplot2::theme(text=ggplot2::element_text(size=font_size))
  }
  if(return_list){
    return(GG)
  }else{
    p <- ggpubr::ggarrange(plotlist = GG, ncol = ncol, nrow = nrow,
                           legend = legend,
                           labels = labels)
    return(p)
  }
}

#' Visualize one (average) inferred component evaluated at data points
#' @export
#' @param fit An object of class \code{lgpfit}.
#' @param idx Index of component to be plotted.
#' @param linealpha line alpha
#' @param linetype line type
#' @param subsamples how many samples to plot
#' @return a ggplot object
plot_component <- function(fit, idx,
                           linealpha = 0.2,
                           linetype = 1,
                           subsamples = 100)
{
  # Check input correctness
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  model <- fit@model
  info  <- model@info
  n     <- model@stan_dat$n
  D     <- model@stan_dat$D
  sum_D <- sum(D)
  if(idx < 1 || idx > sum_D){
    stop("idx must be and integer between 1 and ", sum_D, "!")
  }
  
  # Create plot
  EFF     <- get_function_component_samples(fit)
  ddd     <- dim(EFF)[3]
  EFF     <- EFF[,,1:(ddd-2)]
  names   <- colnames(EFF[1,,])
  title   <- names[idx]
  EFF     <- EFF[,,idx] # EFF is now n_samples x n_data
  
  n_smp   <- dim(EFF)[1]
  i_sub   <- sample.int(n_smp, subsamples)
  EFF     <- EFF[i_sub,] 
  
  TSCL    <- model@scalings$TSCL
  X       <- t(model@stan_dat$X)[1:n,]
  age     <- TSCL$fun_inv(X[,2])
  id      <- X[,1]
  ptitle  <- title
  
  cvn <- model@info$covariate_names
  h   <- plot_component_create_gg(EFF, id, age, idx, D, cvn, X)
  h   <- h + ggplot2::geom_line(linetype = linetype, alpha = linealpha) + 
    ggplot2::ggtitle(label = ptitle)
  
  # Return ggplot object
  return(h)
  
}


#' Subroutine for plot_component
#' 
#' @param EFF a matrix of size n_samples x n_data
#' @param id id covariate
#' @param age age covariate
#' @param idx covariate index
#' @param D an integer vector of length 6
#' @param covariate_names covariate names
#' @param X covariate matrix
#' @return a ggplot object
plot_component_create_gg <- function(EFF, id, age, idx, D, covariate_names, X){

  # Create
  ctype    <- component_index_to_type(D, idx)
  yname    <- " "
  cind     <- component_index_to_covariate_index(D, idx)
  leg      <- ggplot2::theme(legend.position = "none")
  if(ctype==2){
    df     <- plot_component_create_df(EFF, NULL, NULL, age)
    aes    <- ggplot2::aes_string(x = 'timevar', y = 'f', group = 'sample')
  }else if(ctype==3){
    colorvar <- as.factor(model@stan_dat$X_notnan)
    leg <- ggplot2::labs(color = 'Group')
    aes <- ggplot2::aes_string(x = 'timevar', y = 'f', 
                               group = 'groupvar',
                               color = 'colorvar')
    df  <- plot_component_create_df(EFF, id, colorvar, age)
  }else if(ctype == 5 || ctype == 6){
    cname    <- covariate_names[cind]
    leg      <- ggplot2::labs(color = cname)
    groupvar <- as.factor(X[,cind])
    colorvar <- as.factor(X[,cind])
    df       <- plot_component_create_df(EFF, groupvar, colorvar, age)
    df$grp   <- paste(df$groupvar, df$sample)
    aes      <- ggplot2::aes_string(x = 'timevar', y = 'f', 
                                    group = 'grp',
                                    color = 'colorvar')
  }else{
    df     <- plot_component_create_df(EFF, id, NULL, age)
    df$grp <- paste(df$groupvar, df$sample)
    aes    <- ggplot2::aes_string(x = 'timevar', y = 'f', group = 'grp')
  }

  # Create ggplot object
  colnames(df)[1] <- 'timevar'
  h <- ggplot2::ggplot(df, aes) + leg + ggplot2::labs(y = yname)
  return(h)
  
}


#' Create dataframe for componentwise plotting
#' @export
#' @param FFF a matrix of size n_samples x n_data
#' @param groupvar grouping variable
#' @param colorvar coloring variable
#' @param timevar time variable
#' @return a data frame
plot_component_create_df <- function(FFF, groupvar, colorvar, timevar){
  n_smp    <- dim(FFF)[1]
  n        <- dim(FFF)[2]
  f        <- as.numeric(t(FFF))
  sample   <- as.factor(rep(1:n_smp, each = n))
  timevar  <- rep(timevar, n_smp)
  if(is.null(colorvar) && is.null(groupvar)){
    df <- data.frame(timevar, f, sample)
  }else if(is.null(colorvar)){
    groupvar <- as.factor(rep(groupvar, n_smp))
    df <- data.frame(timevar, groupvar, f, sample)
  }else if(is.null(groupvar)){
    colorvar <- as.factor(rep(colorvar, n_smp))
    df <- data.frame(timevar, colorvar, f, sample)
  }else{
    groupvar <- as.factor(rep(groupvar, n_smp))
    colorvar <- as.factor(rep(colorvar, n_smp))
    df <- data.frame(timevar, f, groupvar, colorvar, sample)
  }
  return(df)
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