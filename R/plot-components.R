#' Visualize the (average) inferred components evaluated at data points
#' @export
#' @param fit An object of class \code{lgpfit}.
#' @param corrected Should this plot the covariate-effect corrected components?
#' @param title optional prefix to plot title
#' @param sample_idx If given, only one sample will be plotted, else the average
#' components over all samples.
#' @param linealpha line alpha
#' @param linetype line type
#' @param return_list should this return a list of ggplot objects?
#' @param ncol number of plot columns
#' @param nrow number of plot rows
#' @param legend legend argument for ggarrange, use "none" to remove legends
#' @param labels labels argument for ggarrange
#' @param ylim y axis limits
#' @return alist of ggplot objects or one combined plot
plot_components <- function(fit, 
                            corrected = TRUE, 
                            title = NULL, 
                            sample_idx = NULL,
                            linealpha = 1,
                            linetype = 3,
                            return_list = FALSE,
                            ncol = NULL,
                            nrow = NULL,
                            legend = NULL,
                            labels = "auto",
                            ylim = NULL){
  
  GG <- list()
  D <- fit@model@stan_dat$D
  sum_D <- sum(D)
  for(d in 1:sum_D){
    gg <- plot_component(fit, d, corrected, sample_idx, linealpha, linetype)
    if(is.null(ylim)){
      COMP <- get_inferred_components(fit, corrected, sample_idx)
      fff  <- COMP$EFF
      ylim <- range(fff)
    }
    gg <- gg + ggplot2::ylim(ylim)
    
    if(is.null(legend)){
      gg <- gg + ggplot2::theme(
        legend.key = ggplot2::element_blank(),
        legend.justification = c(1,1), 
        legend.position = c(1,1),
        legend.background = ggplot2::element_rect(fill = NA),
        legend.direction = "horizontal"
        )
    }
    GG[[d]] <- gg
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
#' @param corrected Should this plot a covariate-effect corrected component?
#' @param sample_idx If given, only one sample will be plotted, else the average
#' components over all samples.
#' @param linealpha line alpha
#' @param linetype line type
#' @return a ggplot object
plot_component <- function(fit, idx, corrected = TRUE, 
                           sample_idx = NULL, linealpha = 1,
                           linetype = 3)
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
  idvar   <- info$varInfo$id_variable
  timevar <- info$varInfo$time_variable
  COMP    <- get_inferred_components(fit, corrected, sample_idx)
  title   <- COMP$names[idx]
  EFF     <- COMP$EFF
  relev   <- COMP$relev
  TSCL    <- model@scalings$TSCL
  X       <- t(model@stan_dat$X)[1:n,]
  age     <- TSCL$fun_inv(X[,2])
  id      <- X[,1]
  f       <- EFF[,idx]
  ctype   <- component_index_to_type(D, idx)
  yname   <- get_component_plot_titles(corrected, sample_idx)$yname
  cind    <- component_index_to_covariate_index(D, idx)
  groupvar <- as.factor(id)
  leg      <- ggplot2::theme(legend.position = "none")
  if(ctype==2){
    aes <- ggplot2::aes_string(x = timevar, y = 'f')
    df  <- data.frame(age, f)
  }else if(ctype==3){
    colorvar <- as.factor(model@stan_dat$X_notnan)
    leg <- ggplot2::labs(color = 'Group')
    aes <- ggplot2::aes_string(x = timevar, y = 'f', 
                               group = groupvar,
                               color = colorvar)
    df  <- data.frame(age, f, groupvar, colorvar)
  }else if(ctype == 5 || ctype == 6){
    cname <- get_inferred_components(fit, corrected = TRUE,
                                     sample_idx)$names[idx]
    leg   <- ggplot2::labs(color = cname)
    groupvar <- as.factor(X[,cind])
    colorvar <- as.factor(X[,cind])
    aes <- ggplot2::aes_string(x = timevar, y = 'f', 
                               group = groupvar,
                               color = colorvar)
    df  <- data.frame(age, f, groupvar, colorvar)
  }else{
    aes <- ggplot2::aes_string(x = timevar, y = 'f', group = groupvar)
    df  <- data.frame(age, f, groupvar)
  }
  
  
  # Create ggplot object
  colnames(df)[1] <- timevar
  h   <- ggplot2::ggplot(df, aes)

  # Edit ggplot object
  ptitle <- title
  if(corrected){
    ptitle <- paste('Effect of', ptitle)
  }
  if(!is.null(sample_idx)){
    ptitle <- paste(ptitle, ', MCMC sample ', sample_idx, sep="")
  }
  
  h <- h + leg + 
    ggplot2::geom_line(linetype = linetype, alpha = linealpha) + 
    ggplot2::geom_point() + 
    ggplot2::labs(y = yname) +
    ggplot2::ggtitle(label = ptitle)
  
  # Return ggplot object
  return(h)
  
}


#' Get inferred components
#' @export
#' @param fit An object of class \code{lgpfit}.
#' @param corrected Should this obtain the covariate-effect corrected components?
#' @param sample_idx If given, only one sample will be obtained else the average
#' components over all samples.
#' @return a list
get_inferred_components <- function(fit, corrected = TRUE, 
                                    sample_idx = NULL){
  
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  model <- fit@model
  info  <- model@info
  
  # Get components and relevances
  if(is.null(sample_idx)){
    if(corrected){
      names  <- info$covariate_names
      EFF    <- fit@components_corrected
      relev  <- fit@covariate_relevances$average
    }else{
      names  <- info$component_names
      EFF    <- fit@components
      relev  <- fit@component_relevances$average
    }
  }else{
    if(corrected){
      names  <- info$covariate_names
      EFF    <- extract_components_onesample(fit, sample_idx)$corrected
      relev  <- fit@covariate_relevances$samples[sample_idx,]
    }else{
      names  <- info$component_names
      EFF    <- extract_components_onesample(fit, sample_idx)$components
      relev  <- fit@component_relevances$samples[sample_idx,]
    }
  }

  dp2 <- dim(EFF)[2]
  EFF <- EFF[,1:(dp2-2)]
  colnames(EFF) <- names
  names <- colnames(EFF)
  
  return(list(EFF=EFF,names=names,relev=relev))
}


#' Get inferred components
#' @param corrected Should this obtain the covariate-effect corrected components?
#' @param sample_idx Integer or NULL
#' @return a list
get_component_plot_titles <- function(corrected = TRUE, sample_idx = NULL){
  
  if(is.null(sample_idx)){
    if(corrected){
      ptitle <- "Averaged covariate effects [relevance]"
      yname  <- "(corrected) component"
    }else{
      ptitle <- "Average inferred components [relevance]"
      yname  <- "component"
    }
  }else{
    if(corrected){
      ptitle <- paste("Covariate effects for sample ", sample_idx,
                      " [relevance]", sep="")
      yname  <- "(corrected) component"
    }else{
      ptitle <- paste("Inferred components for sample ", sample_idx,
                      " [relevance]", sep="")
      yname  <- "component"
    }
  }
  return(list(ptitle=ptitle, yname=yname))
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