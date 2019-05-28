#' A spaghetti plot of longitudinal data.
#' @export
#' @param data A data frame. Must contain columns \code{id} and \code{age}.
#' @param highlight_col Name of a categorical covariate to be highlighted with color.
#' @param highlight_lty Name of a categorical covariate to be highlighted with line style.
#' @param highlight_cont Name of a continuous covariate to be highlighted.
#' @param response Name of the response variable (\code{default ="y"}).
#' @param colgrad A colour gradient type, must be either \code{"gray"}
#' (default), or \code{"diverging"}.
#' @param palette Name of a ColorBrewer palette.
#' @param psize Point 
#' @param title optional prefix to plot title
#' @param lwd line width
#' @return a ggplot object
plot_data <- function(data, 
                      highlight_col  = NULL,
                      highlight_lty  = NULL,
                      highlight_cont = NULL,
                      response       = "y", 
                      palette        = "Set1",
                      colgrad        = "gray",
                      psize          = 2,
                      title          = NULL,
                      lwd            = 0.5)
{
  data$id  <- as.factor(data$id)
  iresp    <- which(colnames(data)==response)
  colnames(data)[iresp] <- "Response"
  
  # Are there diseased and healthy individuals?
  if("diseaseAge" %in% colnames(data)){
    Group <- as.numeric(!is.nan(data$diseaseAge))
    data <- cbind(data, Group)
  }
  
  # Extend data frame
  data <- create_data_plot_df(data, highlight_col, highlight_lty, highlight_cont)
  
  # Create plot
  p <- ggplot2::ggplot(data = data, ggplot2::aes_string(x = 'age', y = 'Response', group = 'id'))
  
  # Highlight categorical
  
  # i)
  if(!is.null(highlight_col) && is.null(highlight_lty)){
    p <- p + ggplot2::labs(colour = highlight_col)
    p <- p + ggplot2::scale_color_brewer(palette = palette)
    p <- p + ggplot2::geom_line(ggplot2::aes_string(color = 'Z1'), size = lwd)
    if(psize > 0){
      p <- p + ggplot2::geom_point(ggplot2::aes_string(color = 'Z1'), size = psize)
    }
  }
  
  # ii)
  if(is.null(highlight_col) && !is.null(highlight_lty)){
    p <- p + ggplot2::labs(linetype = highlight_lty)
    p <- p + ggplot2::geom_line(ggplot2::aes_string(linetype = 'Z2'), size = lwd)
    if(psize > 0){
      p <- p + ggplot2::geom_point(size = psize)
    }
  }
  
  # iii)
  if(!is.null(highlight_col) && !is.null(highlight_lty)){
    p <- p + ggplot2::labs(colour = highlight_col, linetype = highlight_lty)
    p <- p + ggplot2::scale_color_brewer(palette = palette)
    p <- p + ggplot2::geom_line(ggplot2::aes_string(color = 'Z1', linetype = 'Z2'), size = lwd)
    if(psize > 0){
      p <- p + ggplot2::geom_point(ggplot2::aes_string(color = 'Z1'), size = psize)
    }
  }
  
  # Highligh continuous
  if(!is.null(highlight_cont)){
    if(is.null(highlight_col) && is.null(highlight_lty)){
      p <- p + ggplot2::geom_line(ggplot2::aes_string(color = 'Value'), size = lwd)
      if(colgrad!="diverging"){
        colgrad <- ggplot2::scale_colour_gradient(low        = "#c8e2ed", 
                                                  high       = "#1e2a30",
                                                  space      = "Lab", 
                                                  na.value   = "firebrick3", 
                                                  guide      = "colourbar",
                                                  aesthetics = "colour")
      }else{
        colgrad <- ggplot2::scale_colour_gradient2()
      }
      p <- p + colgrad
      p <- p + ggplot2::labs(colour = highlight_cont)
    }else{
      stop("Cannot highlight both a categorical and continuous covariate!")
    }
  }
  
  # Normal line
  if(is.null(highlight_col) && is.null(highlight_lty) && is.null(highlight_cont)){
    p <- p + ggplot2::geom_line(size = lwd) 
    if(psize > 0){
      p <- p + ggplot2::geom_point(size = psize) 
    }
  }
  
  # x label
  p <- p + ggplot2::labs(x = "Age")
  
  # Edit plot title
  ptitle <- "Data"
  subtitle <- paste(length(unique(data$id))," individuals, ",
                    length(data$id), " data points", sep="")
  if(!is.null(title)){
    ptitle <- paste(title, ": ", ptitle, sep="")
  }
  p <- p + ggplot2::ggtitle(ptitle, subtitle = subtitle)
  return(p)
}


#' Visualize the (average) inferred components evaluated at data points
#' @export
#' @param fit An object of class \code{lgpfit}.
#' @param corrected Should this plot the covariate-effect corrected components?
#' @param title optional prefix to plot title
#' @param sample_idx If given, only one sample will be plotted, else the average
#' components over all samples.
#' @return a ggplot object
plot_components <- function(fit, corrected = TRUE, title = NULL, sample_idx = NULL){
  
  
  # Check input correctness
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  model <- fit@model
  if(class(model)!="lgpmodel") stop("Class of 'fit@model' must be 'lgpmodel'!")
  info <- model@info
  n    <- model@stan_dat$n
  D    <- model@stan_dat$D
  if(sum(D) <=1){
    stop("A model must have more than 1 component to use this!")
  }
  idvar   <- info$varInfo$id_variable
  timevar <- info$varInfo$time_variable
  
  # Get components and relevances
  if(is.null(sample_idx)){
    if(corrected){
      names  <- info$covariate_names
      ptitle <- "Averaged covariate effects [relevance]"
      yname  <- "(corrected) component"
      EFF    <- fit@components_corrected
      relev  <- fit@covariate_relevances$average
    }else{
      names  <- info$component_names
      ptitle <- "Average inferred components [relevance]"
      yname  <- "component"
      EFF    <- fit@components
      relev  <- fit@component_relevances$average
    }
  }else{
    if(corrected){
      names  <- info$covariate_names
      ptitle <- paste("Covariate effects for sample ", sample_idx, " [relevance]", sep="")
      yname  <- "(corrected) component"
      EFF    <- extract_components_onesample(fit, sample_idx)$corrected
      relev  <- fit@covariate_relevances$samples[sample_idx,]
    }else{
      names  <- info$component_names
      ptitle <- paste("Inferred components for sample ", sample_idx, " [relevance]", sep="")
      yname  <- "component"
      EFF    <- extract_components_onesample(fit, sample_idx)$components
      relev  <- fit@component_relevances$samples[sample_idx,]
    }
  }
  
  
  # Remove f, g, and residual
  dp2 <- dim(EFF)[2]
  EFF <- EFF[,1:(dp2-2)]
  colnames(EFF) <- names
  names <- colnames(EFF)
  
  # Edit plot title
  if(!is.null(title)){
    ptitle <- paste(title, ": ", ptitle, sep="")
  }
  
  # Create plot
  TSCL  <- model@scalings$TSCL
  X     <- t(model@stan_dat$X)
  X     <- X[1:n,]
  n     <- dim(EFF)[1]
  d     <- dim(EFF)[2]
  X_age <- TSCL$fun_inv(X[,2])
  t     <- rep(X_age, d)
  id    <- rep(X[,1], d)
  y     <- rep(0, length(t))
  H     <- rep('foo', length(t))
  for(j in 1:d){
    inds <- ((j-1)*n + 1):(j*n)
    y[inds] <- EFF[,j]
    H[inds] <- names[j]
  }
  component <- as.factor(H)
  relev_ord <- relev[levels(component)]
  headers   <- paste(levels(component), "   [", round(relev_ord, 3),"]",sep="")
  levels(component) <- headers
  age <- t
  f   <- y
  df  <- data.frame(id, age, f, component)
  colnames(df)[1] <- idvar
  colnames(df)[2] <- timevar
  h   <- ggplot2::ggplot(df, ggplot2::aes_string(x = timevar, 
                                                 y = 'f', 
                                                 group = idvar)) + 
    ggplot2::geom_line() + 
    ggplot2::geom_point() + 
    ggplot2::facet_wrap(. ~ component) +
    ggplot2::labs(y = yname) +
    ggplot2::ggtitle(label = ptitle)
  return(h)
}

