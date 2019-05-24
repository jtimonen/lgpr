#' Plot simulated data for each individual separately
#'
#' @export
#' @description This is a wrapper for \code{plot_simdata_by_individual} and
#' \code{plot_simdata_by_component}
#' @param simData a list returned by \code{\link{simulate_data}}
#' @param componentwise should each component be plotted separately?
#' @param color_scheme name of bayesplot color scheme
#' @param nrow an argument for \code{ggplot2::facet_wrap}
#' @param ncol an argument for \code{ggplot2::facet_wrap}
#' @param plot_point should \code{geom_point} be used also for the underlying signal?
#' @param vlinecolors colors to plot the real and observed disease onsets if they exist
#' @param vlinetypes line types to use for the real and observed disease onsets if they exist
#' @param i_test test point indices
#' @param color_test test point color
#' @return a ggplot object
plot_simdata <- function(simData, 
                         componentwise = FALSE,
                         color_scheme  = "red", 
                         nrow          = NULL, 
                         ncol          = NULL,
                         plot_point    = componentwise,
                         vlinecolors   = rep("steelblue4",2),
                         vlinetypes    = c(1,2),
                         i_test        = NULL,
                         color_test    = "steelblue2")
{
  
  # Color
  color_scheme <- bayesplot::color_scheme_get(color_scheme)
  linecolor    <- color_scheme$mid
  pointcolor   <- color_scheme$mid_highlight
  if(componentwise){
    if(!is.null(i_test)){
      stop("cannot specify test points when plotting componentwise")
    }
    h <- plot_simdata_by_component(simData, linecolor, nrow, ncol, 
                                   plot_point, pointcolor)
  }else{
    h <- plot_simdata_by_individual(simData, linecolor, nrow, ncol, 
                                    plot_point, vlinecolors, vlinetypes,
                                    i_test, color_test)
  }
  return(h)
}


#' Plot a simulated longitudinal data set for each individual separately
#'
#' @param simData a list returned by \code{\link{simulate_data}}
#' @param linecolor line color
#' @param nrow an argument for \code{ggplot2::facet_wrap}
#' @param ncol an argument for \code{ggplot2::facet_wrap}
#' @param plot_point should \code{geom_point} be used also for the underlying signal?
#' @param vlinecolors colors to plot the real and observed disease onsets if they exist
#' @param vlinetypes line types to use for the real and observed disease onsets if they exist
#' @param i_test test point indices
#' @param color_test test point color
#' @return a ggplot object
plot_simdata_by_individual <- function(simData, 
                                       linecolor   = "firebrick3", 
                                       nrow        = NULL, 
                                       ncol        = NULL,
                                       plot_point  = FALSE,
                                       vlinecolors = rep("steelblue4",2),
                                       vlinetypes  = c(1,2),
                                       i_test      = NULL,
                                       color_test  = "steelblue2"){
  
  dat  <- simData$data
  comp <- simData$components
  ons  <- simData$onsets
  D3   <- ("diseaseAge" %in% colnames(comp))
  g    <- comp$g
  y    <- dat$y
  yval <- c(g,y)
  leg  <- rep(c("g", "y"), each = length(g))
  id   <- rep(dat$id, 2)
  age  <- rep(dat$age, 2)
  N    <- length(unique(dat$id))
  n    <- length(dat$id)
  
  if(plot_point){
    stop("points must be plotted in this mode! set plot_point = TRUE")
  }
  
  DF      <- data.frame(id, age, yval, leg)
  DF$age  <- as.numeric(age)
  DF$yval <- as.numeric(yval)
  is_test <- rep("y_train", n)
  if(!is.null(i_test)){
    is_test[i_test] <- "y_test"
  }
  it <- rep(is_test, 2)
  DF$is_test <- as.factor(it)
  
  # Create ggplot object
  h <- ggplot2::ggplot(data = DF, ggplot2::aes_string(x='age', y='yval', group='leg'))
  subt <- paste(n, " data points, ", N, " individuals", sep="")
  
  # Faceting
  h <- h + ggplot2::facet_wrap(.~ id, nrow = nrow, ncol = ncol)
  
  # Plot real and observed disease onsets
  ons1 <- simData$onsets
  ons2 <- simData$onsets_observed
  no1  <- sum(!is.nan(ons1))
  if(no1 > 0){
    
    subt <- paste(subt,". Vertical lines are the real and observed disease onset.", sep ="")
    vline1.data <- data.frame(z = ons1, id = c(1:N))
    vline2.data <- data.frame(z = ons2, id = c(1:N))
    
    # Plot real onsets
    h <- h + ggplot2::geom_vline(ggplot2::aes_string(xintercept = "z"), 
                                 na.rm    = TRUE,
                                 data     = vline1.data, 
                                 color    = vlinecolors[1],
                                 linetype = vlinetypes[1])
    
    # Plot observed onsets
    h <- h + ggplot2::geom_vline(ggplot2::aes_string(xintercept = "z"), 
                                 na.rm    = TRUE,
                                 data     = vline2.data, 
                                 color    = vlinecolors[2],
                                 linetype = vlinetypes[2])
    
  }
  
  # Plot signal line and data points
  h <- h + ggplot2::geom_line(ggplot2::aes_string(linetype = 'leg'), color = linecolor) +
    ggplot2::scale_shape_manual(values=c(NA, 16)) +
    ggplot2::scale_linetype_manual(values=c(1, 0))
  
  h <- h + ggplot2::geom_point(ggplot2::aes_string(shape = 'leg', color = "is_test"), na.rm = TRUE)
  h <- h + ggplot2::labs(x = "Age", y = "y")
  
  # Theme and titles
  h <- h + ggplot2::ggtitle("Simulated data", subtitle = subt)
  h <- h + ggplot2::theme_bw()
  h <- h + ggplot2::theme(legend.title=ggplot2::element_blank())
  
  # Point color and type
  if(!is.null(i_test)){
    h <- h + ggplot2::scale_color_manual(values=c(color_test, "black"))
  }
  return(h)
}


#' Plot each component of a simulated longitudinal data set separately
#'
#' @param simData a list returned by \code{\link{simulate_data}}
#' @param linecolor line color
#' @param pointcolor point color
#' @param nrow an argument for \code{ggplot2::facet_wrap}
#' @param ncol an argument for \code{ggplot2::facet_wrap}
#' @param plot_point should also \code{geom_point} be used?
#' @return a ggplot object
plot_simdata_by_component <- function(simData, 
                                      linecolor = "black", 
                                      nrow = NULL, 
                                      ncol = NULL,
                                      plot_point = TRUE,
                                      pointcolor = "black"){
  DF <- create_simdata_plot_df(simData)
  h <- ggplot2::ggplot(DF, ggplot2::aes_string(x='age', y='value', group='id'))
  h <- h + ggplot2::geom_line(color = linecolor)
  if(plot_point){
    h <- h + ggplot2::geom_point(color = pointcolor)
  }
  h <- h + ggplot2::facet_wrap(.~ component, nrow = nrow, ncol = ncol)
  
  dat  <- simData$data
  N    <- length(unique(dat$id))
  n    <- length(dat$id)
  subt <- paste(n, " data points, ", N, " individuals", sep="")
  h <- h + ggplot2::ggtitle("Simulated data", subtitle = subt)
  h <- h + ggplot2::labs(x = "Age", y = "value")
  h <- h + ggplot2::theme_bw()
  h <- h + ggplot2::theme(legend.title=ggplot2::element_blank())
  
  return(h)
}
