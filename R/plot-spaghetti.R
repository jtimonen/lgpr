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
                      highlight_col = NULL,
                      highlight_lty = NULL,
                      highlight_cont = NULL,
                      response = "y", 
                      palette = "Set1",
                      colgrad = "gray",
                      psize = 2,
                      title = NULL,
                      lwd = 0.5)
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
        colgrad <- ggplot2::scale_colour_gradient(low = "#c8e2ed", 
                                                  high = "#1e2a30",
                                                  space = "Lab", 
                                                  na.value = "firebrick3", 
                                                  guide = "colourbar",
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


#' Plot computed predictions componentwise
#' @export
#' @description NOTE: currently assumes that the case individuals come first!
#' @param fit An object of class \code{lgpfit}.
#' @param PRED Predictions computed using \code{lgp_predict}.
#' @param componentwise Should the predictions be plotted componentwise?
#' @param color_scheme Name of bayesplot color scheme or a list with fieds 'dark' and 'light'.
#' @param alpha Ribbon fill opacity.
#' @param alpha_line Line opacity
#' @param plot_pred should the (mean) predictions line be plotted?
#' @param plot_ribbon Should an uncertainty ribbon be plotted?
#' @param plot_ribbon_edge Should an uncertainty ribbon edge be plotted?
#' @param theme ggplot theme
#' @param original_y_scale should the predictions be scaled back to the original data y scale
#' @param title optional prefix to plot title
#' @param ylim y axis limits
#' @param plot_obs_onset should the observed disease onset be plotted by a vertical line
#' @param plot_onset_samples should a distribution of sampled onsets be plotted
#' @param alpha2 alpha of t_onset density
#' @param color_scheme_onset color scheme name for onset density plotting
#' @param ypos_dens y-position of the density plot
#' @param test_data Test data frame
#' @param color_test test point color
#' @param pch_test test point marker
#' @param size_test test point size
#' @return a ggplot object
plot_predictions <- function(fit, 
                             PRED             = NULL, 
                             componentwise    = FALSE,
                             color_scheme     = "red",
                             alpha            = 0.5-0.4*as.numeric(componentwise), 
                             plot_ribbon      = TRUE,
                             title            = NULL,
                             theme            = ggplot2::theme_gray(),
                             original_y_scale = TRUE,
                             plot_ribbon_edge = !componentwise,
                             plot_pred        = TRUE,
                             alpha_line       = alpha,
                             ylim             = NULL,
                             plot_obs_onset   = FALSE,
                             alpha2           = 0.5,
                             color_scheme_onset = "gray",
                             plot_onset_samples = FALSE,
                             ypos_dens        = NULL,
                             test_data        = NULL,
                             color_test       = "steelblue4",
                             pch_test         = 8,
                             size_test        = 2)
{
  # Check input correctness
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  model <- fit@model
  if(class(model)!="lgpmodel") stop("Class of 'fit@model' must be 'lgpmodel'!")
  info    <- model@info
  idvar   <- info$varInfo$id_variable
  timevar <- info$varInfo$time_variable
  respvar <- info$varInfo$response_variable
  
  # Get color
  if(class(color_scheme)=="character"){
    color_scheme <- bayesplot::color_scheme_get(color_scheme)
  }
  if(class(color_scheme_onset)=="character"){
    color_scheme_onset <- bayesplot::color_scheme_get(color_scheme_onset)
  }
  linecolor <- color_scheme$dark
  if(!componentwise){
    fill      <- color_scheme$mid
    hlcolor   <- color_scheme$mid_highlight
  }else{
    fill      <- color_scheme$dark
    hlcolor   <- color_scheme$dark_highlight
  }
  # Create plot data frame
  if(is.null(PRED)){
    if(componentwise){
      stop("Plotting predictions componentwise not implemented in this case!")
    }
    DF <- create_predictions_plot_df1(fit, 
                                      scale_f = original_y_scale)
  }else{
    if(info$sample_F){
      stop("Do not give predictions as input for a model where F was sampled!")
    }
    DF <- create_predictions_plot_df2(model, 
                                      PRED, 
                                      scale_f = original_y_scale, 
                                      componentwise = componentwise)
  }
  
  if(componentwise){
    facet_var    <- "component"
    DF$facet_var <- DF[[facet_var]]
    ptitle       <- "Componentwise predictions"
    ylab         <- " "
  }else{
    facet_var    <- idvar
    DF$facet_var <- DF[[facet_var]]
    ptitle       <- "Model predictions"
    ylab         <- respvar
  }
  
  # Create ggplot object
  if(info$sample_F){
    group_var <- "idx"
    y_var <- "pred"
  }else{
    group_var <- idvar
    y_var <- "mu"
  }
  h  <- ggplot2::ggplot(DF, ggplot2::aes_string(x = timevar, 
                                                y = y_var, 
                                                group = group_var))
  
  # Edit plot title
  if(!is.null(title)){
    ptitle <- paste(title, ": ", ptitle, sep="")
  }
  h <- h + ggplot2::ggtitle(label = ptitle)
  h <- h + ggplot2::labs(y = ylab)
  
  # Lines and ribbons
  if(plot_ribbon && !info$sample_F){
    h <- h + ggplot2::geom_ribbon(ggplot2::aes_string(ymin = 'lower', 
                                                      ymax = 'upper'),
                                  fill  = fill, 
                                  alpha = alpha)
  }
  if(!is.null(hlcolor) && plot_ribbon_edge && !info$sample_F){
    h <- h + ggplot2::geom_line(ggplot2::aes_string(x = timevar, y = 'lower'),
                                color = fill, 
                                alpha = alpha)
    h <- h + ggplot2::geom_line(ggplot2::aes_string(x = timevar, y = 'upper'),
                                color = fill, 
                                alpha = alpha)
  }
  if(plot_pred){
    h <- h + ggplot2::geom_line(color = linecolor, alpha = alpha_line)
  }
  
  # Plot also the data
  if(!componentwise){
    dat  <- model@data
    X_id <- as.numeric(model@stan_dat$X[1,])
    df   <- data.frame(idvar   = as.factor(X_id), 
                       timevar  = as.numeric(dat[[timevar]]), 
                       y        = as.numeric(dat[[respvar]]))
    colnames(df)[1] <- idvar
    colnames(df)[2] <- timevar
    colnames(df)[3] <- respvar
    df$facet_var <- df[[facet_var]]
    
    h   <- h + ggplot2::geom_point(data    = df, 
                                   mapping = ggplot2::aes_string(x = timevar,
                                                                 y = respvar,
                                                                 group = idvar))
  }
  # Faceting
  h <- h + ggplot2::facet_wrap(. ~ facet_var)
  
  # Theme
  h <- h + theme
  
  # Y limits
  if(!is.null(ylim)){
    h <- h + ggplot2::coord_cartesian(ylim = ylim)
  }
  
  # Plot observed onsets
  D <- model@stan_dat$D
  if(D[3]==1 && !componentwise && plot_obs_onset){
    df <- model@data
    davar <- fit@model@info$varInfo$disAge_variable
    t_ons <- get_onset_times(id = df[[idvar]], age = df[[timevar]], disAge = df[[davar]])
    vline.data <- data.frame(zzz = t_ons, facet_var = names(t_ons))
    
    h <- h + ggplot2::geom_vline(ggplot2::aes_string(xintercept = "zzz"), 
                                 na.rm    = TRUE,
                                 data     = vline.data, 
                                 color    = "black",
                                 linetype = 3)
    
  }
  
  # Plot onset samples
  UNCRT <- model@stan_dat$UNCRT
  if(UNCRT==1 && !componentwise && D[3]==1 & plot_onset_samples){

    T_smp <- extract_t_onset_samples(fit)
    n_smp <- dim(T_smp)[1]
    t_smp <- as.numeric(T_smp)
    df    <- model@data
    yvar  <- model@info$varInfo$response_variable
    ystd  <- stats::sd(df[[yvar]])
    dens_scale <- ystd
    if(is.null(ypos_dens)){
      ypos_dens <- min(df[[yvar]]) - 1.5*dens_scale
    }
    
    # This part assumes that case individuals come first
    uid   <- unique(df[[idvar]])
    N     <- length(uid)
    Nc    <- model@stan_dat$N_cases
    t_smp <- c(t_smp, rep(NaN, (N-Nc)*n_smp))
    fv    <- as.factor(rep(uid, each = n_smp))
    dens.data <- data.frame(zzz = t_smp, facet_var = fv, id = fv)
    h <- h + ggplot2::geom_density(ggplot2::aes_string(x = "zzz", y = paste(dens_scale,"*..scaled..",sep="")), 
                                 na.rm    = TRUE,
                                 data     = dens.data, 
                                 color    = color_scheme_onset$mid_highlight,
                                 fill     = color_scheme_onset$mid,
                                 alpha    = alpha2,
                                 inherit.aes = F, 
                                 position = ggplot2::position_nudge(x=0, y=ypos_dens),
                                 trim     = TRUE)
  }
  
  # Plot test data
  if(!is.null(test_data)){
    id_test    <- test_data[[idvar]]
    age_test   <- test_data[[timevar]]
    y_test     <- test_data[[respvar]]
    point.data <- data.frame(xxx = age_test,
                             yyy = y_test,
                             facet_var = id_test)
    
    h <- h + ggplot2::geom_point(mapping  = ggplot2::aes_string(x = "xxx", y = "yyy"), 
                                 na.rm    = TRUE,
                                 inherit.aes = F,
                                 data     = point.data, 
                                 color    = color_test,
                                 pch      = pch_test,
                                 size     = size_test)
  }
  
  return(h)
}


