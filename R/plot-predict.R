#' Plot posterior of f
#' 
#' @export
#' @description This is a wrapper for \code{\link{plot_posterior_predictions}}.
#' @param fit An object of class \code{lgpfit}.
#' @param PRED Predictions computed using \code{lgp_predict}.
#' @param plot_uncertainty Should an uncertainty ribbon be plotted?
#' @param n_sds number of standard deviations for the uncertainty band width
#' @param data_marker pch for data points
#' 
#' @return a ggplot object
plot_posterior_f <- function(fit, 
                             PRED             = NULL, 
                             plot_uncertainty = TRUE,
                             data_marker      = 16,
                             n_sds            = 2){
  
  plot_obs_onset     <- FALSE
  plot_t_effect_samples <- FALSE
  if(fit@model@stan_dat$UNCRT == 1){
    plot_obs_onset     <- TRUE
    plot_t_effect_samples <- TRUE
  }
  if(fit@model@info$sample_F){
    alpha_line <- 0.1
  }else{
    alpha_line <- 1
  }
  h <- plot_posterior_predictions(fit, 
                                  mode               = "posterior", 
                                  PRED               = PRED, 
                                  plot_uncertainty   = plot_uncertainty,
                                  n_sds              = n_sds,
                                  plot_obs_onset     = plot_obs_onset,
                                  plot_t_effect_samples = plot_t_effect_samples,
                                  alpha_line         = alpha_line,
                                  data_marker        = data_marker)
  return(h)
}


#' Plot posterior predictive distribution
#' 
#' @export
#' @description This is a wrapper for \code{\link{plot_posterior_predictions}}.
#' @param fit An object of class \code{lgpfit}.
#' @param PRED Predictions computed using \code{lgp_predict}.
#' @param uncertainty Either "none", "ribbon" or "errorbar".
#' @param test_data Test data set.
#' @param n_sds number of standard deviations for the uncertainty band width
#' @param data_marker pch for data points
#' @return a ggplot object
plot_posterior_y <- function(fit, PRED, 
                             uncertainty = "ribbon", 
                             test_data = NULL,
                             data_marker = 16,
                             n_sds = 2){
  
  if(uncertainty=="ribbon"){
    h <- plot_posterior_predictions(fit, mode = "predictive", PRED = PRED,
                                    test_data = test_data, n_sds = n_sds,
                                    data_marker = data_marker)
  }else if(uncertainty=="errorbar"){
    h <- plot_posterior_predictions(fit, mode = "predictive", PRED = PRED,
                                    error_bar = TRUE, test_data = test_data,
                                    n_sds = n_sds, data_marker = data_marker)
  }else{
    h <- plot_posterior_predictions(fit, mode = "predictive", PRED = PRED, 
                                    plot_uncertainty = FALSE, test_data = test_data,
                                    n_sds = n_sds, data_marker = data_marker) 
  }
  return(h)
}


#' Plot posterior of f or predictive distribution for y
#' 
#' @param fit An object of class \code{lgpfit}.
#' @param mode Must be either "posterior" or "predictive".
#' @param PRED Predictions computed using \code{lgp_predict}.
#' @param color_scheme Name of bayesplot color scheme or a list with fieds 'dark' and 'light'.
#' @param alpha Ribbon fill opacity.
#' @param alpha_line Line opacity.
#' @param alpha2 alpha of t_onset density
#' @param plot_uncertainty Should an uncertainty ribbon be plotted?
#' @param title optional prefix to plot title
#' @param ylim y axis limits
#' @param plot_obs_onset should the observed disease onset/initiation time be plotted
#' by a vertical line
#' @param plot_t_effect_samples should a distribution of sampled effect times be plotted
#' @param color_scheme_t_effect color scheme name for effect time density plotting
#' @param ypos_dens y-position of the density plot
#' @param test_data Test data frame
#' @param color_test test point color
#' @param pch_test test point marker
#' @param size_test test point size
#' @param error_bar should uncertainty be plotted using error bars instead of a ribbon
#' @param n_sds number of standard deviations for the uncertainty band width
#' @param reference_times reference onset times 
#' @param post_t_effect_stat statistic computed from effect time samples (mean or median)
#' @param ons_linetypes onset line types
#' @param ons_linecolors onset line colors
#' @param data_marker data marker type
#' @param data_color data marker color
#' @param original_y_scale should the predictions be scaled back to original data scale
#' @return a ggplot object
plot_posterior_predictions <- function(fit, 
                                       mode,
                                       PRED               = NULL, 
                                       color_scheme       = "red",
                                       color_scheme_t_effect = "gray",
                                       alpha              = 0.5, 
                                       alpha_line         = 1,
                                       alpha2             = 0.5,
                                       plot_uncertainty   = TRUE,
                                       title              = NULL,
                                       ylim               = NULL,
                                       plot_obs_onset     = FALSE,
                                       plot_t_effect_samples = FALSE,
                                       ypos_dens          = NULL,
                                       test_data          = NULL,
                                       color_test         = "deepskyblue2",
                                       pch_test           = 21,
                                       size_test          = 2,
                                       error_bar          = FALSE,
                                       n_sds              = 2,
                                       reference_times   = NULL,
                                       post_t_effect_stat = "none",
                                       original_y_scale   = TRUE,
                                       data_color         = "black",
                                       data_marker        = 21,
                                       ons_linetypes      = c(1,2,3),
                                       ons_linecolors     = c("black", "red", "gray50"))
{
  
  # Input checks and options
  OPT <- plot_predictions_options(fit, color_scheme,
                                  original_y_scale, PRED, test_data, 
                                  color_scheme_t_effect, mode, n_sds)
  DF      <- OPT$DF
  idvar   <- OPT$idvar
  timevar <- OPT$timevar
  respvar <- OPT$respvar
  model   <- fit@model
  info    <- model@info
  
  # Function for formatting the IDs
  formatter <- function(x){formatC(x, width = 2, format = "d", flag = "0")}
  id_pred   <- as.numeric(as.vector(DF[[idvar]]))
  DF$facet_var <- formatter(id_pred)
  # Set title
  if(mode == "posterior"){
    if(info$sample_F){
      ptitle <- " " 
    }else{
      ptitle <- " "
    }
  }else{
    ptitle <- "Posterior predictive distribution"
  }
  ylab <- respvar
  
  # Create ggplot object
  if(info$sample_F){
    group_var <- "idx"
    y_var     <- "pred"
  }else{
    group_var <- idvar
    y_var     <- "mu"
  }
  
  if(error_bar){
    h <- ggplot2::ggplot(DF, ggplot2::aes_string(x = timevar, y = y_var))
  }else{
    h <- ggplot2::ggplot(DF, ggplot2::aes_string(x = timevar, y = y_var, group = group_var))
  }
  
  # Edit plot title
  if(!is.null(title)){
    ptitle <- paste(title, ": ", ptitle, sep="")
  }
  h <- h + ggplot2::ggtitle(label = ptitle)
  h <- h + ggplot2::labs(y = ylab)
  
  # Plot uncertainty ribbon or error bars
  if(plot_uncertainty && !info$sample_F){
    if(error_bar){
      age_train <- as.numeric(model@data[[timevar]])
      age_test  <- c(DF[[timevar]])
      trange    <- range(c(age_test, age_train))
      width     <- 0.035*(trange[2]-trange[1])
      h <- h + ggplot2::geom_errorbar(ggplot2::aes_string(ymin = 'lower', ymax = 'upper'),
                                      color = OPT$hlcolor,
                                      width = width)
    }else{
      h <- h + ggplot2::geom_ribbon(ggplot2::aes_string(ymin = 'lower', ymax = 'upper'),
                                    fill  = OPT$fill, 
                                    alpha = alpha)
    }
  }
  
  # Plot mean prediction
  if(error_bar){
    h <- h + ggplot2::geom_point(color = OPT$hlcolor)
  }else{
    h <- h + ggplot2::geom_line(color = OPT$linecolor, alpha = alpha_line)
  }
  
  # Add onsets
  h <- plot_predictions_add_onsets(fit, h, plot_obs_onset, plot_t_effect_samples,
                                   idvar, timevar, ypos_dens, OPT$cs_onset,
                                   reference_times, post_t_effect_stat,
                                   ons_linetypes, ons_linecolors)
  
  # Plot also the data
  dat  <- model@data
  X_id <- as.numeric(dat[[idvar]])
  df   <- data.frame(idvar   = X_id, 
                     timevar = as.numeric(dat[[timevar]]),
                     y       = as.numeric(dat[[respvar]]),
                     facet_var = formatter(X_id))
  colnames(df)[1] <- idvar
  colnames(df)[2] <- timevar
  colnames(df)[3] <- respvar
  h <- h + ggplot2::geom_point(data = df, mapping = ggplot2::aes_string(
    x = timevar,
    y = respvar),
    color = data_color,
    pch = data_marker,
    inherit.aes = FALSE
  )
  
  
  # Plot test data
  if(!is.null(test_data)){
    id_test    <- formatter(test_data[[idvar]])
    age_test   <- test_data[[timevar]]
    y_test     <- test_data[[respvar]]
    point.data <- data.frame(xxx       = age_test,
                             yyy       = y_test,
                             facet_var = id_test)
    if(pch_test > 20){
      edgecol <- "black"
    }else{
      edgecol <- color_test
    }
    h <- h + ggplot2::geom_point(mapping  = ggplot2::aes_string(x = "xxx", y = "yyy"), 
                                 na.rm    = TRUE,
                                 inherit.aes = FALSE,
                                 data     = point.data, 
                                 color    = edgecol,
                                 pch      = pch_test,
                                 size     = size_test,
                                 fill     = color_test)
  }
  
  # Faceting, theme and limits
  h <- h + ggplot2::facet_wrap(. ~ facet_var)
  if(!is.null(ylim)){h <- h + ggplot2::coord_cartesian(ylim = ylim)}
  
  # Return modified ggplot object
  return(h)
}


#' Do input checks and set options for plotting predictions
#' @param fit An object of class \code{lgpfit}.
#' @param PRED Predictions computed using \code{lgp_predict}.
#' @param color_scheme Name of bayesplot color scheme.
#' @param color_scheme_t_effect Another color scheme.
#' @param original_y_scale Boolean value.
#' @param test_data test data
#' @param mode mode
#' @param n_sds number of standard deviations for the uncertainty band width
#' @return a list
plot_predictions_options <- function(fit, 
                                     color_scheme, 
                                     original_y_scale,
                                     PRED,
                                     test_data,
                                     color_scheme_t_effect,
                                     mode,
                                     n_sds){
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
  if(class(color_scheme_t_effect)=="character"){
    color_scheme_t_effect <- bayesplot::color_scheme_get(color_scheme_t_effect)
  }
  linecolor <- color_scheme$dark_highlight
  fill      <- color_scheme$mid
  hlcolor   <- color_scheme$mid_highlight
  
  # Create plot data frame
  if(is.null(PRED)){
    if(mode=="predictive"){ stop("PRED not given") }
    DF <- create_predictions_plot_df1(fit,  scale_f = original_y_scale,
                                      n_sds = n_sds)
  }else{
    if(info$sample_F){
      stop("Do not give predictions as input for a model where F was sampled!")
    }
    DF <- create_predictions_plot_df2(model, PRED, scale_f = original_y_scale, 
                                      mode = mode, n_sds = n_sds)
  }
  
  # Return
  ret <- list(linecolor = linecolor,
              fill      = fill,
              hlcolor   = hlcolor,
              DF        = DF,
              idvar     = idvar,
              timevar   = timevar,
              respvar   = respvar,
              cs_onset  = color_scheme_t_effect)
  
  return(ret)
}


#' Add disease onset / effect times to predictions plot
#' @description NOTE: currently assumes that diseased individuals come first.
#' @param fit An object of class \code{lgpfit}.
#' @param h a ggplot object
#' @param plot_obs_onset a boolean value
#' @param plot_t_effect_samples a boolean value
#' @param idvar id variable name
#' @param timevar time variable name
#' @param ypos_dens y position of the estimated onset density
#' @param color_scheme_t_effect color scheme
#' @param reference_times reference onset times 
#' @param linetypes onset line types
#' @param linecolors onset line colors
#' @param post_t_effect_stat statistic computed from effect time samples
#' @param alpha2 alpha parameter
#' @return a modified ggplot object
plot_predictions_add_onsets <- function(fit, h, 
                                        plot_obs_onset, 
                                        plot_t_effect_samples,
                                        idvar,
                                        timevar,
                                        ypos_dens,
                                        color_scheme_t_effect,
                                        reference_times,
                                        post_t_effect_stat,
                                        linetypes = c(1,2,3),
                                        linecolors = c("black", "red", "gray50"),
                                        alpha2 = 1)
{
  model <- fit@model
  D     <- model@stan_dat$D
  if(D[3]==1 && plot_obs_onset){
    df    <- model@data
    davar <- model@info$varInfo$disAge_variable
    t_ons <- get_obs_onset_times(id = df[[idvar]], age = df[[timevar]], disAge = df[[davar]])
    vline.data <- data.frame(zzz = t_ons, facet_var = names(t_ons))
    
    
    # First possible reference onsets
    if(!is.null(reference_times)){
      l1 <- length(reference_times)
      l2 <- length(t_ons)
      if(l1!=l2){
        stop("invalid length of reference_times (", l1, "), must be ", l2)
      }
      refline.data <- data.frame(zzz = reference_times, facet_var = names(t_ons))
      h <- h + ggplot2::geom_vline(ggplot2::aes_string(xintercept = "zzz"), 
                                   na.rm    = TRUE,
                                   data     = refline.data, 
                                   color    = linecolors[1],
                                   linetype = linetypes[1])
    }
    
    # Then observed onsets
    h <- h + ggplot2::geom_vline(ggplot2::aes_string(xintercept = "zzz"), 
                                 na.rm    = TRUE,
                                 data     = vline.data, 
                                 color    = linecolors[2],
                                 linetype = linetypes[2])
  }
  
  # Plot effect time samples
  UNCRT <- model@stan_dat$UNCRT
  if(UNCRT==1 && D[3]==1 & plot_t_effect_samples){
    
    T_smp <- extract_t_effect_samples(fit)
    n_smp <- dim(T_smp)[1]
    t_smp <- as.numeric(T_smp)
    df    <- model@data
    yvar  <- model@info$varInfo$response_variable
    ystd  <- stats::sd(df[[yvar]])
    dens_scale <- ystd
    if(is.null(ypos_dens)){
      ypos_dens <- min(df[[yvar]]) - 1.5*dens_scale
    }
    
    # Plot effect time mean or median
    cid_str   <- colnames(T_smp)
    if(post_t_effect_stat=="mean"){
      statistic <- colMeans(T_smp)
    }else if(post_t_effect_stat=="median"){
      statistic <- apply(T_smp, 2, stats::median)
    }else{
      statistic <- NULL
    }
    if(!is.null(statistic)){
      nams <- names(t_ons)
      i_good <- which(nams %in% cid_str)
      nams <- nams[i_good]
      statline.data <- data.frame(zzz = statistic, facet_var = nams)
      h <- h + ggplot2::geom_vline(ggplot2::aes_string(xintercept = "zzz"), 
                                   na.rm    = TRUE,
                                   data     = statline.data, 
                                   color    = linecolors[3],
                                   linetype = linetypes[3])
    }
    
    # Plot kernel density estimate of posterior of effect time
    fv        <- as.factor(rep(cid_str, each = n_smp))
    dens.data <- data.frame(zzz = t_smp, facet_var = fv, id = fv)
    h <- h + ggplot2::geom_density(ggplot2::aes_string(x = "zzz", 
                                                       y = paste(dens_scale,
                                                                 "*..scaled..",sep="")), 
                                   na.rm    = TRUE,
                                   data     = dens.data, 
                                   color    = color_scheme_t_effect$mid_highlight,
                                   fill     = color_scheme_t_effect$mid,
                                   alpha    = alpha2,
                                   inherit.aes = F, 
                                   position = ggplot2::position_nudge(x=0, y=ypos_dens),
                                   trim     = TRUE)
  }
  return(h)
}


#' Create a plotting data frame for ggplot
#' 
#' @description A helper function for \code{plot_predictions}.
#' @param fit An object of class \code{lgpfit}.
#' @param n_sds number of standard deviations for the uncertainty band width
#' @param scale_f Should the predictions be scaled back to the original data scale?
#' @return a data frame
create_predictions_plot_df1 <- function(fit, scale_f = TRUE, n_sds){
  
  model <- fit@model
  info  <- model@info
  SCL   <- model@scalings
  TSCL  <- model@scalings$TSCL
  YSCL  <- model@scalings$YSCL
  
  idvar   <- info$varInfo$id_variable
  timevar <- info$varInfo$time_variable
  
  X    <- t(model@stan_dat$X)
  id   <- X[,1]
  age  <- X[,2]
  age  <- TSCL$fun_inv(age)
  LH   <- model@stan_dat$LH
  
  PRED  <- get_predicted(fit)
  if(!info$sample_F){
    
    mu   <- as.numeric(PRED$pred)
    std  <- as.numeric(PRED$std)
    
    # Scale some things
    upper <- mu + n_sds*std
    lower <- mu - n_sds*std
    if(scale_f){
      mu    <- YSCL$fun_inv(mu)
      upper <- YSCL$fun_inv(upper)
      lower <- YSCL$fun_inv(lower)
    }
    
    # Create data frame
    DF    <- data.frame(id    = as.factor(id),
                        age   = as.numeric(age),
                        mu    = as.numeric(mu),
                        upper = as.numeric(upper),
                        lower = as.numeric(lower)
    )
    colnames(DF)[1] <- idvar
    colnames(DF)[2] <- timevar
    
  }else{
    
    G_smp <- PRED$G_smp
    n_smp <- dim(G_smp)[1]
    n_tot <- dim(G_smp)[2]
    
    age   <- rep(age, n_smp)
    id    <- rep(id, n_smp)
    s_idx <- rep(1:n_smp, each = n_tot)
    if(scale_f && LH==1){
      G_smp  <- YSCL$fun_inv(G_smp)
    }
    g <- as.vector(t(G_smp))
    
    # Create data frame
    DF    <- data.frame(id    = as.factor(id),
                        age   = as.numeric(age),
                        idx   = as.factor(s_idx),
                        pred  = as.numeric(g)
    )
    colnames(DF)[1] <- idvar
    colnames(DF)[2] <- timevar
  }
  
  return(DF)
}


#' Create a plotting data frame for ggplot
#' 
#' @description A helper function for \code{plot_predictions}.
#' @param model An object of class \code{lgpmodel}.
#' @param PRED Predictions computed using \code{lgp_predict}.
#' @param scale_f Should the predictions be scaled back to the original data scale?
#' @param mode mode
#' @param n_sds number of standard deviations for the uncertainty band width
#' @return a data frame
create_predictions_plot_df2 <- function(model, 
                                        PRED, 
                                        scale_f = TRUE,
                                        mode,
                                        n_sds){
  
  LIST    <- PRED$LIST
  X_test  <- PRED$X_test_scaled
  L       <- length(LIST)
  if(L > 1){
    cat("* Averaging the predictions over", L, "samples.\n")
    pred <- average_predictions(LIST)
  }else{
    pred <- LIST[[1]]
  }
  TSCL    <- model@scalings$TSCL
  YSCL    <- model@scalings$YSCL
  info    <- model@info
  idvar   <- info$varInfo$id_variable
  timevar <- info$varInfo$time_variable
  
  ID  <- X_test[,1]
  AGE <- X_test[,2]
  AGE <- TSCL$fun_inv(AGE)
  n   <- length(ID)
  cn  <- model@info$component_names
  
  # Create data frame columns
  MU  <- pred$mu_f
  if(mode=="predictive"){
    S2 <- pred$s2_y
  }else if(mode=="posterior"){
    S2  <- pred$s2_f
  }else{
    stop("mode must be 'posterior' or 'predictive' !")
  }
  id  <- ID
  age <- AGE
  mu  <- as.numeric(MU)
  std <- sqrt(as.numeric(S2))
  
  # Scale some things
  upper <- mu + n_sds*std
  lower <- mu - n_sds*std
  if(scale_f){
    mu    <- YSCL$fun_inv(mu)
    upper <- YSCL$fun_inv(upper)
    lower <- YSCL$fun_inv(lower)
  }
  
  # Create data frame
  DF    <- data.frame(id    = as.factor(id),
                      age   = as.numeric(age),
                      mu    = as.numeric(mu),
                      upper = as.numeric(upper),
                      lower = as.numeric(lower)
  )
  colnames(DF)[1] <- idvar
  colnames(DF)[2] <- timevar
  return(DF)
}

