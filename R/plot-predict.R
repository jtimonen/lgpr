#' Plot posterior of f
#' 
#' @export
#' @param fit An object of class \code{lgpfit}.
#' @param componentwise A boolean value.
#' @param PRED Predictions computed using \code{lgp_predict}.
#' @return a ggplot object
plot_posterior_f <- function(fit, componentwise = FALSE, PRED = NULL){
  if(class(componentwise)!="logical"){
    stop("componentwise must be a logical argument!")
  }
  h <- plot_predictions(fit, componentwise = componentwise, 
                        mode = "posterior", PRED)
  return(h)
}

#' Plot posterior predictive distribution
#' 
#' @export
#' @param fit An object of class \code{lgpfit}.
#' @param PRED Predictions computed using \code{lgp_predict}.
#' @param uncertainty Either "none", "ribbon" or "errorbar".
#' @param test_data Test data set.
#' @return a ggplot object
plot_posterior_y <- function(fit, PRED, uncertainty = "ribbon", test_data = NULL){
  if(is.null(PRED)){
    stop("You must specify the PRED argument!")
  }
  if(uncertainty=="ribbon"){
    h <- plot_predictions(fit, mode = "predictive", PRED = PRED,
                          test_data = test_data)
  }else if(uncertainty=="errorbar"){
    h <- plot_predictions(fit, mode = "predictive", PRED = PRED, 
                          error_bar = TRUE, test_data = test_data)
  }else{
    h <- plot_predictions(fit, mode = "predictive", PRED = PRED, 
                          plot_uncertainty = FALSE, test_data = test_data) 
  }
  return(h)
}


#' Plot posterior of f or predictive distribution for y
#' 
#' @description NOTE: currently plotting the onsets assumes that the case individuals come first!
#' @param fit An object of class \code{lgpfit}.
#' @param mode Must be either "posterior" or "predictive".
#' @param PRED Predictions computed using \code{lgp_predict}.
#' @param componentwise Should the predictions be plotted componentwise?
#' @param color_scheme Name of bayesplot color scheme or a list with fieds 'dark' and 'light'.
#' @param alpha Ribbon fill opacity.
#' @param alpha_line Line opacity
#' @param plot_uncertainty Should an uncertainty ribbon be plotted?
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
#' @param error_bar should uncertainty be plotted using error bars instead of a ribbon
#' @return a ggplot object
plot_predictions <- function(fit, 
                             mode,
                             PRED               = NULL, 
                             componentwise      = FALSE,
                             color_scheme       = "red",
                             alpha              = 0.5-0.4*as.numeric(componentwise), 
                             plot_uncertainty   = TRUE,
                             title              = NULL,
                             theme              = ggplot2::theme_gray(),
                             original_y_scale   = TRUE,
                             alpha_line         = 1,
                             ylim               = NULL,
                             plot_obs_onset     = FALSE,
                             alpha2             = 0.5,
                             color_scheme_onset = "gray",
                             plot_onset_samples = FALSE,
                             ypos_dens          = NULL,
                             test_data          = NULL,
                             color_test         = "deepskyblue2",
                             pch_test           = 21,
                             size_test          = 2,
                             error_bar          = FALSE)
{
  
  # Input checks and options
  OPT <- plot_predictions_options(fit, color_scheme, componentwise,
                                  original_y_scale, PRED, test_data, color_scheme_onset, mode)
  DF      <- OPT$DF
  idvar   <- OPT$idvar
  timevar <- OPT$timevar
  respvar <- OPT$respvar
  model   <- fit@model
  info    <- model@info
  
  # Facet variable and labels
  X_id <- as.numeric(model@stan_dat$X[1,])
  if(componentwise){
    facet_var    <- "component"
    DF$facet_var <- DF[[facet_var]]
    ptitle       <- "Posteriors of components of f"
    ylab         <- " "
  }else{
    facet_var    <- idvar
    id_train     <- unique(X_id)
    DF$facet_var <- DF[[facet_var]]
    id_test      <- unique(as.numeric(DF$facet_var))
    all_ids      <- union(id_train, id_test)
    all_ids      <- sort(all_ids, decreasing = FALSE)
    DF$facet_var <- as.numeric(factor(DF$facet_var, levels = all_ids))
    if(mode == "posterior"){
      ptitle <- "Posterior distribution of f"
    }else{
      ptitle <- "Posterior predictive distribution"
    }
    
    ylab <- respvar
  }
  
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
  
  # Plot also the data
  if(!componentwise){
    dat  <- model@data
    df   <- data.frame(idvar   = as.factor(X_id), 
                       timevar = as.numeric(dat[[timevar]]),
                       y       = as.numeric(dat[[respvar]]))
    colnames(df)[1] <- idvar
    colnames(df)[2] <- timevar
    colnames(df)[3] <- respvar
    df$facet_var    <- df[[facet_var]]
    h <- h + ggplot2::geom_point(data = df, mapping = ggplot2::aes_string(
      x = timevar,
      y = respvar)
    )
  }
  
  # Faceting, theme and limits
  h <- h + ggplot2::facet_wrap(. ~ facet_var)
  h <- h + theme
  if(!is.null(ylim)){h <- h + ggplot2::coord_cartesian(ylim = ylim)}
  
  # Add onsets
  if(!componentwise){
    h <- plot_predictions_add_onsets(fit, h, plot_obs_onset, plot_onset_samples,
                                     idvar, timevar, ypos_dens, OPT$cs_onset)
  }
  
  # Plot test data
  if(!is.null(test_data)){
    id_test    <- test_data[[idvar]]
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
                                 inherit.aes = F,
                                 data     = point.data, 
                                 color    = edgecol,
                                 pch      = pch_test,
                                 size     = size_test,
                                 fill     = color_test)
  }
  
  # Return modified ggplot object
  return(h)
}


#' Do input checks and set options for plotting predictions
#' @param fit An object of class \code{lgpfit}.
#' @param PRED Predictions computed using \code{lgp_predict}.
#' @param componentwise Should the predictions be plotted componentwise?
#' @param color_scheme Name of bayesplot color scheme.
#' @param color_scheme_onset Another color scheme.
#' @param original_y_scale Boolean value.
#' @param test_data test data
#' @param mode mode
#' @return a list
plot_predictions_options <- function(fit, 
                                     color_scheme, 
                                     componentwise,
                                     original_y_scale,
                                     PRED,
                                     test_data,
                                     color_scheme_onset,
                                     mode){
  # Check input correctness
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  model <- fit@model
  if(class(model)!="lgpmodel") stop("Class of 'fit@model' must be 'lgpmodel'!")
  info    <- model@info
  idvar   <- info$varInfo$id_variable
  timevar <- info$varInfo$time_variable
  respvar <- info$varInfo$response_variable
  if(componentwise && !is.null(test_data)){
    stop("componentwise must be FALSE if test data is given!")
  }
  
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
    if(mode=="predictive"){
      stop("Plotting predictive distribution not implemented in this case!")
    }
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
                                      componentwise = componentwise,
                                      mode = mode)
  }
  
  # Return
  ret <- list(linecolor = linecolor,
              fill      = fill,
              hlcolor   = hlcolor,
              DF        = DF,
              idvar     = idvar,
              timevar   = timevar,
              respvar   = respvar,
              cs_onset  = color_scheme_onset)
  
  return(ret)
}


#' Add disease onsets to predictions plot
#' @description NOTE: currently assumes that diseased individuals come first.
#' @param fit An object of class \code{lgpfit}.
#' @param h a ggplot object
#' @param plot_obs_onset a boolean value
#' @param plot_onset_samples a boolean value
#' @param idvar id variable name
#' @param timevar time variable name
#' @param ypos_dens y position of the estimated onset density
#' @param color_scheme_onset color scheme
#' @param alpha2 alpha parameter
#' @return a modified ggplot object
plot_predictions_add_onsets <- function(fit, h, 
                                        plot_obs_onset, 
                                        plot_onset_samples,
                                        idvar,
                                        timevar,
                                        ypos_dens,
                                        color_scheme_onset,
                                        alpha2 = 1)
{
  # Plot observed onsets
  model <- fit@model
  D     <- model@stan_dat$D
  if(D[3]==1 && plot_obs_onset){
    df    <- model@data
    davar <- model@info$varInfo$disAge_variable
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
  if(UNCRT==1 && D[3]==1 & plot_onset_samples){
    
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
  return(h)
}


#' Create a plotting data frame for ggplot
#' 
#' @description A helper function for \code{plot_predictions}.
#' @param fit An object of class \code{lgpfit}.
#' @param scale_f Should the predictions be scaled back to the original data scale?
#' @return a data frame
create_predictions_plot_df1 <- function(fit, scale_f = TRUE){
  
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
    upper <- mu + 2*std
    lower <- mu - 2*std
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
#' @param componentwise Should the predictions be plotted componentwise?
#' @param mode mode
#' @return a data frame
create_predictions_plot_df2 <- function(model, PRED, scale_f = TRUE, componentwise = FALSE, mode){
  
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
  if(componentwise){
    if(mode=="predictive"){
      stop("cannot compute posterior predictive distribution when plotting componentwise")
    }
    MU        <- pred$mu_cmp
    S2        <- pred$s2_cmp
    d         <- dim(MU)[2]
    id        <- rep(ID, d)
    age       <- rep(AGE, d)
    component <- rep(cn, each = n)
  }else{
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
  }
  mu  <- as.numeric(MU)
  std <- sqrt(as.numeric(S2))
  
  # Scale some things
  upper <- mu + 2*std
  lower <- mu - 2*std
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
  if(componentwise){
    DF$component <- as.factor(component)
  }
  return(DF)
}

