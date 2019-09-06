#' A spaghetti plot of longitudinal data.
#' @export
#' @param data A data frame.
#' @param highlight Name of a covariate to be highlighted with color, or id of a subject
#' to be highlighted.
#' @param id_variable Name of id variable.
#' @param time_variable Name of time variable.
#' @param response Name of the response variable.
#' @param psize point size
#' @param lwd line width
#' @param title additional string added to title
#' @return a ggplot object
plot_data <- function(data, 
                      highlight     = NULL,
                      response      = "y",
                      id_variable   = "id",
                      time_variable = "age",
                      psize         = 2,
                      lwd           = 0.5,
                      title         = NULL)
{
  
  if(!is.null(highlight)){
    if(is.numeric(highlight)){
      p <- plot_data_hl_individual(data, highlight, 
                            response, id_variable, time_variable, 
                            psize, lwd)
    }else{
      if(highlight == "disease"){
        highlight <- "diseaseAge"
      }
      if(highlight == "diseaseAge"){
        p <- plot_data_hl_disease(data, "diseaseAge", 
                                  response, id_variable, time_variable, 
                                  psize, lwd)
      }else{
        p <- plot_data_hl_cat(data, highlight, 
                              response, id_variable, time_variable, 
                              psize, lwd)
      }
    }
  }else{
    p <- plot_data_plain(data, response, id_variable, time_variable, 
                         psize, lwd)
  }
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


#' A spaghetti plot of longitudinal data without highlighting.
#' @param data A data frame.
#' @param id_variable Name of id variable.
#' @param time_variable Name of time variable.
#' @param response Name of the response variable.
#' @param psize point size
#' @param lwd line width
#' @return a ggplot object
plot_data_plain <- function(data, 
                            response      = "y",
                            id_variable   = "id",
                            time_variable = "age",
                            psize         = 2,
                            lwd           = 0.5)
{
  data$id  <- as.factor(data[[id_variable]])
  iresp    <- which(colnames(data)==response)
  colnames(data)[iresp] <- "Response"
  
  ord  <- order(data$id, data[[time_variable]])
  data <- data[ord, ]
  
  # Extend data frame
  data <- create_data_plot_df(data, NULL, NULL, NULL)
  
  # Create plot
  p <- ggplot2::ggplot(data = data, ggplot2::aes_string(x = time_variable, 
                                                        y = 'Response', group = "id")) +
    ggplot2::geom_line(size = lwd)
  if(psize > 0){
    p <- p + ggplot2::geom_point(size = psize)
  }
  return(p)
}

#' A spaghetti plot of longitudinal data, highlighting a categorical covariate.
#' @param data A data frame.
#' @param highlight Name of a categorical covariate to be highlighted with color.
#' @param id_variable Name of id variable.
#' @param time_variable Name of time variable.
#' @param response Name of the response variable.
#' @param psize point size
#' @param lwd line width
#' @return a ggplot object
plot_data_hl_cat <- function(data, 
                             highlight     = NULL,
                             response      = "y",
                             id_variable   = "id",
                             time_variable = "age",
                             psize         = 2,
                             lwd           = 0.5)
{
  data$id  <- as.factor(data[[id_variable]])
  iresp    <- which(colnames(data)==response)
  colnames(data)[iresp] <- "Response"
  
  ord  <- order(data[[highlight]], data$id, data[[time_variable]])
  data <- data[ord, ]
  
  # Extend data frame
  data <- create_data_plot_df(data, highlight, NULL, NULL)
  
  # Create plot
  p <- ggplot2::ggplot(data = data, ggplot2::aes_string(x = time_variable, 
                                                        y = 'Response', group = "id"))
  
  # Highlight categories
  p <- p + ggplot2::labs(colour = highlight)
  p <- p + ggplot2::scale_color_brewer(palette = "Set1")
  p <- p + ggplot2::geom_line(ggplot2::aes_string(color = 'Z1'), size = lwd)
  if(psize > 0){
    p <- p + ggplot2::geom_point(ggplot2::aes_string(color = 'Z1'), size = psize)
  }
  return(p)
}


#' A spaghetti plot of longitudinal data, highlighting a continuous covariate.
#' @param data A data frame.
#' @param highlight Name of a continuous covariate to be highlighted with color.
#' @param id_variable Name of id variable.
#' @param time_variable Name of time variable.
#' @param response Name of the response variable.
#' @param psize point size
#' @param lwd line width
#' @param colgrad color gradient
#' @return a ggplot object
plot_data_hl_cont <- function(data, 
                              highlight     = NULL,
                              response      = "y",
                              id_variable   = "id",
                              time_variable = "age",
                              psize         = 2,
                              lwd           = 0.5,
                              colgrad       = ggplot2::scale_colour_gradient2())
{
  
  data$id  <- as.factor(data[[id_variable]])
  iresp    <- which(colnames(data)==response)
  colnames(data)[iresp] <- "Response"
  
  # Extend data frame
  data <- create_data_plot_df(data, NULL, NULL, highlight)
  
  # Create plot
  p <- ggplot2::ggplot(data = data, ggplot2::aes_string(x = time_variable, 
                                                        y = 'Response', group = "id"))
  
  # Highligh continuous
  p <- p + ggplot2::geom_line(ggplot2::aes_string(color = 'Value'), size = lwd)
  p <- p + colgrad + ggplot2::labs(colour = highlight)
  if(psize > 0){
    p <- p + ggplot2::geom_point(ggplot2::aes_string(color = 'Z1'), size = psize)
  }
  return(p)
}



#' A spaghetti plot of longitudinal data, highlighting based on disease group.
#' @param data A data frame.
#' @param highlight Name of the disease-related age variable.
#' @param id_variable Name of id variable.
#' @param time_variable Name of time variable.
#' @param response Name of the response variable.
#' @param psize point size
#' @param lwd line width
#' @return a ggplot object
plot_data_hl_disease <- function(data, 
                                 highlight     = "diseaseAge",
                                 response      = "y",
                                 id_variable   = "id",
                                 time_variable = "age",
                                 psize         = 2,
                                 lwd           = 0.5)
{
  
  # Auxiliary group variable indicating diseased and healthy individuals
  if("diseaseAge" %in% colnames(data)){
    Group <- as.numeric(!is.nan(data[[highlight]]))
    data <- cbind(data, Group)
  }else{
    stop("diseaseAge must be a column of data!")
  }
  
  ord  <- order(data$Group, data$id, data[[time_variable]])
  data <- data[ord, ]
  strs <- c("control", "case")
  data$Group <- strs[data$Group+1]
  
  p <- plot_data_hl_cat(data, "Group", 
                        response, id_variable, time_variable, 
                        psize, lwd)
  return(p)
}


#' A spaghetti plot of longitudinal data, highlighting one individual.
#' @param data A data frame.
#' @param highlight Number indicating the individual to highlight.
#' @param id_variable Name of id variable.
#' @param time_variable Name of time variable.
#' @param response Name of the response variable.
#' @param psize point size
#' @param lwd line width
#' @return a ggplot object
plot_data_hl_individual <- function(data, 
                                    highlight     = 1,
                                    response      = "y",
                                    id_variable   = "id",
                                    time_variable = "age",
                                    psize         = 2,
                                    lwd           = 0.5)
{
  
  # Auxiliary group variable indicating diseased and healthy individuals
  is_hl <- 1 - as.numeric(data[[id_variable]]==highlight)
  ID <- c(paste("id =", highlight), "id = other")
  ID <- ID[is_hl+1]
  data <- cbind(data, ID)
  
  ord  <- order(data$ID, data$id, data[[time_variable]])
  data <- data[ord, ]
  
  p <- plot_data_hl_cat(data, "ID", 
                        response, id_variable, time_variable, 
                        psize, lwd)
  return(p)
}


#' Visualize the (average) inferred components evaluated at data points
#' @export
#' @param fit An object of class \code{lgpfit}.
#' @param corrected Should this plot the covariate-effect corrected components?
#' @param title optional prefix to plot title
#' @param sample_idx If given, only one sample will be plotted, else the average
#' components over all samples.
#' @param linealpha line alpha
#' @return a ggplot object
plot_components <- function(fit, 
                            corrected = TRUE, 
                            title = NULL, 
                            sample_idx = NULL,
                            linealpha = 0.6){
  
  
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
  subt <- "dotted lines connect an individual"
  h   <- ggplot2::ggplot(df, ggplot2::aes_string(x = timevar, 
                                                 y = 'f', 
                                                 group = idvar)) + 
    ggplot2::geom_line(linetype = 3, alpha = linealpha) + 
    ggplot2::geom_point() + 
    ggplot2::facet_wrap(. ~ component) +
    ggplot2::labs(y = yname) +
    ggplot2::ggtitle(label = ptitle, subtitle = subt)
  return(h)
}

