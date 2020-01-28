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

#' Create a plotting data frame for ggplot
#' 
#' @description A helper function for \code{plot_data}.
#' @param data a data frame
#' @param hl_1 highlighting by color
#' @param hl_2 highlighting by linestyle
#' @param hl_cont highlighting continuous
#' @return an extended data frame
create_data_plot_df <- function(data, hl_1, hl_2, hl_cont){
  
  # Color highlight
  hl <- hl_1
  if(!is.null(hl)){
    if(hl == "diseaseAge" || hl == "disease" || hl == "group"){
      hl <- "Group"
    }
    ihl <- which(colnames(data)==hl)
    if(hl!="id"){
      data[[ihl]] <- as.factor(data[[ihl]])
      colnames(data)[ihl] <- "Z1" 
    }else{
      Z1 <- data[[ihl]]
      data <- cbind(data, Z1)
    }
  }
  
  # Linestyle highlight
  hl <- hl_2
  if(!is.null(hl)){
    if(hl == "diseaseAge" || hl == "disease" || hl == "group"){
      hl <- "Group"
    }
    ihl <- which(colnames(data)==hl)
    if(hl!="id"){
      data[[ihl]] <- as.factor(data[[ihl]])
      colnames(data)[ihl] <- "Z2" 
    }else{
      Z2 <- data[[ihl]]
      data <- cbind(data, Z2)
    }
  }
  
  # Continuous covariate highlight
  if(!is.null(hl_cont)){
    ihl <- which(colnames(data)==hl_cont)
    colnames(data)[ihl] <- "Value"
  }
  return(data)
}

