#' Create a plotting data frame for ggplot
#' 
#' @description A helper function for \code{plot_simdata_by_component}.
#' @param simData An object created using \code{simulate_data}.
#' @return a data frame
create_simdata_plot_df <- function(simData){
  data <- simData$data
  comp <- simData$components
  id   <- data$id
  age  <- data$age
  n    <- length(id)
  d    <- dim(comp)[2] - 4
  CMP  <- comp[,1:(d+1)]
  cn   <- colnames(CMP)
  cn   <- cn[1:(length(cn)-1)]
  cn   <- simdata_colnames_pretty(cn)
  cn   <- c(cn, "f")
  id   <- rep(id, d+1)
  age  <- rep(age, d+1)
  component <- rep(cn, each = n)
  value <- as.numeric(as.matrix(CMP))
  DF    <- data.frame(id        = as.factor(id),
                      age       = as.numeric(age),
                      component = as.factor(component),
                      value     = as.numeric(value))
  return(DF)
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

#' Simulated data column names in a prettier form
#'
#' @param cn column names
#' @return names of model components
simdata_colnames_pretty <- function(cn){
  cn  <- gsub("\\.", ", ", cn)
  cn  <- paste("f_",1:length(cn),"(",cn,")", sep="")
  return(cn)
}
