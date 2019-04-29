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
#' @return a data frame
create_predictions_plot_df2 <- function(model, PRED, scale_f = TRUE, componentwise = FALSE){
  MU      <- PRED$F_mu
  S2      <- PRED$F_var
  TSCL    <- model@scalings$TSCL
  YSCL    <- model@scalings$YSCL
  info    <- model@info
  idvar   <- info$varInfo$id_variable
  timevar <- info$varInfo$time_variable
  
  ID  <- PRED$X_test_scaled[,1]
  AGE <- PRED$X_test_scaled[,2]
  AGE <- TSCL$fun_inv(AGE)
  n   <- length(ID)
  d   <- dim(MU)[2] - 1
  cn  <- model@info$component_names
  cn  <- c(cn, "f")
  
  # Create data frame columns
  if(componentwise){
    id        <- rep(ID, d + 1)
    age       <- rep(AGE, d + 1)
    component <- rep(cn, each = n)
    mu        <- as.numeric(MU)
    std       <- sqrt(as.numeric(S2))
    
  }else{
    id  <- ID
    age <- AGE
    mu  <- as.numeric(MU[,d+1])
    std <- sqrt(as.numeric(S2[,d+1]))
  }
  
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
  DF    <- data.frame(id  = as.factor(id),
                      age = as.numeric(age),
                      component = as.factor(component),
                      value  = as.numeric(value)
  )
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

#' Simulated data column names  in a prettier form
#'
#' @param cn column names
#' @return names of model components
simdata_colnames_pretty <- function(cn){
  cn  <- gsub("\\.", ", ", cn)
  cn  <- paste("f_",1:length(cn),"(",cn,")", sep="")
  return(cn)
}
