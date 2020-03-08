
#' Convert the Stan likelihood encoding to a string
#'
#' @param LH an integer
#' @return a string
likelihood_as_str <- function(LH){
  if(LH==1 || LH==-1){
    str <- "Gaussian"
  }else if(LH==0){
    str <- "none"
  }else if(LH==2){
    str <- "Poisson"
  }else if(LH==3){
    str <- "NB"
  }else if(LH==4){
    str <- "binomial"
  }else if(LH==5){
    str <- "ordinal"
  }else{
    str <- "Unknown likelihood"
  }
  return(str)
}


#' Convert likelihood string to Stan encoding
#'
#' @param likelihood a string
#' @return an integer
likelihood_as_int <- function(likelihood){
  likelihood <- tolower(likelihood)
  if(likelihood=="none"){
    LH <- 0
  }else if(likelihood=="gaussian"){
    LH <- 1
  }else if(likelihood=="poisson"){
    LH <- 2
  }else if(likelihood=="nb"){
    LH <- 3
  }else if(likelihood=="binomial"){
    LH <- 4
  }else if(likelihood=="ordinal"){
    LH <- 5
  }else{
    stop("likelihood must be either 'none', 'Gaussian', 
         'Poisson', 'NB', binomial', or 'ordinal'!")
  }
  return(LH)
}
