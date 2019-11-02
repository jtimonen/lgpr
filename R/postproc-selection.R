
#' Selection of relevant components
#'
#' @export
#' @param object An object of class \code{lgpfit}.
#' @param threshold Threshold for relevance sum. 
#' Must be a value between 0 and 1.
#' @return A named list
selection <- function(object, threshold = 0.95)
{
  if(class(object)!="lgpfit") stop("Class of 'object' must be 'lgpfit'!")
  if(threshold > 1 || threshold < 0) {
    stop("'threshold' must be between 0 and 1!")
  }
  
  # Get relevances
  relevances  <- object@relevances
  rel_avg     <- relevances$average
  rel_smp     <- relevances$samples
  
  # Strict selection using average relevances
  info     <- object@model@info
  rel_avg  <- object@relevances$average
  i_sel    <- selection_fixed_threshold(rel_avg, threshold)
  nam      <- c(info$component_names, "noise")
  selected <- nam[i_sel]
  
  # Selection probabilities
  prob     <- selection_prob_fixed_threshold(rel_smp, threshold)
  
  res <- list(selected  = selected,
              prob      = prob,
              threshold = threshold)
  return(res)
}


#' Probabilistic selection of relevant components
#'
#' @export
#' @param object An object of class \code{lgpfit}.
#' @param p a function defining a density over interval [0,1]
#' @param h discretization parameter for computing a quadrature
#' @return Selection probabilities for each component
selection_prob <- function(object, 
                           p = function(x){stats::dbeta(x, 100, 5)},
                             h = 0.01)
{
  if(class(object)!="lgpfit") stop("Class of 'object' must be 'lgpfit'!")
  if(!is.function(p)){stop('p must be a function')}
  
  # Get relevances
  rel_smp     <- object@relevances$samples
  
  info <- object@model@info
  H    <- seq(0, 1, by = h)
  L    <- length(H)
  P    <- rep(0, L)
  prob <- selection_prob_fixed_threshold(rel_smp, threshold = 0)
  D    <- length(prob) - 1 
  prob <- prob[1:D]
  nam  <- names(prob)
  PROB <- matrix(0, L, D)
  for(i in 1:L){
    P[i] <- p(H[i])
    prob <- selection_prob_fixed_threshold(rel_smp, threshold = H[i])
    PROB[i,] <- prob[1:D]
  }
  colnames(PROB) <- nam
  P <- matrix(rep(P, D), L, D, byrow = FALSE)
  res <- h * colSums(PROB * P)
  plt <- selection_prob_plot(PROB, H, P)
  ret <- list(prob = res, plot = plt)
  return(ret)
}


#' Selection probabilities using a fixed threshold
#'
#' @param relevances The \code{relevances$samples} slot of an 
#' \code{lgpfit} object.
#' @param threshold value between 0 and 1
#' @return proportion of times each component was selected
selection_prob_fixed_threshold <- function(relevances, threshold)
{
  n_smp <- dim(relevances)[1]
  n_cmp <- dim(relevances)[2]-1
  names <- colnames(relevances)
  sel   <- matrix(0, n_smp, n_cmp+1)
  for(i in 1:n_smp){
    i_sel <- selection_fixed_threshold(relevances[i,], threshold)
    sel[i,i_sel] <- 1
  }
  colnames(sel) <- names
  return(colMeans(sel))
}


#' Select relevant components
#'
#' @param rel a named vector of component relevances
#' @param threshold value between 0 and 1
#' @return indices of selected components (including "noise" always)
selection_fixed_threshold <- function(rel, threshold)
{
  n_cmp      <- length(rel) - 1
  i_noise    <- n_cmp + 1
  p_noise    <- rel[i_noise]
  rel        <- as.numeric(rel[1:n_cmp])
  s          <- sort(rel, decreasing = TRUE, index.return = TRUE)
  rel        <- s$x
  if(p_noise >= threshold){
    i_sel <- i_noise
    return(i_sel)
  }else{
    for(j in 1:length(rel)){
      h <- p_noise + sum(rel[1:j])
      if(h >= threshold){
        i_sel <- c(i_noise, s$ix[1:j])
        return(i_sel)
      }
    }
    i_sel <- c(1:(n_cmp+1))
    return(i_sel)
  }

}


#' Plot of probabilistic selection of relevant components
#'
#' @param PROB computed probabilities at points H
#' @param H a grid on interval [0,1]
#' @param P threshold probability distribution evaluated at H
#' @return a ggplot object
selection_prob_plot <- function(PROB, H, P)
{
  nam <- colnames(PROB)
  n   <- dim(PROB)[1]
  d   <- dim(PROB)[2]
  pr  <- as.numeric(PROB)
  cp  <- as.factor(rep(nam, each = n))
  h   <- rep(H, d)
  p   <- rep(P, d)
  df  <- data.frame(h, pr, cp)
  colnames(df) <- c("Threshold", "Probability", "Component")
  plt <- ggplot2::ggplot(df, 
        ggplot2::aes_string(x = 'Threshold',
                            y = 'Probability',
                            group = 'Component',
                             color = 'Component'))
  plt <- plt + ggplot2::geom_line() + ggplot2::theme_linedraw()
  plt <- plt + ggplot2::xlim(c(0,1)) + ggplot2::ylim(c(0,1))
  return(plt)
}