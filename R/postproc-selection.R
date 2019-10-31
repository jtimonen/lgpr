
#' Selection of relevant components using a threshold
#' for the total proportion of explained variance
#'
#' @export
#' @param object An object of class \code{lgpfit}.
#' @param threshold A value between 0 and 1
#' @return A named list
selection <- function(object, threshold = 0.95)
{
  if(class(object)!="lgpfit") stop("Class of 'object' must be 'lgpfit'!")
  if(threshold > 1 || threshold < 0) stop("'threshold' must be between 0 and 1!")
  info       <- object@model@info
  rel        <- object@relevances$average
  rel        <- sort(rel, decreasing = TRUE)
  i_noise    <- which(names(rel)=="noise")
  rel_rem    <- rel[-i_noise]
  rel        <- c(rel[i_noise], rel_rem)
  names(rel) <- c("noise", names(rel_rem))
  for(j in 1:length(rel)){
    h <- sum(rel[1:j])
    if(h > threshold){
      selected <- names(rel[1:j])
      res <- list(selected = selected, ev_sum = h, 
                  prop_ev = rel, threshold = threshold)
      return(res)
    }
  }
  res <- list(selected = names(rel), ev_sum = 1,
              prop_ev = rel, threshold = threshold)
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
                           p = function(x){stats::dbeta(x, 30, 3)},
                             h = 0.01)
{
  if(class(object)!="lgpfit") stop("Class of 'object' must be 'lgpfit'!")
  if(!is.function(p)){stop('p must be a function')}
  info <- object@model@info
  H    <- seq(0, 1, by = h)
  L    <- length(H)
  P    <- rep(0, L)
  prob <- selection_prob_fixed_threshold(object, threshold = 0)
  nam  <- names(prob)
  D    <- length(nam)
  PROB <- matrix(0, L, D)
  for(i in 1:L){
    P[i] <- p(H[i])
    PROB[i,] <- selection_prob_fixed_threshold(object, threshold = H[i])
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
#' @param object An object of class \code{lgpfit}.
#' @param threshold value between 0 and 1
#' @return selection probabilities for each component
selection_prob_fixed_threshold <- function(object, threshold = 0.95)
{
  if(class(object)!="lgpfit") stop("Class of 'object' must be 'lgpfit'!")
  if(threshold > 1 || threshold < 0) stop("'threshold' must be between 0 and 1!")
  info  <- object@model@info
  names <- info$component_names
  REL   <- object@relevances$samples
  n_cmp <- length(names)
  n_smp <- dim(REL)[1]
  SEL   <- matrix(0, n_smp, n_cmp)
  for(i in 1:n_smp){
    rel <- REL[i,]
    selected <- selection_helper(rel, threshold)
    for(j in 1:n_cmp){
      if(names[j] %in% selected){
        SEL[i,j] <- 1
      }
    }
  }
  colnames(SEL) <- names
  return(colMeans(SEL))
}


#' Select relevant components (one MCMC sample)
#'
#' @param rel a named vector of component relevances
#' @param threshold value between 0 and 1
#' @return names of selected components (including "noise" always)
selection_helper <- function(rel, threshold)
{
  rel        <- sort(rel, decreasing = TRUE)
  i_noise    <- which(names(rel)=="noise")
  rel_rem    <- rel[-i_noise]
  rel        <- c(rel[i_noise], rel_rem)
  names(rel) <- c("noise", names(rel_rem))
  for(j in 1:length(rel)){
    h <- sum(rel[1:j])
    if(h > threshold){
      SEL <- names(rel)[1:j]
      return(SEL)
    }
  }
  SEL <- names(rel)
  return(SEL)
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