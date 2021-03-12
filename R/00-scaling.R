
#' Create a standardizing transform
#'
#' @param x variable measurements (might contain \code{NA} or \code{NaN})
#' @param name variable name
#' @return an object of class \linkS4class{lgpscaling}
#' @family variable scaling functions
create_scaling <- function(x, name) {
  check_length_geq(x, 2)
  loc <- mean(x, na.rm = TRUE)
  scale <- stats::sd(x, na.rm = TRUE)
  if (scale == 0) {
    msg <- paste0("the variable <", name, "> has zero variance!")
    stop(msg)
  }
  new("lgpscaling", loc = loc, scale = scale, var_name = name)
}


#' Apply variable scaling
#'
#' @param scaling an object of class \linkS4class{lgpscaling}
#' @param x object to which apply the scaling (numeric)
#' @param inverse whether scaling should be done in inverse direction
#' @return a similar object as \code{x}
#' @family variable scaling functions
apply_scaling <- function(scaling, x, inverse = FALSE) {
  loc <- scaling@loc
  scale <- scaling@scale
  if (inverse) {
    x <- scale * x + loc
  } else {
    x <- (x - loc) / scale
  }
  return(x)
}
