#' Check that data contains a variable with a certain name
#'
#' @param var_name the variable to be searched for
#' @param data an object of class \code{data.frame}
#' @return \code{TRUE} if the variable is found
check_in_data <- function(var_name, data) {
  d_names <- colnames(data)
  ok <- (var_name %in% d_names)
  if (!ok) {
    str <- paste(d_names, collapse = ", ")
    msg <- paste0(
      "The variable '", var_name, "' not found in <data>! ",
      " Found data columns = {", str, "}."
    )
    stop(msg)
  }
  return(TRUE)
}

#' Create the function that does a standardizing transform and its inverse
#'
#' @param y variable measurements
#' @param name variable name
#' @return an object of class \linkS4class{lgpscaling}
create_scaling <- function(y, name) {
  if (length(y) < 2) {
    stop("length of <y> must be at least 2!")
  }
  m <- mean(y)
  std <- stats::sd(y)
  if (std == 0) {
    stop("the varible measurements have zero variance!")
  }
  fun <- function(x) {
    (x - m) / std
  }
  fun_inv <- function(x) {
    x * std + m
  }
  new("lgpscaling", fun = fun, fun_inv = fun_inv, var_name = name)
}
