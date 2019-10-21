#' Warp input (Cpp)
#' 
#' @param t an input
#' @param a an input
#' @param b an input
#' @param c an input
#' @export
stanfunc_warp_input <- function(t, a, b, c) {
  return(STANFUNC_warp_input(t, a, b, c))
}