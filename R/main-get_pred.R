#' Extract model predictions and function posteriors
#'
#' @inheritParams pred
#' @return an object of class \linkS4class{GaussianPrediction} or
#' \linkS4class{Prediction}
get_pred <- function(fit, draws = NULL, reduce = NULL, verbose = TRUE) {
  check_type(fit, "lgpfit")
  if (!contains_postproc(fit)) {
    warning("fit contains no postprocessing information, need to call postproc")
    fit <- postproc(fit)
  }
  pr <- fit@postproc_results
  return(dollar(pr, "pred"))
  return(pred)
}
