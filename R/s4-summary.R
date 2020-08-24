
#' @rdname prior_summary
setMethod(
  "prior_summary", signature(object = "lgpmodel"),
  function(object, digits = 3) {
    prior_to_df(object@stan_input, digits = digits)
  }
)

#' @rdname prior_summary
setMethod(
  "prior_summary", signature(object = "lgpfit"),
  function(object, digits = 3) {
    prior_to_df(object@model@stan_input, digits = digits)
  }
)
