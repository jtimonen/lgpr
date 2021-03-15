setGeneric(
  "parameter_info",
  function(object, digits) standardGeneric("parameter_info")
)

setGeneric(
  "component_info", function(object) standardGeneric("component_info")
)

setGeneric(
  "covariate_info", function(object) standardGeneric("covariate_info")
)

setGeneric(
  "component_names", function(object) standardGeneric("component_names")
)

setGeneric(
  "get_model",
  function(object) standardGeneric("get_model")
)

setGeneric(
  "get_stanfit",
  function(object) standardGeneric("get_stanfit")
)

setGeneric(
  "get_draws",
  function(object, draws = NULL, reduce = NULL, ...) standardGeneric("get_draws")
)
