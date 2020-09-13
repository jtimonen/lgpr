#' Character representations of different S4 objects
#'
#' @param x an object of some S4 class
#' @return a character representation of the object
#' @name as_character
NULL

#' @rdname as_character
setMethod("as.character", "lgpexpr", function(x) {
  paste0(x@fun, "(", x@covariate, ")")
})

#' @rdname as_character
setMethod("as.character", "lgpterm", function(x) {
  facs <- x@factors
  desc <- sapply(facs, as.character)
  desc <- paste(desc, collapse = "*")
  return(desc)
})

#' @rdname as_character
setMethod("as.character", "lgpformula", function(x) {
  return(x@call)
})


#' Operations on formula terms and expressions
#'
#' @description
#' \itemize{
#'   \item Sum two \linkS4class{lgprhs}'s to yield an \linkS4class{lgprhs}
#'   \item Sum two \linkS4class{lgpterm}'s to yield an \linkS4class{lgprhs}
#'   \item Sum an \linkS4class{lgprhs} and an \linkS4class{lgpterm}
#'   to yield an \linkS4class{lgprhs}
#'   \item Multiply two \linkS4class{lgpterm}'s to yield
#'   an \linkS4class{lgpterm}
#' }
#' @param e1 The first sum, term or expression
#' @param e2 The second sum, term or expression
#' @name operations
NULL

#' @rdname operations
setMethod(
  "+", signature(e1 = "lgprhs", e2 = "lgprhs"),
  function(e1, e2) {
    new("lgprhs", summands = c(e1@summands, e2@summands))
  }
)

#' @rdname operations
setMethod(
  "+", signature(e1 = "lgpterm", e2 = "lgpterm"),
  function(e1, e2) {
    new("lgprhs", summands = list(e1, e2))
  }
)

#' @rdname operations
setMethod(
  "+", signature(e1 = "lgprhs", e2 = "lgpterm"),
  function(e1, e2) {
    e1 + new("lgprhs", summands = list(e2))
  }
)

#' @rdname operations
setMethod(
  "*", signature(e1 = "lgpterm", e2 = "lgpterm"),
  function(e1, e2) {
    new("lgpterm", factors = c(e1@factors, e2@factors))
  }
)

#' Printing S4 object info using the show generic
#'
#' @name show
#' @param object an object of some S4 class
#' @return the object invisibly
NULL

#' @rdname show
setMethod("show", "lgpformula", function(object) {
  cat(as.character(object))
  invisible(object)
})

#' @rdname show
setMethod("show", "lgpsim", function(object) {
  msg <- class_info("lgpsim")
  cat(msg)
  invisible(object)
})

#' @rdname show
setMethod("show", "lgpmodel", function(object) {
  msg <- class_info("lgpmodel")
  cat(msg)
  cat("\n")
  model_summary(object)
})

#' @rdname show
setMethod("show", "lgpfit", function(object) {
  msg <- class_info("lgpfit")
  cat(msg)
  cat("\n")
  fit_summary(object)
})

#' Class info for show methods
#'
#' @param class_name class name
#' @return a string
class_info <- function(class_name) {
  str <- paste0(
    "An object of class ", class_name, ". See ?",
    class_name, " for more info."
  )
  return(str)
}


#' Visualize an artificial longitudinal data set
#'
#' @export
#' @description Creates plots where each observation unit has a separate panel.
#' @param simdata an object of class \linkS4class{lgpsim}
#' @param signal_name signal label to show in legend
#' @param signal_column name of the signal column in \code{simdata$components}
#' @param x_name name of x-axis variable
#' @param y_name name of y-axis variable
#' @param group_by grouping variable
#' @param ... keyword arguments to \code{\link{plot_panel}}
#' @return a \code{ggplot object}
plot_sim <- function(simdata, signal_name = "signal", signal_column = "f",
                     x_name = "age", y_name = "y", group_by = "id", ...) {
  data <- simdata@data
  df_points <- data[c(group_by, x_name, y_name)]
  df_lines <- df_points
  df_lines[[y_name]] <- simdata@components[[signal_column]]
  colnames(df_lines)[3] <- signal_name
  
  teff_true <- null_if_all_nan(simdata@effect_times$true)
  teff_obs <- null_if_all_nan(simdata@effect_times$observed)
  
  h <- plot_panel(
    df_data = df_points,
    df_signal = df_lines,
    teff_true = teff_true,
    teff_obs = teff_obs,
    ...
  )
  
  info <- paste0(
    "Vertical lines are the real effect time (solid) \n",
    "and observed disease onset / initiation time (dashed)."
  )
  h <- h + ggplot2::ggtitle("Simulated data", subtitle = info)
  return(h)
}
