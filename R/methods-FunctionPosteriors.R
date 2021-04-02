#' @describeIn FunctionPosteriors Print a summary about the object.
setMethod("show", "FunctionPosteriors", function(object) {
  comp_names <- component_names(object@model)
  n_eval_points <- nrow(object@x)
  DIMS <- c(object@num_paramsets, n_eval_points)
  desc <- class_info_fp("FunctionPosteriors", comp_names, DIMS)
  cat(desc)
})


#' @describeIn FunctionDraws Print a summary about the object.
setMethod("show", "FunctionDraws", function(object) {
  df <- object@components
  comps <- levels(dollar(df, "component"))
  paramsets <- levels(dollar(df, "paramset"))
  eval_points <- levels(dollar(df, "eval_point"))
  num_comps <- length(comps)
  DIMS <- c(length(paramsets), length(eval_points))
  desc <- class_info_fp("FunctionDraws", num_comps, DIMS)
  cat(desc)
})

# Class info for show methods of function posterior objects
class_info_fp <- function(class_name, comp_names, D) {
  num_comps <- length(comp_names)
  comp_str <- paste(comp_names, collapse = ", ")
  desc <- class_info(class_name)
  desc <- paste0(desc, "\n - ", num_comps, " components: ", comp_str)
  desc <- paste0(desc, "\n - ", D[1], " parameter set(s)")
  desc <- paste0(desc, "\n - ", D[2], " evaluation points")
  return(desc)
}

#' @describeIn FunctionPosteriors Access the main data frames.
#' @param paramset index of parameter set to pick (NULL = all)
#' @param component name of component to pick (NULL = all)
#' @param covariates covariates to include (NULL = none)
setMethod(
  "get_df", "FunctionPosteriors",
  function(object, paramset = NULL, component = NULL,
           covariates = NULL) {
    df <- object@f
    S <- object@num_paramsets
    J <- get_num_comps(object@model)
    num_repeat_x <- S * (J + 1)
    if (!is.null(paramset)) {
      df <- df[dollar(df, "paramset") == paramset, , drop = FALSE]
      rownames(df) <- NULL
      num_repeat_x <- num_repeat_x / S
    }
    if (!is.null(component)) {
      df <- df[dollar(df, "component") == component, , drop = FALSE]
      rownames(df) <- NULL
      num_repeat_x <- num_repeat_x / (J + 1)
    }
    if (nrow(df) == 0) stop("zero rows selected!")
    if (!is.null(covariates)) {
      x_add <- object@x[, covariates, drop = FALSE]
      x_add <- rep_df(x_add, times = num_repeat_x)
      if (nrow(x_add) != nrow(df)) stop("differing number of rows!")
      df <- cbind(x_add, df)
    }
    return(df)
  }
)

#' @describeIn FunctionPosteriors Visualization of function posteriors.
#' Optional arguments (\code{...}) are passed to
#' \code{\link{plot_FunctionPosteriors}}.
#' @param x a \linkS4class{FunctionPosteriors} object to visualize
#' @param y unused argument
setMethod(
  "plot",
  signature = c("FunctionPosteriors", "missing"),
  function(x, y, ...) {
    plot_FunctionPosteriors(x, ...)
  }
)

#' Visualize function posterior distributions
#'
#' @param fp An object of class \linkS4class{FunctionPosteriors}.
#' @param paramset index of parameter set to pick (NULL = all, but only
#' posterior means will be plotted)
#' @param component name of component to pick (NULL = all)
#' @param t_name name of the x-axis variable
#' @param group_by name of the grouping factor (use \code{group_by=NULL}
#' to avoid grouping)
#' @param color_by name of coloring factor (use \code{color_by=NULL}
#' to avoid coloring)
#' @param MULT_STD a multiplier for standard deviation
#' @param alpha line opacity
#' @param alpha_err ribbon opacity
#' @param no_err hide error bar even when it would normally be plotted?
#' @param no_line hide line even when it would normally be plotted?
#' @return a \code{\link[ggplot2]{ggplot}} object
#' @param verbose Can this print any messages?
plot_FunctionPosteriors <- function(fp,
                                    paramset = NULL,
                                    component = NULL,
                                    group_by = "id",
                                    color_by = NULL,
                                    t_name = "age",
                                    MULT_STD = 2.0,
                                    alpha_err = 0.2,
                                    alpha = 1.0,
                                    no_err = FALSE,
                                    no_line = FALSE,
                                    verbose = TRUE) {

  # Get data frame for ggplot
  covariates <- c(t_name, group_by, color_by)
  df <- get_df(fp, paramset, component, covariates)

  # Add possible factor crossing with paramset
  fac1 <- dollar(df, "paramset")
  if (is.null(group_by)) {
    df[["GGPLOT_GROUP"]] <- fac1
  } else {
    fac2 <- dollar(df, group_by)
    df[["GGPLOT_GROUP"]] <- as.factor(paste(fac1, fac2, sep = "x"))
  }

  # Create aesthetics and add possible coloring factor
  aes <- ggplot2::aes_string(
    x = t_name, y = "mean", group = "GGPLOT_GROUP", color = color_by
  )

  # Create ggplot object
  plt <- ggplot2::ggplot(df, aes) +
    ggplot2::facet_grid(. ~ component)
  do_ribbon <- FALSE
  if (!no_err) {
    S <- length(levels(dollar(df, "paramset")))
    if (S == 1) {
      do_ribbon <- TRUE
      if (!is.null(color_by)) {
        fac3 <- dollar(df, color_by)
        if (class(fac3) != "factor") {
          color_by_rib <- NULL
          msg <- paste0(
            "Coloring variable (color_by=", color_by, ") is not a factor, ",
            "so not coloring the ribbon."
          )
          if (verbose) cat(msg, "\n")
        }
        else {
          color_by_rib <- color_by
        }
      } else {
        color_by_rib <- color_by
      }
      aes_rib <- ggplot2::aes_string(
        x = t_name,
        ymin = paste0("mean - ", MULT_STD, " * sd"),
        ymax = paste0("mean + ", MULT_STD, " * sd"),
        group = group_by,
        color = color_by_rib,
        fill = color_by_rib
      )
      plt <- plt + ggplot2::geom_ribbon(
        mapping = aes_rib,
        alpha = alpha_err,
        color = NA # no ribbon edge
      )
    } else {
      msg <- paste0(
        "Number of parameter sets > 1, so plotting only the ",
        "conditional mean given each parameter set"
      )
      if (verbose) cat(msg, "\n")
    }
  }
  if (!no_line) {
    plt <- plt + ggplot2::geom_line(alpha = alpha)
  }
  ylab <- if (do_ribbon) paste0("mean +- ", MULT_STD, " x std") else "mean"
  plt <- plt + ggplot2::ylab(ylab)
  return(plt)
}
