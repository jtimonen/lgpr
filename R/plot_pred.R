#' Visualizing model predictions or inferred covariate effects
#'
#' @description
#' \itemize{
#'   \item Predictions at data points can be visualized using \code{plot_pred}.
#'     Out-of-sample predictions can be visualized by giving the \code{pred}
#'     and \code{x} arguments.
#'   \item The total signal \code{f} or any of its
#'   additive components can be plotted using \code{plot_f}.
#' }
#' @param fit An object of class \linkS4class{lgpfit}.
#' @param pred An object of class \linkS4class{GaussianPrediction} or
#' \linkS4class{Prediction}.
#' @param x A data frame of prediction points. Must be specified if
#' \code{pred} is not \code{NULL}.
#' @param t_name name of the x-axis variable
#' @param group_by name of the grouping variable
#' @param draws only has effect if \code{pred} is \code{NULL}
#' @param reduce only has effect if \code{pred} is \code{NULL}
#' @param MULT_STD a multiplier for standard deviation
#' @param ... additional arguments to \code{\link{plot_api_g}} or
#' \code{\link{plot_api_c}}
#' @return a \code{\link[ggplot2]{ggplot}} object
#' @family main plot functions
#' @name plot_pred
NULL

#' @export
#' @rdname plot_pred
plot_pred <- function(fit,
                      pred = NULL,
                      x = NULL,
                      group_by = "id",
                      t_name = "age",
                      draws = NULL,
                      reduce = NULL,
                      MULT_STD = 2.0,
                      ...) {
  check_type(fit, "lgpfit")
  df_data <- create_plot_df(fit, t_name, group_by)
  if (is.null(pred)) {
    pred <- get_pred(fit, draws, reduce)
    df_base <- create_plot_df(fit, t_name, group_by)[, 1:2]
  } else {
    df_base <- data.frame(dollar(x, group_by), dollar(x, t_name))
    colnames(df_base) <- c(group_by, t_name)
  }
  PRED <- plot_pred.create(fit, pred, df_base, MULT_STD)
  h <- plot_api_g(
    df_data = df_data,
    df_fit =  dollar(PRED, "df_line"),
    df_fit_err = dollar(PRED, "df_ribbon"),
    ...
  )
  return(h)
}

#' @export
#' @rdname plot_pred
plot_f <- function(fit,
                   component_idx,
                   color_by = NA,
                   x_name = "age",
                   group_by = "id",
                   draws = NULL,
                   ...) {
  check_type(fit, "lgpfit")
  h <- 0
  return(h)
}


#' Helper function for plot_pred
#'
#' @inheritParams plot_pred
#' @param df_base a data frame with two columns
#' @return a list
plot_pred.create <- function(fit, pred, df_base, MULT_STD) {

  # Create fit line(s)
  df_line <- plot_pred.create.df_line(fit, pred, df_base)
  df_ribbon <- plot_pred.create.df_ribbon(fit, pred, df_base, MULT_STD)
  info <- NULL

  # Return
  list(
    df_line = df_line,
    df_ribbon = df_ribbon,
    info = info
  )
}

#' @rdname plot_pred.create
plot_pred.create.df_line <- function(fit, pred, df_base) {
  line_arr <- if (is_f_sampled(fit)) pred@h else pred@y_mean
  num_draws <- nrow(line_arr)
  df_wide <- make_draw_df(line_arr)
  df_long <- to_long_format(df_wide)
  colnames(df_long) <- c("_draw_", "y")
  df_base <- rep_df(df_base, times = num_draws)
  out <- cbind(df_base, df_long)
  if (num_draws == 1) out[["_draw_"]] <- NULL
  return(out)
}

#' @rdname plot_pred.create
plot_pred.create.df_ribbon <- function(fit, pred, df_base, MULT_STD) {
  check_positive(MULT_STD)
  line_arr <- if (is_f_sampled(fit)) pred@h else pred@y_mean
  num_draws <- nrow(line_arr)
  if (num_draws > 1 || is_f_sampled(fit)) {
    return(NULL)
  }
  m <- as.vector(line_arr)
  std <- as.vector(pred@y_std)
  upper <- m + MULT_STD * std
  lower <- m - MULT_STD * std
  df <- data.frame(upper, lower)
  cbind(df_base, df)
}
