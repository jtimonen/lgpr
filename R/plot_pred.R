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
#' @param x A data frame of prediction points, containing at least the
#' variables used as covariates in the model. This must be specified if
#' \code{pred} is not \code{NULL}.
#' @param pred An object of class \linkS4class{GaussianPrediction} or
#' \linkS4class{Prediction}.
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
                      x = NULL,
                      pred = NULL,
                      group_by = "id",
                      t_name = "age",
                      draws = NULL,
                      reduce = mean,
                      MULT_STD = 2.0,
                      ...) {
  check_type(fit, "lgpfit")
  df_data <- create_plot_df(fit, t_name, group_by)
  input <- plot_pred.create_input(fit, pred, x, draws, reduce, group_by, t_name)
  pred <- dollar(input, "pred")
  df_base <- dollar(input, "df_base")

  # Create plot
  plot_api_g(
    df_data = df_data,
    df = plot_pred.create.df_line(fit, pred, df_base),
    df_err = plot_pred.create.df_ribbon(fit, pred, df_base, MULT_STD),
    ...
  )
}


#' @export
#' @param comp_idx Index of component to plot. The total sum is plotted
#' if this is \code{NULL}.
#' @param color_by name of coloring factor
#' @rdname plot_pred
plot_f <- function(fit,
                   x = NULL,
                   pred = NULL,
                   group_by = "id",
                   t_name = "age",
                   draws = NULL,
                   reduce = mean,
                   MULT_STD = 2.0,
                   comp_idx = NULL,
                   color_by = NA,
                   ...) {
  check_type(fit, "lgpfit")
  color_fac <- plot_f.get_color_fac(fit, pred, x, color_by)
  input <- plot_pred.create_input(fit, pred, x, draws, reduce, group_by, t_name)
  pred <- dollar(input, "pred")
  df_base <- dollar(input, "df_base")
  bn <- colnames(df_base)
  df_base <- cbind(df_base, color_fac)
  colnames(df_base) <- c(bn, color_by)
  ylab <- if (is.null(comp_idx)) "f" else paste0("f[", comp_idx, "]")

  # Create plot
  h <- plot_api_c(
    df = plot_f.create.df_line(fit, pred, df_base, comp_idx),
    df_err = plot_f.create.df_ribbon(fit, pred, df_base, comp_idx, MULT_STD),
    ...
  )
  h + ggplot2::ylab(ylab)
}

#' Visualize all model components
#'
#' @export
#' @description
#' This calls \code{\link{plot_f}} for all model components.
#' @inheritParams plot_f
#' @param color_by Names of coloring factors. Can have length 1 or equal to
#' the number of components. See the \code{color_by} argument of
#' \code{\link{plot_f}}.
#' @param no_err Should the error ribbons be skipped even though they
#' otherwise would be shown? Can have length 1 or equal to number of
#' components + 1. See the \code{no_err} argument of \code{\link{plot_api_c}}.
#' @param ylim a vector of length 2 (upper and lower y-axis limits), or NULL
#' @param ... additional arguments to \code{\link{plot_api_c}}
#' @param draw if this is TRUE, the plot grid is drawn using
#' \code{\link[gridExtra]{arrangeGrob}}
#' @param nrow number of grid rows
#' @param ncol number of grid columns
#' @param gg_add additional ggplot obejct to add to each plot
#' @return a list of ggplot objects invisibly
#' @family main plot functions
plot_components <- function(fit,
                            x = NULL,
                            pred = NULL,
                            group_by = "id",
                            t_name = "age",
                            draws = NULL,
                            reduce = mean,
                            MULT_STD = 2.0,
                            color_by = NA,
                            no_err = FALSE,
                            ylim = NULL,
                            draw = TRUE,
                            nrow = NULL,
                            ncol = NULL,
                            gg_add = NULL,
                            ...) {
  check_type(fit, "lgpfit")
  cn <- component_names(fit)
  J <- length(cn)
  out <- list()
  check_length_1_or(color_by, J + 1)
  check_length_1_or(no_err, J + 1)
  for (j in seq_len(J + 1)) {
    cb <- if (length(color_by) == 1) color_by else color_by[j]
    ne <- if (length(no_err) == 1) no_err else no_err[j]
    comp_idx <- if (j > J) NULL else j
    title <- if (j > J) "sum" else cn[j]
    p <- plot_f(
      fit = fit,
      x = x,
      pred = pred,
      group_by = group_by,
      t_name = t_name,
      draws = draws,
      reduce = reduce,
      MULT_STD = MULT_STD,
      comp_idx = comp_idx,
      color_by = cb,
      no_err = ne,
      ...
    )
    p <- if (is.null(ylim)) p else p + ggplot2::ylim(ylim[1], ylim[2])
    p <- p + ggplot2::ggtitle(title) + gg_add
    out[[j]] <- p
  }

  # Draw
  if (draw) gridExtra::grid.arrange(grobs = out, nrow = nrow, ncol = ncol)
  invisible(out)
}

#' Helper function for plot_pred and plot_f
#'
#' @inheritParams plot_pred
#' @return a list with names \code{pred} and \code{df_base}
plot_pred.create_input <- function(fit, pred, x, draws, reduce,
                                   group_by, t_name) {
  if (is.null(pred)) {
    pred <- get_pred(fit, draws, reduce)
    df_base <- create_plot_df(fit, t_name, group_by)[, 1:2]
  } else {
    df_base <- data.frame(dollar(x, group_by), dollar(x, t_name))
    colnames(df_base) <- c(group_by, t_name)
  }
  list(pred = pred, df_base = df_base)
}

#' Helper functions for plot_pred
#'
#' @inheritParams plot_pred
#' @param df_base a data frame with two columns
#' @return a data frame
#' @name plot_pred.create
NULL

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

#' Helper function for plot_f
#'
#' @inheritParams plot_f
#' @return a factor
plot_f.get_color_fac <- function(fit, pred, x, color_by) {
  n <- if (is.null(pred)) get_stan_input(fit)$num_obs else nrow(x)
  if (is.na(color_by)) {
    color_fac <- as.factor(rep(1, n))
    return(color_fac)
  }
  if (is.null(pred)) {
    dat <- get_data(fit)
    color_fac <- dollar(dat, color_by)
  } else {
    color_fac <- dollar(x, color_by)
  }
  if (!is.factor(color_fac)) {
    # Color by whether or not the coloring variable is NA or NaN
    color_fac <- as.numeric(!is.na(color_fac))
    color_fac <- as.factor(c("N/A", "available")[color_fac + 1])
  }
  return(color_fac)
}


#' Helper functions for plot_f
#'
#' @inheritParams plot_f
#' @param df_base a data frame with two columns
#' @return a data frame
#' @name plot_f.create
NULL

#' @rdname plot_f.create
plot_f.create.df_line <- function(fit, pred, df_base, comp_idx) {
  line_arr <- plot_f.create.line_arr(fit, pred, comp_idx)
  num_draws <- nrow(line_arr)
  df_wide <- make_draw_df(line_arr)
  df_long <- to_long_format(df_wide)
  colnames(df_long) <- c("_draw_", "y")
  df_base <- rep_df(df_base, times = num_draws)
  out <- cbind(df_base, df_long)
  if (num_draws == 1) out[["_draw_"]] <- NULL
  return(out)
}

#' @rdname plot_f.create
plot_f.create.df_ribbon <- function(fit, pred, df_base, comp_idx,
                                    MULT_STD) {
  check_positive(MULT_STD)
  line_arr <- plot_f.create.line_arr(fit, pred, comp_idx)
  num_draws <- nrow(line_arr)
  if (num_draws > 1 || is_f_sampled(fit)) {
    return(NULL)
  }
  m <- as.vector(line_arr)
  f_std <- if (is.null(comp_idx)) pred@f_std else pred@f_comp_std[[comp_idx]]
  std <- as.vector(f_std)
  upper <- m + MULT_STD * std
  lower <- m - MULT_STD * std
  df <- data.frame(upper, lower)
  cbind(df_base, df)
}

#' @rdname plot_f.create
plot_f.create.line_arr <- function(fit, pred, comp_idx) {
  if (is.null(comp_idx)) {
    line_arr <- if (is_f_sampled(fit)) pred@f else pred@f_mean
  } else {
    f_cmp <- if (is_f_sampled(fit)) pred@f_comp else pred@f_comp_mean
    line_arr <- f_cmp[[comp_idx]]
  }
  return(line_arr)
}
