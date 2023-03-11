#' Visualizing model predictions or inferred covariate effects
#'
#' @description
#' \itemize{
#'   \item Function draws at data points can be visualized using
#'   \code{plot_pred}. If the \code{pred} argument is \code{NULL}, it
#'   is computed using the \code{\link{pred}} function with \code{x=NULL}.
#'   \item The total signal \code{f} or any of its
#'   additive components can be plotted using \code{plot_f}.
#' }
#' @param fit An object of class \linkS4class{lgpfit}.
#' @param x Deprecated argument. This is now taken from the \code{pred}
#' object to ensure compatibility.
#' @param pred An object of class \linkS4class{GaussianPrediction} or
#' \linkS4class{Prediction}. If \code{pred=NULL}, the \code{\link{pred}}
#' function is called with the given \code{reduce} and \code{draws} arguments.
#' @param t_name name of the x-axis variable
#' @param group_by name of the grouping variable (use \code{group_by=NA}
#' to avoid grouping)
#' @param draws Only has effect if \code{pred=NULL}.
#' @param reduce Only has effect if \code{pred=NULL}.
#' @param MULT_STD a multiplier for standard deviation
#' @param ... additional arguments to \code{\link{plot_api_g}} or
#' \code{\link{plot_api_c}}
#' @param verbose Can this print any messages?
#' @return a \code{\link[ggplot2]{ggplot}} object
#' @family main plot functions
#' @name plot_pred
NULL

#' @export
#' @rdname plot_pred
plot_pred <- function(fit,
                      pred = NULL,
                      group_by = "id",
                      t_name = "age",
                      MULT_STD = 2.0,
                      verbose = TRUE,
                      draws = NULL,
                      reduce = function(x) base::mean(x),
                      # deprecated
                      x = NULL,
                      ...) {
  # Process arguments, pred can't be NULL after this
  df_data <- create_plot_df(fit, t_name, group_by)
  pred <- plot_pred.process_args(fit, pred, x, draws, reduce, verbose)
  x <- pred@x

  # Create input to plotting function
  input <- plot_pred.create_input(pred, x, draws, reduce, group_by, t_name)
  pred <- dollar(input, "pred")
  df_base <- dollar(input, "df_base")
  y_name <- get_y_name(fit)

  # Create y transform
  lh <- get_obs_model(fit)
  if (is_bin_or_bb(lh)) {
    y_fun <- function(x) divide_by_num_trials(x, fit)
    y_lab <- paste0(y_name, " / num_trials")
  } else {
    y_fun <- function(x) x
    y_lab <- y_name
  }

  # Create plot
  h <- plot_api_g(
    df_data = df_data,
    df = plot_pred.create.df_line(fit, pred, df_base),
    df_err = plot_pred.create.df_ribbon(fit, pred, df_base, MULT_STD),
    y_transform = y_fun,
    ...
  )
  h + ggplot2::ylab(y_lab)
}


#' @export
#' @param comp_idx Index of component to plot. The total sum is plotted
#' if this is \code{NULL}.
#' @param color_by name of coloring factor
#' @rdname plot_pred
plot_f <- function(fit,
                   pred = NULL,
                   group_by = "id",
                   t_name = "age",
                   MULT_STD = 2.0,
                   verbose = TRUE,
                   draws = NULL,
                   reduce = function(x) base::mean(x),
                   comp_idx = NULL,
                   color_by = NA,
                   # deprecated
                   x = NULL,
                   ...) {
  # Process arguments, pred can't be NULL after this
  pred <- plot_pred.process_args(fit, pred, x, draws, reduce, verbose)
  x <- pred@x
  color_fac <- plot_f.get_color_fac(fit, pred, x, color_by)

  # Create input to plotting function
  input <- plot_pred.create_input(pred, x, draws, reduce, group_by, t_name)
  pred <- dollar(input, "pred")
  df_base <- dollar(input, "df_base")

  # Add column to df_base
  bn <- colnames(df_base)
  df_base <- cbind(df_base, color_fac)
  color_by_name <- if (is.na(color_by)) NA else paste0(color_by, "__")
  colnames(df_base) <- c(bn, color_by_name)
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
                            pred = NULL,
                            group_by = "id",
                            t_name = "age",
                            MULT_STD = 2.0,
                            verbose = TRUE,
                            draws = NULL,
                            reduce = function(x) base::mean(x),
                            color_by = NA,
                            no_err = FALSE,
                            ylim = NULL,
                            draw = TRUE,
                            nrow = NULL,
                            ncol = NULL,
                            gg_add = NULL,
                            # deprecated
                            x = NULL,
                            ...) {
  # Process args, pred can't be NULL after this
  pred <- plot_pred.process_args(fit, pred, x, draws, reduce, verbose)

  # Setup
  cn <- component_names(fit@model)
  J <- length(cn)
  out <- list()
  check_length_1_or(color_by, J + 1)
  check_length_1_or(no_err, J + 1)

  # Loop thorough components
  for (j in seq_len(J + 1)) {
    cb <- if (length(color_by) == 1) color_by else color_by[j]
    ne <- if (length(no_err) == 1) no_err else no_err[j]
    comp_idx <- if (j > J) NULL else j
    title <- if (j > J) "sum" else cn[j]
    p <- plot_f(
      fit = fit,
      pred = pred,
      group_by = group_by,
      t_name = t_name,
      draws = draws,
      reduce = reduce,
      MULT_STD = MULT_STD,
      comp_idx = comp_idx,
      color_by = cb,
      no_err = ne,
      verbose = verbose,
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


# Process arguments to plot_pred or plot_f
plot_pred.process_args <- function(fit, pred_in, x, draws, reduce, verbose) {
  check_type(fit, "lgpfit")
  vrb <- verbose
  if (!is.null(x)) {
    warning("The <x> argument is deprecated and has no effect..")
  }
  if (is.null(pred_in)) {
    log_info(
      "Missing 'pred' argument, computing or extracting it...",
      verbose
    )
    pred_out <- pred(fit, draws = draws, reduce = reduce, verbose = vrb)
  } else {
    pred_out <- pred_in
  }
  return(pred_out)
}

# Create a list with names \code{pred} and \code{df_base}
plot_pred.create_input <- function(pred, x, draws, reduce, group_by, t_name) {
  x_grp <- create_grouping_factor(x, group_by) # util
  df_base <- data.frame(x_grp, dollar(x, t_name))
  group_by <- if (is.na(group_by)) "group__" else group_by
  colnames(df_base) <- c(group_by, t_name)
  list(pred = pred, df_base = df_base)
}

# Helper functions for plot_pred
#
# @param df_base a data frame with two columns
# @return a data frame
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


# Helper functions for plot_f
#
# @param df_base a data frame with two columns
# @return a data frame
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

plot_f.create.line_arr <- function(fit, pred, comp_idx) {
  if (is.null(comp_idx)) {
    line_arr <- if (is_f_sampled(fit)) pred@f else pred@f_mean
  } else {
    f_cmp <- if (is_f_sampled(fit)) pred@f_comp else pred@f_comp_mean
    line_arr <- f_cmp[[comp_idx]]
  }
  return(line_arr)
}

# Create a factor for coloring
plot_f.get_color_fac <- function(fit, pred, x, color_by) {
  n <- if (is.null(pred)) get_num_obs(fit) else nrow(x)
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
