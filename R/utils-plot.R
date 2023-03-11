#' Plot longitudinal data and/or model fit so that each subject/group has
#' their own panel
#'
#' @description Data frames specified in arguments \code{df_data},
#' \code{df_signal}, \code{df_fit}, and \code{df_fit_err} must have a format
#' where
#' \itemize{
#'   \item the first column is the grouping factor (usually id)
#'   \item the second column is the x-axis variable (usually age)
#'   \item a column named \code{y} must contain the y-axis variable
#'   (not for \code{df_fit_err})
#'   \item a column named \code{lower} (\code{upper}) must contain the lower
#'   (upper) bound of error bar (only for \code{df_fit_err})
#'   \item a column named \code{draw} must be a factor that
#'   specifies the posterior draw using which the fit has been computed
#'   (only for \code{df_fit})
#' }
#' @param df_data A data frame containing the observations.
#' @param df_signal A data frame containing the true signal. Omitted if
#' \code{NULL}.
#' @param df A data frame containing the model fit, or a list of data
#' frames. The list version can be used for example so that each list element
#' corresponds to the fit computed using one parameter draw. Omitted if
#' \code{NULL}.
#' @param df_err A data frame containing error bars. Omitted if \code{NULL}.
#' Must be \code{NULL} if \code{df_fit} is a list.
#' @param teff_signal A named vector containing true effect times used to
#' generate the signal. Omitted if \code{NULL}.
#' @param teff_obs A named vector containing observed effect times. Omitted if
#' \code{NULL}.
#' @param i_test Indices of test points.
#' @param color_signal Line color for true signal.
#' @param color Line color for model fit.
#' @param color_err Color of the error ribbon.
#' @param color_vlines Two line colors for vertical lines
#' (true and obs. effect time).
#' @param alpha Line opacity for model fit.
#' @param alpha_err Opacity of the error ribbon.
#' @param nrow number of rows, an argument for
#' \code{\link[ggplot2]{facet_wrap}}
#' @param ncol number of columns, an argument for
#' \code{\link[ggplot2]{facet_wrap}}
#' @param y_transform A function to be applied to the third column of
#' \code{df_data}.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @family internal plot API functions
plot_api_g <- function(df_data,
                       df_signal = NULL,
                       df = NULL,
                       df_err = NULL,
                       teff_signal = NULL,
                       teff_obs = NULL,
                       i_test = NULL,
                       color_signal = color_palette(2)[1],
                       color = color_palette(2)[2],
                       color_err = colorset("red", "light_highlight"),
                       color_vlines = colorset("gray", "mid_highlight"),
                       alpha = 1.0,
                       alpha_err = 0.5,
                       nrow = NULL,
                       ncol = NULL,
                       y_transform = function(x) x) {
  check_type(df_data, "data.frame")
  check_type(y_transform, "function")

  # Create the plot
  df_data[, 3] <- call_fun(y_transform, df_data[, 3])
  h <- plot_api_g_create(df_data, i_test, nrow, ncol)
  group_by <- colnames(df_data)[1]
  x_name <- colnames(df_data)[2]

  # Add true and observed effect times
  col <- color_vlines
  h <- plot_api_g_add_effect_times(h, teff_signal, group_by, col, 1, 0.7)
  h <- plot_api_g_add_effect_times(h, teff_obs, group_by, col, 2, 0.7)

  # Add other layers
  h <- add_g_layer_ribbon(
    h, df_err, x_name, group_by,
    color_err, alpha_err
  )
  h <- add_g_layer_fit(h, df, x_name, group_by, color, alpha)
  h <- add_g_layer_signal(h, df_signal, x_name, group_by, color_signal)

  # Add data and return
  h <- h + ggplot2::geom_point()
  return(h)
}

#' Plot a generated/fit model component
#'
#' @description Data frames specified in arguments \code{df},
#' and \code{df_err} must have a format where
#' \itemize{
#'   \item The first column is the grouping factor (usually id).
#'   \item The second column is the x-axis variable (usually age).
#'   \item The third column is the coloring factor. If name of the third
#'   column is \code{NA}, coloring is not done.
#'   \item A column named \code{y} must contain the y-axis variable
#'   (not for \code{df_err}).
#'   \item A column named \code{lower} (\code{upper}) must contain the lower
#'   (upper) bound of error bar (only for \code{df_err}).
#'   \item The posterior draw using which the fit has been computed can be
#'   specified with a factor named \code{_draw_} (only for \code{df}).
#' }
#' @param df a data frame
#' @param df_err a data frame
#' @param alpha line opacity
#' @param alpha_err ribbon opacity
#' @param no_err hide error bar even when it would normally be plotted?
#' @param no_line hide line even when it would normally be plotted?
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @family internal plot API functions
plot_api_c <- function(df,
                       df_err = NULL,
                       alpha = 1.0,
                       alpha_err = 0.2,
                       no_err = FALSE,
                       no_line = FALSE) {
  # Create the plot
  group_by <- colnames(df)[1]
  x_name <- colnames(df)[2]
  color_by <- colnames(df)[3]
  color_by <- if (is.na(color_by)) NULL else color_by
  aes <- ggplot2::aes_string(x = x_name, y = "y")
  h <- ggplot2::ggplot(data = df, mapping = aes)

  # Add layers
  if (!no_err) {
    h <- add_c_layer_ribbon(h, df_err, x_name, group_by, color_by, alpha_err)
  }
  if (!no_line) {
    h <- add_c_layer_line(h, df, x_name, group_by, color_by, alpha)
  }
  return(h)
}

# Initialize a grouped plot
plot_api_g_create <- function(df_data, i_test, nrow, ncol) {
  # Create dummy variable for type (train/test)
  show_test <- !is.null(i_test)
  if (show_test) {
    n <- dim(df_data)[1]
    a <- rep("train", n)
    a[i_test] <- "test"
    type_name <- "Type"
    df_data[[type_name]] <- as.factor(a)
  } else {
    type_name <- NULL
  }

  # Create aes
  aes <- plot_api_g_create_aes(df_data, type_name)
  h <- ggplot2::ggplot(df_data, aes)

  # Create faceting
  group_by <- colnames(df_data)[1]
  h <- add_g_layer_faceting(h, group_by, nrow, ncol)

  # Add shape and color scale
  if (show_test) {
    h <- h + ggplot2::scale_shape_manual(values = c(16, 16))
    col3 <- color_palette(3)[3]
    cols <- c(col3, "black")
    h <- h + ggplot2::scale_color_manual(values = cols)
  }
  return(h)
}

# Helper function for plot_grouped (creates a ggplot aes object)
#
# @param data a data frame with three columns
# @param type_name name of a factor that is to be distinguished by
# plot marker and/or type
plot_api_g_create_aes <- function(data, type_name) {
  group_by <- colnames(data)[1]
  x_name <- colnames(data)[2]
  y_name <- colnames(data)[3]
  if (!is.null(type_name)) {
    aes <- ggplot2::aes_string(
      x = x_name,
      y = y_name,
      group = group_by,
      pch = type_name,
      color = type_name
    )
  } else {
    aes <- ggplot2::aes_string(x = x_name, y = y_name, group = group_by)
  }
  return(aes)
}

# Add real or observed disease effect times to a grouped plot
# @param teff a named vector
plot_api_g_add_effect_times <- function(h, teff, group_by,
                                        linecolor, linetype, lwd) {
  if (!is.null(teff)) {
    nams <- names(teff)
    z <- as.numeric(teff)
    df <- data.frame(z, nams)
    colnames(df) <- c("z", group_by)
    aes <- ggplot2::aes_string(xintercept = "z")
    h <- h + ggplot2::geom_vline(aes,
      na.rm = TRUE,
      data = df,
      color = linecolor,
      linetype = linetype,
      lwd = lwd
    )
  }
  return(h)
}

# Add line layer to a componentwise plot
add_c_layer_line <- function(h, df, x_name, group_by, color_by, alpha) {
  if (!is.null(df)) {
    check_type(df, "data.frame")
    if ("_draw_" %in% colnames(df)) {
      df <- add_factor_crossing(df, group_by, "_draw_", "group_x_draw")
      group_by <- "group_x_draw"
    }
    aes_fit <- ggplot2::aes_string(
      x = x_name,
      y = "y",
      color = color_by,
      group = group_by
    )
    h <- h + ggplot2::geom_line(
      data = df,
      mapping = aes_fit,
      inherit.aes = FALSE,
      alpha = alpha
    )
    if (!is.null(color_by)) {
      cfac <- dollar(df, color_by)
      num_colors <- length(levels(cfac))
      h <- h + scale_color(num_colors)
    }
  }
  return(h)
}

# Add ribbon layer to a componentwise plot
add_c_layer_ribbon <- function(h, df, x_name, group_by, color_by, alpha) {
  if (!is.null(df)) {
    check_type(df, "data.frame")
    aes_rib <- ggplot2::aes_string(
      x = x_name,
      ymin = "lower",
      ymax = "upper",
      group = group_by,
      color = color_by,
      fill = color_by
    )
    h <- h + ggplot2::geom_ribbon(
      data = df,
      mapping = aes_rib,
      inherit.aes = FALSE,
      alpha = alpha,
      color = NA # no ribbon edge
    )
    if (!is.null(color_by)) {
      cfac <- dollar(df, color_by)
      num_colors <- length(levels(cfac))
      h <- h + scale_fill(num_colors)
    }
  }
  return(h)
}

# Add ribbon layer to a grouped plot
add_g_layer_ribbon <- function(h, df, x_name, group_by, color, alpha) {
  if (!is.null(df)) {
    check_type(df, "data.frame")
    aes_rib <- ggplot2::aes_string(
      x = x_name,
      ymin = "lower",
      ymax = "upper",
      group = group_by
    )
    h <- h + ggplot2::geom_ribbon(
      data = df,
      mapping = aes_rib,
      inherit.aes = FALSE,
      color = NA,
      fill = color,
      alpha = alpha
    )
  }
  return(h)
}

# Add model fit (line) layer to a grouped plot
add_g_layer_fit <- function(h, df, x_name, group_by, color, alpha) {
  if (!is.null(df)) {
    check_type(df, "data.frame")
    if ("_draw_" %in% colnames(df)) {
      df <- add_factor_crossing(df, group_by, "_draw_", "group_x_draw")
      group_by <- "group_x_draw"
    }
    aes_fit <- ggplot2::aes_string(
      x = x_name,
      y = "y",
      group = group_by
    )
    h <- h + ggplot2::geom_line(
      data = df,
      mapping = aes_fit,
      inherit.aes = FALSE,
      color = color,
      alpha = alpha
    )
  }
  return(h)
}

# Add true signal (line) layer to a grouped plot
add_g_layer_signal <- function(h, df, x_name, group_by, color) {
  if (!is.null(df)) {
    check_type(df, "data.frame")
    aes_sig <- ggplot2::aes_string(
      x = x_name,
      y = "y",
      group = group_by
    )
    h <- h + ggplot2::geom_line(
      data = df,
      mapping = aes_sig,
      inherit.aes = FALSE,
      color = color
    )
  }
  return(h)
}

# Add faceting to a grouped plot
add_g_layer_faceting <- function(h, facet_by, nrow, ncol) {
  f <- stats::as.formula(paste("~", facet_by))
  labfun <- ggplot2::label_both
  h <- h + ggplot2::facet_wrap(f, nrow = nrow, ncol = ncol, labeller = labfun)
  return(h)
}


# Plot colors
#
# @param main Color name. Must be a valid scheme name for
# \code{\link[bayesplot]{color_scheme_get}}.
# @param variant Must be one of {"light", "light_highlight", "mid",
# "mid_highlight", "dark", "dark_highlight"}.
# @return A hex value of the color.
colorset <- function(main, variant = "mid") {
  scheme <- bayesplot::color_scheme_get(scheme = main)
  col <- scheme[[variant]]
  if (is.null(col)) {
    stop("Invalid color!")
  }
  return(col)
}

# A color palette for lines etc
color_palette <- function(n) {
  c1 <- colorset("brightblue", "mid_highlight")
  c2 <- colorset("red", "mid_highlight")
  c3 <- colorset("orange", "mid_highlight")
  c4 <- colorset("green", "mid_highlight")
  c5 <- colorset("gray", "dark_highlight")
  palette <- c(c1, c2, c3, c4, c5)
  palette[1:n]
}

# A color palette for fills
fill_palette <- function(n) {
  c1 <- colorset("brightblue", "mid")
  c2 <- colorset("red", "mid")
  c3 <- colorset("orange", "mid")
  c4 <- colorset("green", "mid")
  c5 <- colorset("gray", "dark")
  palette <- c(c1, c2, c3, c4, c5)
  palette[1:n]
}

# A color scale for lines etc
scale_color <- function(n) {
  if (n > 5) {
    return(NULL)
  }
  values <- color_palette(n)
  ggplot2::scale_color_manual(values = values)
}

# A color scale for fills
scale_fill <- function(n) {
  if (n > 5) {
    return(NULL)
  }
  values <- fill_palette(n)
  ggplot2::scale_fill_manual(values = values)
}


#' Visualize input warping function with several steepness parameter values
#'
#' @param wrp a vector of values of the warping steepness parameter
#' @param x a vector of input values
#' @param color line color
#' @param alpha line alpha
#' @return a \code{ggplot} object
plot_inputwarp <- function(wrp,
                           x,
                           color = colorset("red", "dark"),
                           alpha = 0.5) {
  x <- sort(x)
  L <- length(x)
  S <- length(wrp)
  W <- matrix(0, S, L)
  for (i in seq_len(S)) {
    w <- warp_input(x, a = wrp[i])
    W[i, ] <- w
  }
  af <- as.factor(rep(1:S, each = L))
  df <- data.frame(rep(x, S), as.vector(t(W)), af)
  colnames(df) <- c("x", "w", "idx")

  # Create ggplot object
  aes <- ggplot2::aes_string(x = "x", y = "w", group = "idx")
  plt <- ggplot2::ggplot(df, aes)

  # Add titles and labels
  plt <- plt + ggplot2::labs(x = "Input", y = "Warped input") +
    ggplot2::ggtitle("Input-warping function")
  plt <- plt + ggplot2::ylim(-1.0, 1.0)

  # Plot the actual lines of interest
  plt <- plt + ggplot2::geom_line(color = color, alpha = alpha)
  return(plt)
}
