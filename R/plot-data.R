#' Plotting longitudinal data
#'
#' @export
#' @param data A data frame.
#' @param x_name Name of x-axis variable.
#' @param y_name Name of the y-axis variable.
#' @param group_by Name of grouping variable (must be a factor).
#' @param color_by Name of coloring variable (must be a factor).
#' @param facet_by Name of the faceting variable (must be a factor).
#' @param highlight Value of category of the \code{group_by} variable
#' that is highlighted. Can only be used if \code{color_by} is \code{NULL}.
#' @param main main plot title
#' @param sub plot subtitle
#' @param fit An object of class \linkS4class{lgpfit}. If this is provided,
#' also the model fit is plotted.
#' @return a \code{ggplot} object
plot_data <- function(data,
                      x_name = "age",
                      y_name = "y",
                      group_by = "id",
                      facet_by = NULL,
                      color_by = NULL,
                      highlight = NULL,
                      main = NULL,
                      sub = NULL,
                      fit = NULL) {
  titles <- plot_data_titles(main, sub, data, group_by)
  df <- data[c(x_name, y_name, group_by, color_by, facet_by)]
  df <- plot_data_add_highlight_factor(df, group_by, highlight)
  if (!is.null(highlight)) {
    color_by <- paste0(group_by, "_")
  }
  aes <- plot_data_aes(x_name, y_name, group_by, color_by)
  h <- ggplot2::ggplot(df, aes)
  h <- h + ggplot2::geom_line() + ggplot2::geom_point()
  h <- h + ggplot2::ggtitle(label = titles$main, subtitle = titles$sub)
  if (!is.null(facet_by)) {
    f <- stats::as.formula(paste("~", facet_by))
    h <- h + ggplot2::facet_wrap(f)
  }
  num_colors <- plot_data_num_colors(df, color_by)
  if (num_colors <= 4) {
    values <- color_palette(num_colors)
    h <- h + ggplot2::scale_color_manual(values = values)
  }
  h <- plot_data_add_fit(h, fit)

  return(h)
}

#' Add model fit to data plot
#'
#' @inheritParams plot_data
#' @param h the current \code{ggplot} object
#' @return a modified \code{ggplot} object
plot_data_add_fit <- function(h, fit){
  if (!is.null(fit)) {
    print(fit)
  }
  return(h)
}


#' Create aes for plot_data
#'
#' @inheritParams plot_data
#' @return an \code{aes} object
plot_data_aes <- function(x_name, y_name, group_by, color_by) {
  if (is.null(color_by)) {
    aes <- ggplot2::aes_string(x = x_name, y = y_name, group = group_by)
  } else {
    aes <- ggplot2::aes_string(
      x = x_name, y = y_name, group = group_by,
      color = color_by
    )
  }
  return(aes)
}

#' Create titles for plot_data
#'
#' @inheritParams plot_data
#' @return a list
plot_data_titles <- function(main, sub, data, group_by) {
  g <- data[[group_by]]
  N <- length(unique(g))
  n <- dim(data)[1]
  if (is.null(main)) {
    main <- paste0(n, " data points")
  }
  if (is.null(sub)) {
    sub <- paste0(N, " different levels for ", group_by)
  }
  list(main = main, sub = sub, N = N, n = n)
}


#' Get number of distinct colors needsd by plot_data
#'
#' @inheritParams plot_data
#' @return a list
plot_data_num_colors <- function(data, color_by) {
  if (is.null(color_by)) {
    N <- 1
  } else {
    g <- data[[color_by]]
    N <- length(unique(g))
  }
  return(N)
}

#' Add factor to data frame for highlighting in plot
#'
#' @inheritParams plot_data
#' @param df data frame
#' @return a list
plot_data_add_highlight_factor <- function(df, group_by, highlight) {
  if (!is.null(highlight)) {
    g <- df[[group_by]]
    hl <- 1 + as.numeric(g == highlight)
    levels <- c("other", highlight)
    name <- paste0(group_by, "_")
    df[[name]] <- as.factor(levels[hl])
  }
  return(df)
}
