#' Visualize a model fit
#'
#' @description
#' \itemize{
#'   \item In \code{plot_fit}, the data and model fit are plotted for
#'   each individual separately using \code{\link{plot_api_g}}
#'   \item In \code{plot_fit_component}, one component of the model fit
#'   is plotted using \code{\link{plot_api_c}}
#' }
#' @inheritParams get_draws
#' @inheritParams create_plot_df
#' @inheritParams plot_sim
#' @param draws see the \code{draws} argument of \code{\link{get_f}}
#' @return a \code{\link[ggplot2]{ggplot}} object
#' @family main plot functions
#' @name plot_fit
NULL

#' @export
#' @rdname plot_fit
plot_fit <- function(fit,
                     group_by = "id",
                     x_name = "age",
                     draws = NULL,
                     ...) {
  check_type(fit, "lgpfit")
  df_data <- create_plot_df(fit, x_name, group_by)
  DF <- plot_fit_helper(fit, df_data, draws)
  h <- plot_api_g(
    df_data = df_data,
    df_fit =  dollar(DF, "df_fit"),
    df_fit_err = dollar(DF, "df_ribbon"),
    fit_alpha = dollar(DF, "fit_alpha"),
    ...
  )
  h <- h + ggplot2::ggtitle("Model fit", subtitle = DF$info)
  return(h)
}

#' @export
#' @rdname plot_fit
plot_fit_component <- function(fit,
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
