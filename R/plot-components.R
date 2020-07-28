
#' Helper function for plotting components
#'
#' @inheritParams plot_component
#' @param ncol number of plot columns
#' @param nrow number of plot rows
#' @param legend legend argument for ggarrange, use \code{"none"} to remove legends
#' @param labels labels argument for ggarrange
#' @param ylim y axis limits
#' @param font_size font size for plots
#' @param theme ggplot theme
#' @param legend_dir direction of legend
#' @param xlabel x-axis label
#' @param ylabel y-axis label
#' @param return_list should this return a list of ggplot objects
#' instead of doing ggarrange
#' @return an object returned by \code{ggpubr::ggarrange} or a list
plot_components <- function(MMM, SSS, model, time_is_xvar,
                            X_test = NULL,
                            sum_highlight = NULL,
                            linealpha = 1,
                            linetype = 1,
                            fill_alpha = 0.3,
                            marker = NULL,
                            ncol = NULL,
                            nrow = NULL,
                            legend = NULL,
                            labels = NULL,
                            ylim = NULL,
                            font_size = 9,
                            theme = ggplot2::theme_linedraw(),
                            legend_dir = "horizontal",
                            xlabel = NULL,
                            ylabel = " ",
                            viridis_option = "viridis",
                            return_list = FALSE) {
  GG <- list()
  sum_D <- sum(model@stan_dat$D) + 1
  for (d in 1:sum_D) {
    gg <- plot_component(
      MMM, SSS, model, d, time_is_xvar,
      linealpha, linetype, fill_alpha,
      X_test, marker, sum_highlight, viridis_option
    )
    if (is.null(ylim)) {
      if (!is.null(SSS)) {
        ylim <- c(min(MMM - SSS), max(MMM + SSS))
      } else {
        ylim <- range(MMM)
      }
    }
    gg <- gg + theme + ggplot2::ylim(ylim)
    if (!is.null(xlabel)) {
      gg <- gg + ggplot2::xlab(xlabel)
    }
    if (!is.null(ylabel)) {
      gg <- gg + ggplot2::ylab(ylabel)
    }
    if (is.null(legend)) {
      gg <- gg + ggplot2::theme(
        legend.justification = c(0.95, 0.05),
        legend.position = c(0.95, 0.05),
        legend.key = ggplot2::element_rect(fill = "gray95"),
        legend.background = ggplot2::element_rect(fill = "gray95"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        legend.direction = legend_dir
      )
    }
    GG[[d]] <- gg + ggplot2::theme(text = ggplot2::element_text(size = font_size))
  }
  if (return_list) {
    p <- GG
  } else {
    p <- ggpubr::ggarrange(
      plotlist = GG, ncol = ncol, nrow = nrow,
      legend = legend,
      labels = labels
    )
  }
  return(p)
}

#' Helper function for plotting one component
#' @param MMM a n array of size n_samples x n_data x n_components
#' @param SSS a n array of size n_samples x n_data x n_components
#' @param model an object of class 'lgpmodel'
#' @param idx Index of component to be plotted.
#' @param time_is_xvar is the time variable the x-axis variable
#' @param linealpha line alpha
#' @param linetype line type
#' @param fill_alpha fill alpha for geom_ribbons
#' @param X_test optional matrix of test points
#' @param marker point type
#' @param sum_highlight name of a categorical covariate to be highlighted
#' @param viridis_option the option argument of \code{ggplot2::scale_colour_viridis_c}
#' by colour in the sum plot
#' @return a ggplot object
plot_component <- function(MMM, SSS, model, idx, time_is_xvar,
                           linealpha, linetype, fill_alpha,
                           X_test, marker, sum_highlight, viridis_option) {
  sdat <- model@stan_dat
  if (is.null(X_test)) {
    X <- t(sdat$X)
    Xnn <- sdat$X_notnan
  } else {
    X <- X_test
    Xnn <- as.numeric(!is.nan(X_test[, 3]))
  }
  D <- sdat$D
  SCL <- model@scalings
  sum_D <- sum(D) + 1
  S <- dim(MMM)[1]
  n <- dim(MMM)[2]
  d <- dim(MMM)[3]
  cpn <- model@info$component_names
  cvn <- model@info$covariate_names
  if (sum_D != d) {
    stop("dim(MMM)[3] must be ", sum_D)
  }
  if (idx < 1 || idx > sum_D) {
    stop("idx must be between 1 and ", sum_D)
  }
  if (idx < sum_D) {
    ctype <- component_index_to_type(D, idx)
    cind <- component_index_to_covariate_index(D, idx)
  } else {
    ctype <- 7 # sum
  }

  plot_eb <- !is.null(SSS) && !(ctype %in% c(1, 7))
  MMM_i <- MMM[, , idx]
  f <- as.numeric(t(MMM_i))
  if (plot_eb) {
    UUU_i <- MMM_i + SSS[, , idx]
    LLL_i <- MMM_i - SSS[, , idx]
    ub <- as.numeric(t(UUU_i))
    lb <- as.numeric(t(LLL_i))
  }
  sample <- rep(1:S, each = n)
  id <- rep(X[, 1], S)


  if (ctype == 4) {
    icnt <- idx_to_cont_index(D, idx)
    cscl <- model@scalings$CSCL[[icnt]]
  }

  if (ctype %in% c(1, 2, 5, 6, 7)) {
    time_is_xvar <- TRUE
  }
  if (time_is_xvar) {
    xvar <- rep(SCL$TSCL$fun_inv(X[, 2]), S)
    xvar_name <- model@info$varInfo$time_variable
  } else {
    xvar <- rep(X[, cind], S)
    xvar_name <- cvn[cind]
    if (ctype == 4) {
      xvar <- cscl$fun_inv(xvar)
    }
  }



  # Create df and aes
  if (ctype == 1) {
    leg <- ggplot2::theme(legend.position = "none")
    grpvar <- as.factor(paste(id, sample))
    df <- data.frame(xvar, f, grpvar)
    aes <- ggplot2::aes_string(
      x = "xvar", y = "f",
      group = "grpvar"
    )
    aes_eb <- ggplot2::aes_string(ymin = "lb", ymax = "ub")
  } else if (ctype == 2) {
    leg <- ggplot2::theme(legend.position = "none")
    grpvar <- as.factor(sample)
    df <- data.frame(xvar, f, grpvar)
    aes <- ggplot2::aes_string(x = "xvar", y = "f", group = "grpvar")
    aes <- ggplot2::aes_string(
      x = "xvar", y = "f",
      group = "grpvar"
    )
    aes_eb <- ggplot2::aes_string(ymin = "lb", ymax = "ub")
  } else if (ctype == 3) {
    colorvar <- as.factor(rep(Xnn, S))
    leg <- ggplot2::labs(color = "Group", fill = "Group")
    grpvar <- as.factor(paste(id, sample))
    df <- data.frame(xvar, f, grpvar, colorvar)
    aes <- ggplot2::aes_string(
      x = "xvar", y = "f",
      group = "grpvar", color = "colorvar"
    )
    aes_eb <- ggplot2::aes_string(ymin = "lb", ymax = "ub", fill = "colorvar")
  } else if (ctype == 4) {
    leg <- ggplot2::theme(legend.position = "none")
    grpvar <- as.factor(paste(id, sample))
    colorvar <- cscl$fun_inv(rep(X[, cind], S))
    df <- data.frame(xvar, f, grpvar, colorvar)
    aes <- ggplot2::aes_string(
      x = "xvar", y = "f",
      group = "grpvar",
      color = "colorvar"
    )
    aes_eb <- ggplot2::aes_string(ymin = "lb", ymax = "ub")
  } else if (ctype == 5 || ctype == 6) {
    leg <- ggplot2::labs(color = cvn[cind], fill = cvn[cind])
    colorvar <- rep(X[, cind], S)
    grpvar <- as.factor(paste(colorvar, sample))
    colorvar <- as.factor(colorvar)
    df <- data.frame(xvar, f, grpvar, colorvar)
    aes <- ggplot2::aes_string(
      x = "xvar", y = "f",
      group = "grpvar", color = "colorvar"
    )
    aes_eb <- ggplot2::aes_string(ymin = "lb", ymax = "ub", fill = "colorvar")
  } else if (ctype == 7) {
    if (!is.null(sum_highlight)) {
      if (sum_highlight == "group") {
        colorvar <- as.factor(rep(Xnn, S))
        leg <- ggplot2::labs(color = "Group", fill = "Group")
        grpvar <- as.factor(paste(id, sample))
        df <- data.frame(xvar, f, grpvar, colorvar)
        aes <- ggplot2::aes_string(
          x = "xvar", y = "f",
          group = "grpvar",
          color = "colorvar"
        )
        aes_eb <- ggplot2::aes_string(
          ymin = "lb", ymax = "ub",
          fill = "colorvar"
        )
      } else {
        cind <- which(cvn == sum_highlight)
        if (length(cind) == 0) {
          stop("invalid sum_highlight (", sum_highlight, ")")
        }
        leg <- ggplot2::labs(color = cvn[cind], fill = cvn[cind])
        colorvar <- rep(X[, cind], S)
        grpvar <- as.factor(paste(id, sample))
        colorvar <- as.factor(colorvar)
        df <- data.frame(xvar, f, grpvar, colorvar)
        aes <- ggplot2::aes_string(
          x = "xvar", y = "f",
          group = "grpvar",
          color = "colorvar"
        )
        aes_eb <- ggplot2::aes_string(
          ymin = "lb", ymax = "ub",
          fill = "colorvar"
        )
      }
    } else {
      leg <- ggplot2::theme(legend.position = "none")
      grpvar <- as.factor(paste(id, sample))
      df <- data.frame(xvar, f, grpvar)
      aes <- ggplot2::aes_string(x = "xvar", y = "f", group = "grpvar")
      aes_eb <- ggplot2::aes_string(ymin = "lb", ymax = "ub")
    }
  } else {
    stop("invalid ctype!")
  }

  # Create title
  if (ctype != 7) {
    title <- parse(text = cpn[idx])
  } else {
    title <- parse(text = "f[sum]")
  }

  # Create ggplot object
  h <- ggplot2::ggplot(df, aes)
  if (plot_eb && (fill_alpha > 0)) {
    h <- h + ggplot2::geom_ribbon(aes_eb,
      alpha = fill_alpha,
      lty = 0
    )
  }
  if (ctype == 4) {
    h <- h + ggplot2::scale_colour_viridis_c(option = viridis_option)
  }

  h <- h + leg
  h <- h + ggplot2::geom_line(linetype = linetype, alpha = linealpha)
  h <- h + ggplot2::ggtitle(title)
  h <- h + ggplot2::xlab(xvar_name)

  # Add marker
  if (!is.null(marker)) {
    h <- h + ggplot2::geom_point(pch = marker)
  }
  return(h)
}


#' Component index to component type
#' @param D integer vector of length 6
#' @param idx integer
#' @return an integer
component_index_to_type <- function(D, idx) {
  all <- c(
    rep(1, D[1]), rep(2, D[2]), rep(3, D[3]),
    rep(4, D[4]), rep(5, D[5]), rep(6, D[6])
  )
  return(all[idx])
}

#' Component index to covariate index
#' @param D integer vector of length 6
#' @param idx integer
#' @return an integer
component_index_to_covariate_index <- function(D, idx) {
  covr_idx <- 2 - sum(D[1:2]) + idx
  return(covr_idx)
}


#' Component index to how manyth continuous covariate it is
#' @param D integer vector of length 6
#' @param idx an integer
#' @return an integer
idx_to_cont_index <- function(D, idx) {
  type <- component_index_to_type(D, idx)
  if (type != 4) {
    stop("This is not a continuous covariate")
  }
  all <- c(
    rep(1, D[1]), rep(2, D[2]), rep(3, D[3]),
    rep(4, D[4]), rep(5, D[5]), rep(6, D[6])
  )
  n_moved <- 1
  for (j in 1:D[4]) {
    if (all[idx - j] == 4) {
      n_moved <- n_moved + 1
    }
  }
  return(n_moved)
}
