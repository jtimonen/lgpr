#' Compute predictions for a fitted model
#'
#' @export
#' @description Compute predictions for a fitted model. Only possible for models
#' with Gaussian likelihood.
#' @param fit An object of class \code{lgpfit}.
#' @param X_test The test points where the predictions should be computed.
#' @param samples The predictions can be computed either by using only the
#' posterior mean \cr (\code{samples="mean"}), median (\code{samples="median"}),
#' or MAP (\code{samples="map"}) parameters, or for all parameter samples
#' (\code{samples="all"}). This can also be a set of indices, for example
#' \code{samples=c(1:10)} gives predictions for the parameter samples 1...10.
#' @param print_progress Should progress be printed (if there is more than one
#' sample)?
#' @param print_params Should the parameter values be printed? (only works if
#' \code{samples} is mean or median.)
#' @return A list.
#' @seealso
#' \itemize{
#'   \item For creating an \code{lgpfit} object, see \code{\link{lgp_fit}}.
#'   \item For creating an \code{lgpmodel} object, see \code{\link{lgp_model}}.
#' }
lgp_predict <- function(fit,
                        X_test,
                        samples = "map",
                        print_progress = TRUE,
                        print_params = FALSE) {

  # Run some input checks and scale data correctly
  PP <- predict_preproc(fit, X_test, samples)
  X_data <- PP$X_data
  X_test <- PP$X_test
  y_data <- PP$y_data
  samples <- PP$samples
  info <- PP$info
  D <- PP$D
  cnames <- PP$cnames
  LIST <- list()
  TSCL <- fit@model@scalings$TSCL

  cat("* Computing predictions ")

  if (is.character(samples)) {

    # Predictions using a single parameter estimate
    if (samples %in% c("mean", "median", "map")) {
      if (samples %in% c("mean", "median")) {
        samples_str <- paste("posterior", samples)
      } else {
        samples_str <- "MAP"
      }
      cat("using ", samples_str, " parameters. \n", sep = "")
      params <- hyperparam_estimate(fit, samples)
      if (print_params) {
        print(params)
      }
      LIST[[1]] <- compute_predictions(
        X_data, y_data, X_test, params, D, info,
        cnames, TSCL
      )
    } else {
      stop("\ninvalid 'samples' input for lgp_predict! (", samples, ")")
    }
  } else {

    # Predictions using multiple samples
    cat("using ", length(samples), " parameter sample(s). \n\n", sep = "")
    PAR <- hyperparam_samples(fit, samples)
    pnames <- colnames(PAR)
    ns <- dim(PAR)[1]

    # Progress bar top
    if (print_progress) {
      str <- paste("|   ", seq(10, 100, by = 10), "%", sep = "")
      top <- paste(formatC(str, width = 4), collapse = " ")
      top <- paste(top, "|")
      barlen <- nchar(top) - 1
      iprint <- ceiling(seq(1, ns, length.out = barlen))
      if (ns >= 100) {
        cat(top, "\n ")
      }
    }

    for (i_smp in 1:ns) {
      # Predict with current parameter sample
      params <- as.numeric(PAR[i_smp, ])
      names(params) <- pnames
      LIST[[i_smp]] <- compute_predictions(
        X_data, y_data, X_test, params, D,
        info, cnames, TSCL
      )

      # Print progress
      if (print_progress && (ns > 1)) {
        if (i_smp %in% iprint) {
          cat("=")
        }
        if (i_smp == ns) {
          cat("\n\n")
        }
      }
    }
  }

  # Return
  ret <- list(LIST = LIST, X_test_scaled = X_test)
  return(ret)
}


#' Compute predictions and log-posterior predictive density at test points
#'
#' @export
#' @description This is a convenience function that wraps
#' \code{\link{lgp_predict}}, \code{\link{compute_lppd}} and
#' \code{\link{plot_posterior_y}}.
#' @param fit an object of class \code{lgpfit}
#' @param test_data a test data matrix
#' @param verbose Should this print progress?
#' @param samples Sample indices or a keyword "mean", "median", "map", or "all".
#' @param plot should this return also a plot of the data and predictions?
#' @return a ggplot object or lppd
lgp_test <- function(fit, test_data, plot = FALSE, verbose = TRUE,
                     samples = "mean") {

  # predict only at test points
  info <- fit@model@info
  yname <- info$response_name
  xnames <- info$covariate_names
  idvar <- info$varInfo$id_variable
  tvar <- info$varInfo$time_variable
  xnames <- union(c(idvar, tvar), xnames)
  dat <- data.frame(test_data)
  X_test <- dat[, xnames]
  y_test <- dat[, yname]
  PRED <- lgp_predict(fit, X_test, samples = samples, print_progress = verbose)
  LPPD <- compute_lppd(PRED, y_test)
  mlppd <- mean(as.numeric(LPPD))
  ret <- list(lppd = LPPD, mlppd = mlppd)
  if (plot) {
    p <- plot_posterior_y(fit, PRED,
      test_data = test_data,
      uncertainty = "errorbar"
    )
    subt <- paste("Mean log-posterior predictive density:", round(mlppd, 5))
    p <- p + ggplot2::ggtitle(
      label = "Predictive distribution of y at test points",
      subtitle = subt
    )
    ret$plot <- p
  }
  return(ret)
}
