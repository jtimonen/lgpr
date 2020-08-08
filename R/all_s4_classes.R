#' An S4 class to represent an lgp term with a single covariate
#'
#' @slot covariate name of a covariate
#' @slot type covariate type
#' @slot fun kernel function name
lgpterm <- setClass("lgpterm",
  slots = c(
    covariate = "character",
    type = "character",
    fun = "character"
  )
)

#' Convert an \code{lgpterm} to character
#'
#' @export
#' @param object an object of class \code{lgpterm}
#' @rdname as.character_lgpterm
setMethod(
  f = "as.character", signature = "lgpterm",
  definition = function(x) {
    desc <- paste0(x@fun, "(", x@covariate, ")")
    return(desc)
  }
)

#' An S4 class to represent a product of lgp terms
#'
#' @slot factors a list of at most two \code{\link{lgpterm}}s
lgpproduct <- setClass("lgpproduct", slots = c(factors = "list"))


#' Convert an \code{lgpproduct} to character
#'
#' @export
#' @param object an object of class \code{lgpproduct}
#' @rdname as.character_lgpproduct
setMethod(
  f = "as.character", signature = "lgpproduct",
  definition = function(x) {
    facs <- x@factors
    L <- length(facs)
    c1 <- as.character(facs[[1]])
    if (L == 1) {
      desc <- paste0("(1st order):   ", c1)
    } else {
      c2 <- as.character(facs[[2]])
      desc <- paste0("(interaction): ", c1, " * ", c2)
    }
    return(desc)
  }
)


#' An S4 class to represent the right-hand side of an lgp formula
#'
#' @slot summands a list of one or more \code{\link{lgpproduct}}s
lgpsum <- setClass("lgpsum", slots = c(summands = "list"))

#' Convert an \code{lgpsum} to character
#'
#' @export
#' @param object an object of class \code{lgpsum}
#' @rdname as.character_lgpsum
setMethod(
  f = "as.character", signature = "lgpsum",
  definition = function(x) {
    s <- x@summands
    L <- length(s)
    desc <- ""
    for (j in seq_len(L)) {
      desc <- paste0(desc, "Term ", j, " ", as.character(s[[j]]), "\n")
    }
    return(desc)
  }
)

#' Show a summary of an \code{lgpsum}
#'
#' @export
#' @param object an object of class \code{lgpsum}
#' @rdname show_lgpsum
setMethod(
  f = "show", signature = "lgpsum",
  definition = function(object) {
    cat(as.character(object))
    invisible(object)
  }
)

#' An S4 class to represent an lgp formula
#'
#' @slot components an object of class \code{\link{lgpsum}}
#' @slot response name of the response variable
lgpformula <- setClass("lgpformula",
  slots = c(
    call = "character",
    response = "character",
    components = "lgpsum"
  )
)

#' Show a summary of an \code{lgpformula}
#'
#' @export
#' @param object an object of class \code{lgpformula}
#' @rdname show_lgpformula
setMethod(
  f = "show", signature = "lgpformula",
  definition = function(object) {
    desc <- paste0("Formula: ", object@call, "\n")
    desc <- paste0("Response variable: ", object@response, "\n")
    desc <- paste0(desc, as.character(object@components))
    cat(desc)
    invisible(object)
  }
)


#' An S4 class to represent an lgp model
#'
#' @slot data The original unmodified data frame.
#' @slot stan_dat The data to be given as input to \code{rstan::sampling}.
#' @slot scalings Preprocessing scaling functions and their inverse operations.
#' @slot info Model info.
lgpmodel <- setClass("lgpmodel",
  slots = c(
    data = "data.frame",
    stan_dat = "list",
    scalings = "list",
    info = "list"
  )
)


#' Show a summary of an \code{lgpmodel}
#'
#' @export
#' @param object an object of class \code{lgpmodel}
#' @rdname show_lgpmodel
#' @return nothing
setMethod(
  f = "show",
  signature = "lgpmodel",
  definition = function(object) {
    print_prior <- TRUE
    if (object@stan_dat$UNCRT == 1) {
      str2 <- "uncertain"
    } else {
      str2 <- "fixed"
    }
    if (object@stan_dat$HMGNS == 0) {
      str3 <- "heterogeneous"
    } else {
      str3 <- "homogeneous"
    }
    D <- object@stan_dat$D
    cat("\n ---------- LGPMODEL SUMMARY ----------\n\n")
    str_dis <- paste("  Disease component \n",
      "    - Effect time: ", str2, "\n",
      "    - Effect type: ", str3, "\n",
      sep = ""
    )

    minfo <- model_info(object, print = FALSE)
    cat(minfo)
    if (D[3] == 1) {
      cat(str_dis)
    }
    cat("\n")
    if (print_prior) {
      prior_info <- prior_stan_to_readable(object@stan_dat)
      cat(prior_info)
    }
    invisible(object)
  }
)

#' An S4 class to represent the output of the \code{lgp_fit} function
#'
#' @slot stan_fit The \code{stanfit} object returned by \code{rstan::sampling}.
#' @slot model The \code{lgpmodel} object returned by \code{lgp_model}.
#' @slot relevances Inferred component relevances.
#' @slot selection Component selection info.
#' @slot pkg_version Package version number.
#' @slot diagnostics  A data frame with columns
#' \code{c("Rhat", "Bulk_ESS", "Tail_ESS")}.
lgpfit <- setClass("lgpfit",
  slots = c(
    stan_fit = "stanfit",
    model = "lgpmodel",
    relevances = "list",
    selection = "list",
    pkg_version = "character",
    diagnostics = "data.frame"
  )
)


#' Show a summary of results of the \code{lgp} function
#'
#' @export
#' @param object an object of class \code{lgpfit}
#' @rdname show_lgpfit
#' @return nothing
setMethod(
  f = "show",
  signature = "lgpfit",
  definition = function(object) {
    cat("\n ---------- LGPFIT SUMMARY ----------\n\n")

    # Runtime info
    tryCatch(
      {
        runtime <- get_runtime(object)
        cat("* Average runtime per chain: ",
          runtime$warmup, " s (warmup) and ",
          runtime$sampling, " s (sampling)\n",
          sep = ""
        )
      },
      error = function(e) {
        cat("* Unable to show runtime info. Reason:\n")
        print(e)
      }
    )

    # Convergence info
    tryCatch(
      {
        diag <- object@diagnostics
        Rhat <- diag$Rhat
        imax <- which(Rhat == max(Rhat))
        cat("* Largest R-hat value is ", round(max(Rhat), 4),
          " (", paste(rownames(diag)[imax], collapse = ", "),
          ")\n",
          sep = ""
        )
      },
      error = function(e) {
        cat("* Unable to show convergence info. Reason:\n")
        print(e)
      }
    )

    # Relevance info
    tryCatch(
      {
        sel <- object@selection
        rel_method <- object@relevances$method
        r3 <- sel$prob
        tr <- sel$threshold
        cat("* Used relevance method = ", rel_method, "\n", sep = "")
        cat("* Used selection threshold = ", tr, "\n", sep = "")

        rel <- object@relevances$average
        cn <- names(rel)
        r1 <- round(rel, 3)
        r2 <- cn %in% sel$selected
        DF <- data.frame(r1, r2, r3)
        cat("\n")

        colnames(DF) <- c("Relevance", "Selected", "Prob.")
        print(DF)
        cat("\n")
      },
      error = function(e) {
        cat("* Unable to show component relevance info. Reason:\n")
        print(e)
      }
    )
  }
)

#' Visualize a fitted `lgpfit` object
#'
#' @export
#' @param fit an object of class \code{lgpfit}
#' @param color_scheme bayesplot color scheme
#' @param x does nothing
#' @param y does nothing
#' @rdname plot_lgpfit
#' @return a ggplot object
setMethod(
  f = "plot",
  signature = "lgpfit",
  definition = function(fit, x = 1, y = 1, color_scheme = "red") {
    h <- plot_relevances(fit, color_scheme = color_scheme)
    return(h)
  }
)
