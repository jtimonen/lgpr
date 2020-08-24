#' Printing S4 object info using the show generic
#'
#' @param object an object of some S4 class
#' @return the object invisibly
#' @name show

#' @rdname show
setMethod(
  f = "show", signature = "lgpformula",
  definition = function(object) {
    cat(as.character(object))
    invisible(object)
  }
)

#' @rdname show
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

#' @rdname show
setMethod(
  f = "show", signature = "lgpmodel",
  definition = function(object) {
    cat(as.character(object))
    invisible(object)
  }
)

#' @rdname show
setMethod(
  f = "show", signature = "lgpsim",
  definition = function(object) {
    cat(as.character(object))
    invisible(object)
  }
)
