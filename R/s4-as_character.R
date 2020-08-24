#' Character representations of different S4 objects
#'
#' @param x an object of some S4 class
#' @return a character representation of the object
#' @name as_character
NULL

#' @rdname as_character
setMethod(
  f = "as.character", signature = "lgpexpr",
  definition = function(x) {
    paste0(x@fun, "(", x@covariate, ")")
  }
)

#' @rdname as_character
setMethod(
  f = "as.character", signature = "lgpterm",
  definition = function(x) {
    term_as_character(x, verbose = TRUE)
  }
)

#' @rdname as_character
setMethod(
  f = "as.character", signature = "lgprhs",
  definition = function(x) {
    s <- x@summands
    L <- length(s)
    desc <- ""
    for (j in seq_len(L)) {
      desc <- paste0(
        desc, "Term ", j, " expressions:   ",
        as.character(s[[j]]), "\n"
      )
    }
    return(desc)
  }
)


#' @rdname as_character
setMethod(
  f = "as.character", signature = "lgpformula",
  definition = function(x) {
    desc <- "An object of class 'lgpformula'.\n\n"
    desc <- paste0(desc, "Call:                 ", x@call, "\n")
    desc <- paste0(desc, "Response variable:    ", x@y_name, "\n")
    desc <- paste0(desc, as.character(x@terms))
    return(desc)
  }
)

#' @rdname as_character
setMethod(
  f = "as.character", signature = "lgpmodel",
  definition = function(x) {
    str1 <- covariate_names(x, type = "continuous")
    str2 <- covariate_names(x, type = "categorical")
    desc <- "An object of class lgpmodel.\n"
    desc <- paste0(desc, "\nFormula: ", x@model_formula@call)
    desc <- paste0(desc, "\nContinuous: {", str1, "}")
    desc <- paste0(desc, "\nCategorical: {", str2, "}")
    return(desc)
  }
)

#' @rdname as_character
setMethod(
  f = "as.character", signature = "lgpsim",
  definition = function(x) {
    desc <- "An object of class lgpsim with slots "
    str <- paste(slotNames(x), collapse = ", ")
    desc <- paste0(desc, "{", str, "}")
    return(desc)
  }
)

#' @rdname as_character
setMethod(
  f = "as.character", signature = "lgpfit",
  definition = function(x) {
    desc <- "\n ---------- LGPFIT SUMMARY ----------\n\n"
    
    # Runtime info
    tryCatch(
      {
        runtime <- get_runtime(x)
        desc <- paste0(desc, "* Average runtime per chain: ",
                       runtime$warmup, " s (warmup) and ",
                       runtime$sampling, " s (sampling)\n")
      },
      error = function(e) {
        cat("* Unable to show runtime info. Reason:\n")
        print(e)
      }
    )
    
    # Convergence info
    tryCatch(
      {
        diag <- x@diagnostics
        Rhat <- diag$Rhat
        imax <- which(Rhat == max(Rhat))
        rhat_str <- round(max(Rhat), 4)
        desc <- paste0(desc, "* Largest R-hat value is ", rhat_str, " (", 
                       paste(rownames(diag)[imax], collapse = ", "), ")\n")
      },
      error = function(e) {
        cat("* Unable to show convergence info. Reason:\n")
        print(e)
      }
    )
    
    # Relevance info
    tryCatch(
      {
        sel <- x@selection
        rel_method <- x@relevances$method
        r3 <- sel$prob
        tr <- sel$threshold
        cat("* Used relevance method = ", rel_method, "\n", sep = "")
        cat("* Used selection threshold = ", tr, "\n", sep = "")
        
        rel <- x@relevances$average
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
    return(desc)
  }
)

#' Character format of an lgpterm
#'
#' @param x an object of class \linkS4class{lgpterm}
#' @param verbose should the format be more verbose
#' @return a string
term_as_character <- function(x, verbose = TRUE) {
  facs <- x@factors
  desc <- sapply(facs, as.character)
  desc <- paste(desc, collapse = ", ")
  return(desc)
}
