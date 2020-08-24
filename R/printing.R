#' Is an object printable
#'
#' @param object an object
is_printable <- function(object) {
  d <- dim(object)
  if (is.null(d)) {
    return(TRUE)
  } else if (prod(d) == 0) {
    return(FALSE)
  }
  TRUE
}

#' Print a list in a more compact format
#'
#' @param input a named list
print_list <- function(input) {
  nam <- names(input)
  printed <- c()
  skipped <- c()
  for (name in nam) {
    f <- input[[name]]
    if (is_printable(f)) {
      printed <- c(printed, name)
    } else {
      skipped <- c(skipped, name)
    }
  }

  print(input[printed])
  str <- paste(skipped, collapse = ", ")
  msg <- paste0(
    "Did not print fields with at least one zero dimension:\n    ",
    str, "\n"
  )
  cat(msg)
  invisible(input)
}

#' Print the Stan input of an lgpmodel
#'
#' @param model an object of class \linkS4class{lgpmodel}
print_stan_input <- function(model) {
  print_list(model@stan_input)
  invisible(model)
}
