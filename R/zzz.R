.onAttach <- function(...) {
  get_desc <- function(pkg_name) {
    Lib <- dirname(system.file(package = pkg_name))
    pkgdesc <- suppressWarnings(
      utils::packageDescription(pkg_name, lib.loc = Lib)
    )
    if (length(pkgdesc) > 1) {
      out <- paste0(" ", pkgdesc$Version)
    } else {
      out <- ""
    }
    return(out)
  }

  v1 <- get_desc("lgpr")
  v2 <- get_desc("rstan")
  msg <- paste0("Attached lgpr", v1, ", using rstan", v2, ".\n")
  packageStartupMessage(msg)
}
