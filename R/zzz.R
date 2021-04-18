.onAttach <- function(...) {
  msg <- create_startup_message()
  packageStartupMessage(msg)
}

# Create package startup message
create_startup_message <- function() {
  v_lgpr <- create_desc("lgpr")
  v_rstan <- create_desc("rstan")
  msg <- paste0(
    "Attached lgpr", v_lgpr, ", using rstan", v_rstan,
    ". Type ?lgpr to get started."
  )
  return(msg)
}

# Create package description
create_desc <- function(pkg_name) {
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
