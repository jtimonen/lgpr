.onAttach <- function(...) {
  Lib <- dirname(system.file(package = "lgpr"))
  pkgdesc <- suppressWarnings(utils::packageDescription("lgpr",
    lib.loc = Lib
  ))
  if (length(pkgdesc) > 1) {
    msg <- paste0("This is lgpr (version ", pkgdesc$Version, "). ")
    packageStartupMessage(msg)
  }
}
