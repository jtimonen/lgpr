.onAttach <- function(...) {
  Lib <- dirname(system.file(package = "lgpr"))
  pkgdesc <- suppressWarnings(utils::packageDescription("lgpr",
    lib.loc = Lib
  ))
  if (length(pkgdesc) > 1) {
    msg <- paste0(
      "This is lgpr (version ", pkgdesc$Version, ").\n",
      " - Please note the new syntax of the lgp() etc. functions\n",
      " - Old syntax is avalable in version <= 0.33.3."
    )
    packageStartupMessage(msg)
  }
}
