.onAttach <- function(...) {
  Lib <- dirname(system.file(package = "lgpr"))
  pkgdesc <- suppressWarnings(utils::packageDescription("lgpr",
    lib.loc = Lib
  ))
  if (length(pkgdesc) > 1) {
    builddate <- gsub(";.*$", "", pkgdesc$Packaged)
    packageStartupMessage(paste("This is lgpr, version ",
      pkgdesc$Version, ".",
      sep = ""
    ))
  }
}
