.onAttach <- function(...) {
  Lib <- dirname(system.file(package = "lgpr"))
  pkgdesc <- suppressWarnings(utils::packageDescription("lgpr",
    lib.loc = Lib
  ))
  if (length(pkgdesc) > 1) {
    builddate <- gsub(";.*$", "", pkgdesc$Packaged)
    build_str <- paste0("(", builddate, ")")
    packageStartupMessage(paste("This is lgpr, version ",
      pkgdesc$Version,
      " ",
      build_str,
      ".",
      sep = ""
    ))
  }
}
