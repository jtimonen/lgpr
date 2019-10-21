.onLoad <- function(libname, pkgname) { # nocov start
  sf <- rstan::expose_stan_functions(
    stanmodels$lgp, 
    cacheDir = file.path('src', 'stan_functions'),
                        cleanupCacheDir = TRUE)
} # nocov end

.onAttach <- function(...) {
  Lib <- dirname(system.file(package = "lgpr"))
  pkgdesc <- suppressWarnings(utils::packageDescription("lgpr",
                                                        lib.loc = Lib))
  if (length(pkgdesc) > 1) {
    builddate <- gsub(';.*$', '', pkgdesc$Packaged)
    packageStartupMessage(paste("Hello, this is lgpr (version ", 
                                pkgdesc$Version,").", sep = ""))
  }
}