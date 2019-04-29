# On load
.onLoad <- function(libname, pkgname) { # nocov start
  if(length(stanmodels) > 0){
    modules <- paste0("stan_fit4", names(stanmodels), "_mod")
    for (m in modules) {
      loadModule(m, what = TRUE)
    }
  }
} # nocov end


# On attach
.onAttach <- function(libname, pkgname) {
  lgprLib <- dirname(system.file(package = "lgpr"))
  descr   <- suppressWarnings(utils::packageDescription("lgpr",
                                                        lib.loc = lgprLib))
  if (length(descr) > 1) {
    builddate <- gsub(';.*$', '', descr$Packaged)
    if(length(builddate)<4){
      packageStartupMessage(paste("This is lgpr (version ", descr$Version,")", sep = ""))
    }else{
      packageStartupMessage(paste("This is lgpr (version ", descr$Version,
                                  ", packaged: ", builddate, ")", sep = ""))
    }
  }
}
