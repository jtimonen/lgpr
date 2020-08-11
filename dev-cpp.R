# This file is only for development purposes and not a part of the package
# itself. Run this in the package root directory to build the stan functions
# into C++ code.

# Set paths
STAN_HOME <- file.path("inst", "stan")
FILES <- c(
  "chunks/functions-utils.stan",
  "chunks/functions-kernels_base.stan",
  "chunks/functions-kernels_single.stan",
  "chunks/functions-kernels_many.stan",
  "chunks/functions-posterior.stan",
  "chunks/functions-prior.stan"
)

# Create stan model containing only a functions block with all the functions
two_spaces <- "  " 
f_list <- lapply(file.path(STAN_HOME, FILES), FUN = readLines)
functions <- paste(unlist(f_list), collapse = paste0("\n", two_spaces))
functions <- paste0(two_spaces, functions)
model_code <- paste(c("functions {", functions, "}"), collapse = "\n")
model_code <- paste0(model_code, "\n")

# Write Stan code to file
fn <- file.path(STAN_HOME, 'all_functions.stan')
cat(model_code, file = fn)

# Transpile to C++ code and delete the Stan file
cpp_code <- rstan::expose_stan_functions(fn, verbose = TRUE, dryRun = TRUE)
file.remove(fn)

# Write the C++ code to src/stanFunctions.cpp
cat(cpp_code, file = file.path('src', 'stanFunctions.cpp'))

# Update R/RcppExports.R and src/RcppExports.cpp
Rcpp::compileAttributes(pkgdir = '.', verbose = TRUE)

# Add things
add1 <- paste0('RcppExport SEXP _rcpp_module_boot_stan_fit4lgp_mod(); \n',
               'RcppExport SEXP _rcpp_module_boot_stan_fit4lgp_lc_mod(); \n')

# at the the start of static const R_CallMethodDef CallEntries[] function
add2 <- paste0('{"_rcpp_module_boot_stan_fit4lgp_mod", ', 
               '(DL_FUNC) &_rcpp_module_boot_stan_fit4lgp_mod, 3} \n',
               '{"_rcpp_module_boot_stan_fit4lgp_lc_mod", ',
               '(DL_FUNC) &_rcpp_module_boot_stan_fit4lgp_lc_mod, 3}\n')

cat("\n === Add following to RcppExports.cpp, ")
cat("at the start of R_CallMethodDef CallEntries[] == \n")
cat(paste0(add1, add2))

