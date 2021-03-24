# Source all test helper functions
parent_dir <- 'tests/testthat/helpers/'
files <- dir(parent_dir)
for (fn in files) {
  fn <- file.path(parent_dir, fn)
  msg <- paste0(" - ", fn, "\n")
  cat(msg)
  source(fn)
}
