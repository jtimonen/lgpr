# A larger test suite
#  - tests correctness of inference
#  - measures runtime
library(bench)
library(lgpr)
NUM_CORES <- 4

suite_path <- file.path("tests/suite")
files <- dir(suite_path)
for (f in  files) {
  fp <- file.path(suite_path, f)
  cat("FILE:", fp, "\n")
  source(fp)
  
  fit <- run_test(cores = NUM_CORES)
}
