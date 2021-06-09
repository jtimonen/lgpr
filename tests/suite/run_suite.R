#!/usr/bin/env Rscript

# A larger test suite for lgpr
# -------------------------------------
#  - tests correctness of inference
#  - measures runtime
#
# Rscript run_suite.R  <num_iter> <suite_path> <idx_test> <verbose>
#
# Example uses:
# - Rscript run_suite.R
# - Rscript run_suite.R 500
# - Rscript run_suite.R 500 tests/suite 3
library(rmarkdown)
library(lgpr)
library(rstan)
library(lme4) # for sleepstudy data
library(nlme) # for Orthodont data

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  NUM_ITER <- 2000
  msg <- paste0("Defaulting to NUM_ITER=", NUM_ITER, "\n")
  cat(msg)
} else {
  NUM_ITER <- as.numeric(args[1])
}
if (length(args) == 2) {
  suite_path <- args[2]
} else {
  suite_path <- file.path(".")
}

# Which tests to run (0 = all)
if (length(args) == 3) {
  IDX <- as.numeric(args[3])
} else {
  IDX <- 0
}


# Which tests to run (0 = all)
if (length(args) == 4) {
  verbose <- as.logical(args[4])
} else {
  verbose <- FALSE
}

# Common settings for all tests
NUM_CHAINS <- 4
NUM_CORES <- 4
REFRESH <- 0
STAN_SEED <- 123
DRAW_INDS <- round(seq(1, NUM_ITER * NUM_CHAINS / 2, length.out = 10))

# Set paths
models_path <- file.path(suite_path, "models")
Rmd_path <- file.path(suite_path, "Rmd")
out_path <- file.path(suite_path, "out")
rds_path <- file.path(out_path, "rds")
dir.create(out_path)
dir.create(rds_path)
files <- dir(models_path)
if (IDX != 0) {
  files <- files[IDX]
}

# Source helper files
source(file.path(suite_path, "common.R"))

# Run the test suite
INFO <- c()
HR <- "-----------------------------------------------------------------------"
HR <- paste0("\u001b[1m\u001b[36m", HR, "\u001b[0m\n")
cat(HR)
for (f in files) {

  # Setup
  r_file <- file.path(models_path, f)
  base_name <- strsplit(f, "[.]")[[1]][1]
  rds_file <- file.path(rds_path, paste0(base_name, ".rds"))
  Rmd_file <- file.path(Rmd_path, paste0(base_name, ".Rmd"))
  html_file <- paste0(base_name, ".html")
  MSG <- paste0("\u001b[1m\u001b[36mRunning: ", base_name, "\u001b[0m\n")
  cat(MSG)
  start_time <- Sys.time()
  source(r_file)

  # Run model fitting
  res_fit <- run_example(
    verbose,
    iter = NUM_ITER, chains = NUM_CHAINS, cores = NUM_CORES,
    refresh = REFRESH, seed = STAN_SEED
  )
  fit <- res_fit$fit
  t_fit <- res_fit$time

  # Save the fit object
  saveRDS(fit, file = rds_file)
  size_disk <- file_size_kb(rds_file)

  # Time pred
  t_pred <- run_pred(fit, verbose)

  # Run other post-fitting tasks and knit result
  render_start_time <- Sys.time()
  rmarkdown::render(
    input = Rmd_file, output_file = html_file,
    output_dir = out_path,
    quiet = !verbose
  )
  t_post <- as.double(Sys.time() - render_start_time, units = "secs")

  # Store info
  t_total <- as.double(Sys.time() - start_time, units = "secs")
  info_f <- get_info(fit, base_name, t_fit, t_pred, t_post, t_total, size_disk)
  INFO <- rbind(INFO, info_f)
}
MSG <- paste0("\u001b[1m\u001b[36mFinished. \u001b[0m\n")
cat(MSG)
cat(HR)
cat("\n")

INFO <- round_results(INFO, 2L, 3L)
print(INFO)
