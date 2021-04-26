# A larger test suite
#  - tests correctness of inference
#  - measures runtime
library(rmarkdown)
library(lgpr)
library(rstan)

# Common settings for all tests
NUM_ITER <- 2000
NUM_CHAINS <- 4
NUM_CORES <- 4
REFRESH <- 0
STAN_SEED <- 123
DRAW_INDS <- round(seq(1, NUM_ITER * NUM_CHAINS / 2, length.out = 100))
verbose <- FALSE

# Set paths
suite_path <- file.path("tests", "suite")
models_path <- file.path(suite_path, "models")
Rmd_path <- file.path(suite_path, "Rmd")
out_path <- file.path(suite_path, "out")
rds_path <- file.path(out_path, "rds")
dir.create(out_path)
dir.create(rds_path)
files <- dir(models_path)

# Source helper files
source(file.path(suite_path, "common.R"))

# Run the test suite
INFO <- c()
for (f in files) {

  # Setup
  r_file <- file.path(models_path, f)
  base_name <- strsplit(f, "[.]")[[1]][1]
  rds_file <- file.path(rds_path, paste0(base_name, ".rds"))
  Rmd_file <- file.path(Rmd_path, paste0(base_name, ".Rmd"))
  html_file <- paste0(base_name, ".html")
  cat("Running:", base_name, "\n")
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
  size_disk <- file_size_mb(rds_file)

  # Time pred
  t_pred <- run_pred(fit, verbose)

  # Run other post-fitting tasks and knit result
  render_start_time <- Sys.time()
  rmarkdown::render(
    input = Rmd_file, output_file = html_file,
    output_dir = out_path,
    quiet = !verbose
  )
  t_post <- as.numeric(Sys.time() - render_start_time)

  # Store info
  t_total <- as.numeric(Sys.time() - start_time)
  info_f <- get_info(fit, base_name, t_fit, t_pred, t_post, t_total, size_disk)
  INFO <- rbind(INFO, info_f)
}
cat("Finished.\n\n")

INFO <- round_results(INFO, 2L, 3L)
print(INFO)
