library(lgpr)
source("helpers/SW.R")

# -------------------------------------------------------------------------
N_ITER <- 22 # should be even
N_CHAINS <- 1

context("Covariate selection")

test_that("selection can be done when sample_f = FALSE", {
  SW({
    fit <- example_fit(iter = N_ITER, chains = N_CHAINS)
  })
  sel_before <- select(fit)
  expect_equal(length(sel_before$Component), 5)

  # After clearing postproc info
  fit <- clear_postproc(fit)
  expect_output(
    {
      expect_message(
        {
          select(fit, verbose = TRUE)
        },
        "No existing postprocessing information stored"
      )
    },
    "Computing analytic function posteriors"
  )
  sel_after <- select(fit, verbose = FALSE)
  expect_equal(sel_before, sel_after)
})

test_that("selection can be done when sample_f = TRUE", {
  SW({
    fit <- example_fit(
      iter = N_ITER,
      chains = N_CHAINS,
      likelihood = "binomial"
    )
  })
  sel_before <- select(fit, threshold = 0.3)
  expect_equal(length(sel_before$Component), 5)

  # After clearing postproc info
  fit <- clear_postproc(fit)
  sel_after <- select(fit, verbose = FALSE, threshold = 0.3)
  expect_equal(sel_before, sel_after)
})

test_that("other select() functions can be used", {
  SW({
    fit <- example_fit(iter = N_ITER, chains = N_CHAINS)
  })
  sel <- select_freq(fit)
  expect_equal(length(sel$Component), 5)

  sel <- select.integrate(fit, verbose = FALSE)
  expect_equal(dim(sel$selected), c(101, 5))

  sel <- select_freq.integrate(fit, verbose = FALSE)
  expect_equal(dim(sel$freq), c(101, 5))
})
