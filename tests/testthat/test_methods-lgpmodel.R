library(lgpr)

# -------------------------------------------------------------------------
dat <- testdata_001

context("Methods for lgpmodel objects")

test_that("param_summary works for models with uncertain effect time", {
  et_info <- list(backwards = TRUE, lower = 3, upper = 23, zero = 10)
  my_prior <- list(
    effect_time_info = et_info,
    wrp = igam(14, 3)
  )
  mod <- create_model(y ~ gp(age) + unc(id) * gp_vm(dis_age),
    dat,
    prior = my_prior
  )
  expect_output(print(param_summary(mod)), "- 10) ~ uniform")
})

test_that("param_summary works for models with heterogeneous effects", {
  mod <- create_model(y ~ gp(age) + het(id) * gp_vm(dis_age), dat,
    prior = list(wrp = igam(14, 3))
  )
  expect_output(print(param_summary(mod)), " ~ beta")
})
