context("main function lgp")
library(lgpr)

# Warnings are suppressed because sampling won't converge in such a short time 

test_that("lgp runs", { 
  
  
  expect_identical(suppressWarnings({
    as.vector(lgp(formula = y ~ id + age, 
                  data    = simulate_data(N = 4, 10*c(1,2,3,4,5))$data, 
                  iter    = 100,
                  chains  = 1, 
                  refresh = 1000)@model@stan_dat$D)
  }),
  c(1,1,0,0,0,0)
  )
  
  expect_identical(suppressWarnings({
    as.vector(lgp(formula = y ~ id + age + z, 
                  data    = simulate_data(N          = 4, 
                                          t_data     = 10*c(1,2,3,4,5),
                                          covariates = c(2))$data,
                  iter    = 100,
                  chains  = 1, 
                  refresh = 1000)@model@stan_dat$D)
    }),
    c(1,1,0,0,1,0)
  )
})

test_that("lgp runs without id*age component", {
  
  expect_identical(suppressWarnings({
    as.vector(lgp(formula = y ~ age + z,
                  data    = simulate_data(N          = 4,
                                          t_data     = 10*c(1,2,3,4,5),
                                          covariates = c(2))$data,
                  iter    = 100,
                  chains  = 1,
                  time_variable = "age",
                  id_variable = "id",
                  refresh = 1000)@model@stan_dat$D)
    }),
    c(0,1,0,0,1,0)
  )
})

test_that("lgp runs without age*id component but with shared age", {
  
  expect_identical(suppressWarnings({
    as.vector(lgp(formula = y ~ age + z,
                  data    = simulate_data(N          = 4,
                                          t_data     = 10*c(1,2,3,4,5),
                                          covariates = c(2))$data,
                  iter    = 100,
                  chains  = 1,
                  refresh = 1000)@model@stan_dat$D)
    }),
    c(0,1,0,0,1,0)
  )
  
})


test_that("lgp can sample F", {
  
  expect_identical(suppressWarnings({
    as.vector(lgp(formula = y ~ id + age + diseaseAge + x + z + offset,
                  data    = simulate_data(N          = 4,
                                          t_data     = 10*c(1,2,3,4,5),
                                          covariates = c(0,1,2,3))$data,
                  iter    = 100,
                  chains  = 1,
                  refresh = 1000,
                  offset_vars = c("offset"),
                  sample_F = T)@model@stan_dat$D)
    }),
    c(1,1,1,1,1,1)
  )
  
})

test_that("lgp can used without the vm kernel", {
  
  expect_identical(suppressWarnings({
    as.vector(lgp(formula = y ~ id + age + diseaseAge + group,
                  data    = simulate_data(N          = 4,
                                          t_data     = 10*c(1,2,3,4,5),
                                          covariates = c(0,1,2,4),
                                          dis_fun    = "gp_ns")$data,
                  iter    = 100,
                  chains  = 1,
                  refresh = 1000,
                  offset_vars = c("group"),
                  variance_mask = FALSE)@model@stan_dat$D)
  }),
  c(1,1,1,0,0,1)
  )
  
})
