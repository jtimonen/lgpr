library(lgpr)

# -------------------------------------------------------------------------

context("Plotting")

test_that("color_set works", {
  col1 <- colorset("blue")
  col2 <- colorset("red", "dark")
  expect_equal(nchar(col1), 7)
  expect_equal(nchar(col2), 7)
  expect_error(colorset("definitely_not_a_color_name"))
  expect_error(colorset("green", "asdf"), "Invalid color")
})

test_that("plot_color_palette works", {
  p1 <- plot_color_palette(1)
  p2 <- plot_color_palette(5)
  p3 <- plot_color_palette(6)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_error(color_palette_plot(7))
})

# Create test data
age <- c(10, 20, 30, 10, 20, 30)
id <- as.factor(c(1, 1, 1, 2, 2, 2))
sex <- as.factor(c("Male", "Male", "Male", "Female", "Female", "Female"))
dis_age <- c(12 - age[1:3], NA, NA, NA)
y <- c(1, 2, 4, 5, 0, 1)
dat <- data.frame(id, age, dis_age, sex, y)

test_that("plot_data can plot data frames", {
  p1 <- plot_data(dat)
  p2 <- plot_data(dat, highlight = 1, group_by = "id")
  p3 <- plot_data(dat, highlight = "Male", group_by = "sex", facet_by = "sex")
  p4 <- plot_data(dat, color_by = "sex")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  reason <- "Invalid <highlight> argument 11! The possible values of id are"
  expect_error(plot_data(dat, highlight = 11), reason)
})

test_that("simulated data can be plotted", {
  dat <- simulate_data(
    N = 4,
    t_data = seq(6, 36, by = 6),
    covariates = c(1, 2)
  )
  p <- sim_plot(dat, i_test = c(1, 2, 3), ncol = 4)
  expect_s3_class(p, "ggplot")
})

test_that("simulated data with disease effect can be plotted", {
  dat <- simulate_data(
    N = 4,
    t_data = seq(6, 36, by = 6),
    covariates = c(0, 2),
    t_observed = "after_1"
  )
  p <- sim_plot(dat, i_test = c(1, 2, 3), ncol = 4)
  expect_s3_class(p, "ggplot")
})


test_that("model fit can be visualized with data", {
  sim <- simulate_data(
    N = 4,
    t_data = seq(6, 36, by = 6),
    covariates = c(0, 2),
    t_observed = "after_1"
  )
  data <- sim@data
  suppressWarnings({
    fit <- lgp(y ~ gp(age) + zerosum(z) * gp(age),
      data = data, chains = 1, iter = 500, refresh = 0
    )
    p <- plot_fit(fit, data)
    expect_s3_class(p, "ggplot")
  })
})
