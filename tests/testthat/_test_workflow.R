library(lgpr)

# -------------------------------------------------------------------------

context("Tests for entire workflow (instead of single function)")

test_that("everything works when the dependent variable is not named 'y'", {
  dat <- testdata_001
  colnames(dat)[6] <- "depvar"
  suppressWarnings({
    fit <- lgp(depvar ~ age + sex,
      data = dat, iter = 100, chains = 1,
      refresh = 0
    )
  })
  t <- seq(1, 50, by = 2)
  x_pred <- new_x(dat, t, x = "age")
  p <- pred(fit, x_pred, reduce = mean, verbose = FALSE)
  expect_equal(get_y_name(fit), "depvar")
  h1 <- plot_components(fit, x = x_pred, pred = p)
  h2 <- plot_f(fit, x = x_pred, pred = p)
  h3 <- plot_pred(fit, x = x_pred, pred = p, t_name = "age")
  expect_equal(length(h1), 3)
  expect_s3_class(h2, "ggplot")
  expect_s3_class(h3, "ggplot")
})
