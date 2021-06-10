library(lgpr)

# -------------------------------------------------------------------------

context("Tests for entire workflow (instead of single function)")

test_that("everything works when the dependent variable is not named 'y'", {
  dat <- testdata_001
  colnames(dat)[6] <- "depvar"
  suppressWarnings({
    fit <- lgp(depvar ~ age + sex,
      data = dat, iter = 100, chains = 1,
      refresh = 0, quiet = TRUE
    )
  })
  t <- seq(1, 50, by = 2)
  x_pred <- new_x(dat, t, x = "age")
  p <- pred(fit, x_pred, reduce = mean, verbose = FALSE)
  expect_equal(get_y_name(fit), "depvar")
  h1 <- plot_components(fit, pred = p)
  h2 <- plot_f(fit, pred = p)
  h3 <- plot_pred(fit, pred = p, t_name = "age")
  expect_equal(length(h1), 3)
  expect_s3_class(h2, "ggplot")
  expect_s3_class(h3, "ggplot")
})



test_that("everything works when id variable is not named 'id'", {
  dat <- testdata_001
  dat$SUBJECT <- dat$id
  dat$id <- NULL
  suppressWarnings({
    fit <- lgp(y ~ gp(age) + gp_vm(dis_age),
      data = dat, iter = 100, chains = 1,
      refresh = 0, quiet = TRUE
    )
  })

  # get_teff_obs works correctly
  expect_error(get_teff_obs(dat, x_ns = "dis_age"), "'id' not found in <data>")
  teff <- get_teff_obs(dat, x_ns = "dis_age", group_by = "SUBJECT")
  expect_equal(as.numeric(teff), c(21, 21, NaN, NaN))

  # need to give group_by = "SUBJECT" to every function now
  t <- seq(1, 50, by = 2)
  x_pred <- new_x(dat, t, x = "age", group_by = "SUBJECT", x_ns = "dis_age")
  p <- pred(fit, x_pred, reduce = mean, verbose = FALSE)
  h1 <- plot_components(fit, pred = p, group_by = "SUBJECT")
  h2 <- plot_f(fit, pred = p, group_by = "SUBJECT")
  h3 <- plot_pred(fit, pred = p, t_name = "age", group_by = "SUBJECT")
  expect_equal(length(h1), 3)
  expect_s3_class(h2, "ggplot")
  expect_s3_class(h3, "ggplot")
})
