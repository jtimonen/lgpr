library(lgpr)
set.seed(123)

# -------------------------------------------------------------------------

context("Data splitting")
dat <- testdata_001

test_that("data can be split by id", {
  s <- split_by_factor(dat, var_name = "id", test = 1)
  e <- c("train", "test", "i_train", "i_test")
  expect_equal(names(s), e)
  expect_equal(length(s$i_test), 6)
  expect_equal(length(s$i_train), 24 - 6)
})

test_that("data can be split by time point", {
  s <- split_within_factor(dat, var_name = "id", idx_test = 3)
  expect_equal(length(s$i_test), 4)
  expect_equal(length(s$i_train), 24 - 4)
})

test_that("data can be split by time point randomly", {
  s <- split_within_factor_random(dat, var_name = "id", k_test = 2)
  expect_equal(length(s$i_test), 8)
  expect_equal(length(s$i_train), 24 - 8)
})

test_that("data can be split totally randomly", {
  s <- split_random(dat, p_test = 0.5)
  expect_equal(length(s$i_test), 12)
  expect_equal(length(s$i_train), 24 - 12)
})


# -------------------------------------------------------------------------

context("Misc data utilities")

test_that("adding a factor works correctly", {
  country <- c("FIN", "FIN", "EST", "EST")
  expect_error(add_factor(dat, country), "<x> must have names")
  names(country) <- c(1, 3, 2, 4)
  newdat <- add_factor(dat, country)
  expect_true("country" %in% names(newdat))
  expect_error(add_factor(newdat, country), "already contains")
})

test_that("adding a factor crossing works correctly", {
  df_new <- add_factor_crossing(testdata_002, "sex", "group", "xxx")
  expected <- c("Female*Case", "Male*Case", "Female*Control", "Male*Control")
  expect_equal(levels(df_new$xxx), expected)
})

test_that("adding disease age works correctly", {
  t_init <- c(10, 10, 10, 10)
  names(t_init) <- c(1, 2, 3, 4)
  expect_error(add_dis_age(dat, t_init), "already contains")

  newdat <- dat
  newdat$dis_age <- NULL
  newdat <- add_dis_age(newdat, t_init)
  expect_true("dis_age" %in% names(newdat))
})

test_that("adjusted_c_hat works correctly", {
  reason <- "lengths of y and norm_factors must match"
  expect_error(adjusted_c_hat(c(1, 1), norm_factors = c(1, 2, 3)), reason)
  reason <- "<norm_factors> must have only positive values"
  expect_error(adjusted_c_hat(c(1, 1, 0), norm_factors = c(1, 2, 0)), reason)
  reason <- "<y> must have only integer values"
  expect_error(adjusted_c_hat(c(1, 1.2, 0), norm_factors = c(1, 2, 3)), reason)
  c_hat <- adjusted_c_hat(c(1, 1, 1), norm_factors = c(1, 2, 3))
  diff <- abs(c_hat - c(0.0, 0.6931472, 1.0986123))
  expect_lt(max(diff), 1e-6)
})

test_that("computing observed effect times from data works correctly", {
  age <- c(10, 20, 30, 10, 20, 30)
  dage <- c(-10, 0, 10, NaN, NaN, NaN)
  id <- as.factor(c(4, 4, 4, 6, 6, 6))
  df <- data.frame(id, age, dage)
  t0 <- get_teff_obs(df, x_ns = "dage")
  expect_equal(names(t0), c("4", "6"))
  expect_equal(as.numeric(t0), c(20, NaN))
})

test_that("new_x works correctly without x_ns", {
  dat <- testdata_001
  x_new <- new_x(dat, x_values = seq(-1, 1, by = 0.5))
  expect_equal(dim(x_new), c(20, 3))
  expect_equal(colnames(x_new), c("id", "age", "sex"))
})

test_that("new_x works correctly with x_ns", {
  dat <- testdata_001
  x_new <- new_x(dat, x_values = seq(-1, 1, by = 0.5), x_ns = "dis_age")
  expect_equal(dim(x_new), c(20, 4))
  expect_equal(colnames(x_new), c("id", "age", "dis_age", "sex"))
  num_nans <- sum(is.nan(x_new$dis_age))
  expect_equal(num_nans, 10)
})


# -------------------------------------------------------------------------

context("Data plotting")

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
