library(lgpr)

# -------------------------------------------------------------------------

context("Plotting data")


# Create test data
age <- c(10, 20, 30, 10, 20, 30)
id <- as.factor(c(1, 1, 1, 2, 2, 2))
sex <- as.factor(c("Male", "Male", "Male", "Female", "Female", "Female"))
dis_age <- c(12 - age[1:3], NA, NA, NA)
y <- c(1, 2, 4, 5, 0, 1)
dat <- data.frame(id, age, dis_age, sex, y)

test_that("plot_data works", {
  p1 <- plot_data(dat)
  p2 <- plot_data(dat, highlight = 1, group_by = "id")
  p3 <- plot_data(dat, highlight = "Male", group_by = "sex", facet_by = "sex")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})
