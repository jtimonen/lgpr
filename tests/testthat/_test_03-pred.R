

r <- relevances(fit)
expect_equal(length(r), 3)
s <- select(fit)
expect_equal(length(s$Component), 3)
t <- seq(0, 40, by = 1)
x_pred <- new_x(dat, t)
p <- pred(fit, x_pred, reduce = NULL, draws = c(1:3), verbose = FALSE)
plt1 <- plot_pred(fit, x_pred, p) # [0,1] scale
plt2 <- plot_f(fit, x_pred, p)
expect_s3_class(plt1, "ggplot")
expect_s3_class(plt2, "ggplot")
