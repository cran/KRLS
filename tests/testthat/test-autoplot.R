test_that("autoplot.krls returns a ggplot when ggplot2 is available", {
  skip_if_not_installed("ggplot2")
  set.seed(20260509L); n <- 60
  X <- matrix(rnorm(n * 2), n, 2); colnames(X) <- c("x1", "x2")
  y <- sin(X[, "x1"]) + 0.3 * X[, "x2"] + rnorm(n, sd = 0.2)
  fit <- krls(X = X, y = y, print.level = 0)

  p <- ggplot2::autoplot(fit)
  expect_s3_class(p, "ggplot")
})

test_that("autoplot errors cleanly when derivatives are absent", {
  skip_if_not_installed("ggplot2")
  set.seed(20260509L); n <- 30
  X <- matrix(rnorm(n * 2), n, 2); colnames(X) <- c("x1", "x2")
  y <- X[, 1] + rnorm(n, sd = 0.3)
  fit_no_d <- krls(X = X, y = y, derivative = FALSE, vcov = FALSE, print.level = 0)
  expect_error(ggplot2::autoplot(fit_no_d), "marginal effects")
})
