test_that("formula interface gives the same fit as matrix interface", {
  set.seed(20260509L)
  n <- 60
  X <- matrix(rnorm(n * 2), n, 2); colnames(X) <- c("x1", "x2")
  y <- X[, "x1"] + X[, "x2"] + rnorm(n, sd = 0.3)
  df <- data.frame(y = y, X)

  fit_mat <- krls(X = X, y = y, print.level = 0)
  fit_for <- krls(y ~ x1 + x2, data = df, print.level = 0)

  expect_s3_class(fit_for, "krls")
  expect_equal(fit_for$R2, fit_mat$R2, tolerance = 1e-10)
  expect_equal(as.numeric(fit_for$avgderivatives),
               as.numeric(fit_mat$avgderivatives),
               tolerance = 1e-10)
  expect_equal(colnames(fit_for$X), c("x1", "x2"))
})

test_that("formula interface drops intercept and handles factor expansion", {
  set.seed(20260509L)
  n <- 80
  df <- data.frame(
    x = rnorm(n),
    g = factor(sample(c("a", "b", "c"), n, replace = TRUE))
  )
  df$y <- df$x + ifelse(df$g == "b", 0.5, 0) + rnorm(n, sd = 0.3)

  fit <- krls(y ~ x + g, data = df, print.level = 0)
  expect_false("(Intercept)" %in% colnames(fit$X))
  # factor(g) expands to gb, gc dummies (a is reference)
  expect_true("gb" %in% colnames(fit$X))
  expect_true("gc" %in% colnames(fit$X))
})

test_that("formula interface rejects malformed inputs", {
  df <- data.frame(y = rnorm(20), x1 = rnorm(20), x2 = rnorm(20))
  expect_error(krls(y ~ x1 + x2),
               "data.* required")
  expect_error(krls(~ x1 + x2, data = df),
               "two-sided")

  df_na <- df; df_na$y[1] <- NA
  expect_error(krls(y ~ x1 + x2, data = df_na),
               "y contains missing data")
  df_na2 <- df; df_na2$x1[1] <- NA
  expect_error(krls(y ~ x1 + x2, data = df_na2),
               "X contains missing data")
})
