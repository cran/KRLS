test_that("krls() returns the expected structure on linear data", {
  d <- make_linear_data()
  fit <- krls(X = d$X, y = d$y, print.level = 0)

  expect_s3_class(fit, "krls")
  expect_named(fit, c("K", "coeffs", "Looe", "fitted", "X", "y",
                      "sigma", "lambda", "R2", "derivatives",
                      "avgderivatives", "var.avgderivatives",
                      "vcov.c", "vcov.fitted", "binaryindicator"),
               ignore.order = TRUE)
  expect_length(fit$fitted, length(d$y))
  expect_gt(fit$R2, 0.9)
  expect_true(fit$lambda > 0)
})

test_that("krls() approximately recovers the truth on a nonlinear example", {
  d <- make_nonlinear_data(N = 200)
  fit <- krls(X = d$X, y = d$y, print.level = 0)

  # Fitted values should be highly correlated with the true mean.
  truth <- d$X[, 1]^3 + 0.5 * d$X[, 2]
  expect_gt(cor(as.numeric(fit$fitted), truth), 0.95)
})

test_that("krls() rejects malformed inputs", {
  d <- make_linear_data()

  # Constant y
  expect_error(krls(X = d$X, y = rep(1, length(d$y)), print.level = 0),
               "constant|does not vary")

  # NA in X
  Xbad <- d$X; Xbad[1, 1] <- NA
  expect_error(krls(X = Xbad, y = d$y, print.level = 0), "missing")

  # NA in y
  ybad <- d$y; ybad[1] <- NA
  expect_error(krls(X = d$X, y = ybad, print.level = 0), "missing")

  # derivative=TRUE with vcov=FALSE
  expect_error(krls(X = d$X, y = d$y, derivative = TRUE, vcov = FALSE,
                    print.level = 0), "vcov")
})

test_that("krls(print.level = 0) is silent", {
  d <- make_linear_data()
  out <- capture.output(fit <- krls(X = d$X, y = d$y, print.level = 0))
  expect_length(out, 0L)
})

test_that("R 4.4+ deprecation warning is gone (regression test for 1.0-0)", {
  d <- make_linear_data()
  expect_no_warning(krls(X = d$X, y = d$y, print.level = 0))
})
