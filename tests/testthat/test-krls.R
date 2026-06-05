test_that("krls() returns the expected structure on linear data", {
  d <- make_linear_data()
  fit <- krls(X = d$X, y = d$y, print.level = 0)

  expect_s3_class(fit, "krls")
  expect_named(fit, c("K", "coeffs", "Looe", "fitted", "X", "X_proc",
                      "prep", "y",
                      "sigma", "lambda", "R2", "derivatives",
                      "avgderivatives", "var.avgderivatives",
                      "vcov.c", "vcov.fitted", "binaryindicator",
                      "lambda_method"),
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
  # cat_columns = integer(0) acknowledges the binary x2 column the
  # fixture treats as continuous, silencing the v1.7 categorical
  # nudge so this test still surfaces any R 4.4-deprecation warning
  # if one comes back.
  expect_no_warning(krls(X = d$X, y = d$y, cat_columns = integer(0),
                         print.level = 0))
})

test_that("user-supplied tol is honored by lambdasearch (regression for 1.3-0)", {
  d <- make_nonlinear_data()
  # A very loose tol terminates the golden-section search early at a
  # different lambda than the default tight tol; if krls() drops tol
  # on the floor, both fits get the same lambda.
  fit_loose <- krls(X = d$X, y = d$y, tol = 1e9, print.level = 0)
  fit_tight <- krls(X = d$X, y = d$y, print.level = 0)
  expect_false(isTRUE(all.equal(fit_loose$lambda, fit_tight$lambda)))
})

test_that("lambda_method = 'gcv' returns a plausible lambda (exact path)", {
  d <- make_nonlinear_data(N = 100)
  fit_loo <- krls(d$X, d$y, lambda_method = "loo", print.level = 0)
  fit_gcv <- krls(d$X, d$y, lambda_method = "gcv", print.level = 0)

  expect_equal(fit_loo$lambda_method, "loo")
  expect_equal(fit_gcv$lambda_method, "gcv")
  expect_true(fit_gcv$lambda > 0 && is.finite(fit_gcv$lambda))
  # GCV fit should still recover the truth.
  expect_gt(fit_gcv$R2, 0.85)
  # (We don't pin the GCV lambda to LOO's: the two objectives can
  # legitimately differ by orders of magnitude depending on the data,
  # particularly when the LOO loss is flat across a wide lambda range.)
})

test_that("default lambda_method is 'loo' (backward compatibility)", {
  d <- make_nonlinear_data(N = 80)
  fit_default <- krls(d$X, d$y, print.level = 0)
  fit_loo     <- krls(d$X, d$y, lambda_method = "loo", print.level = 0)
  expect_equal(fit_default$lambda_method, "loo")
  expect_equal(fit_default$lambda, fit_loo$lambda, tolerance = 1e-12)
})

test_that("aggressive eigtrunc that keeps a single eigenvector works", {
  d <- make_nonlinear_data()
  # eigtrunc = 0.99999 keeps only the largest eigenvector; the
  # single-column subset must not drop to a numeric vector or
  # solveforc() / vcovmatc construction fails inside multdiag().
  expect_no_error(
    krls(X = d$X, y = d$y, eigtrunc = 0.99999, print.level = 0)
  )
})
