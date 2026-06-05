.fit_toy <- function() {
  set.seed(20260509L)
  n <- 60
  X <- matrix(rnorm(n * 2), n, 2); colnames(X) <- c("x1", "x2")
  y <- sin(X[, "x1"]) + 0.3 * X[, "x2"] + rnorm(n, sd = 0.2)
  krls(X = X, y = y, print.level = 0)
}

test_that("tidy.krls returns one row per predictor with the documented columns", {
  fit <- .fit_toy()
  out <- tidy.krls(fit)
  expect_named(out, c("term", "estimate", "std.error", "statistic", "p.value",
                      "conf.low", "conf.high",
                      "q25", "median", "q75", "binary"))
  expect_equal(nrow(out), 2L)
  # AMEs should match the underlying field
  expect_equal(out$estimate, as.numeric(fit$avgderivatives), tolerance = 1e-12)
  # CI brackets the estimate
  expect_true(all(out$conf.low <= out$estimate))
  expect_true(all(out$estimate <= out$conf.high))
})

test_that("glance.krls returns the documented summary columns", {
  fit <- .fit_toy()
  out <- glance.krls(fit)
  expect_equal(nrow(out), 1L)
  expect_named(out, c("nobs", "n_predictors", "r.squared", "loo_mse",
                      "lambda", "sigma", "eff_df",
                      "approx", "nystrom_m", "inference"))
  # Exact-path fits report approx = "none" with a NA nystrom_m.
  expect_equal(out$approx, "none")
  expect_true(is.na(out$nystrom_m))
  expect_equal(out$nobs, 60L)
  expect_equal(out$n_predictors, 2L)
  expect_true(out$r.squared > 0 && out$r.squared <= 1)
  expect_true(is.finite(out$eff_df))
})

test_that("augment.krls joins .fitted/.resid/derivative columns to data", {
  fit <- .fit_toy()
  set.seed(20260509L); n <- 60
  X <- matrix(rnorm(n * 2), n, 2); colnames(X) <- c("x1", "x2")
  y <- sin(X[, "x1"]) + 0.3 * X[, "x2"] + rnorm(n, sd = 0.2)
  df <- data.frame(y = y, X)

  out <- augment.krls(fit, data = df)
  expect_true(all(c(".fitted", ".resid", ".dy_d_x1", ".dy_d_x2") %in% names(out)))
  expect_equal(nrow(out), n)
  # Without supplied data, augment uses the X matrix from the fit
  out2 <- augment.krls(fit)
  expect_true(all(c(".fitted", ".resid", ".dy_d_x1", ".dy_d_x2") %in% names(out2)))
})

test_that("tidy/augment require derivative = TRUE", {
  set.seed(20260509L); n <- 30
  X <- matrix(rnorm(n * 2), n, 2); colnames(X) <- c("x1", "x2")
  y <- X[, 1] + rnorm(n, sd = 0.3)
  fit_no_d <- krls(X = X, y = y, derivative = FALSE, vcov = FALSE, print.level = 0)
  expect_error(tidy.krls(fit_no_d), "marginal effects")
})
