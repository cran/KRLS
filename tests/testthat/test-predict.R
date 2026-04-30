test_that("predict() reproduces fitted values when newdata = X", {
  d <- make_linear_data()
  fit <- krls(X = d$X, y = d$y, derivative = FALSE, vcov = FALSE,
              print.level = 0)
  pr <- predict(fit, newdata = d$X, se.fit = FALSE)
  expect_equal(as.numeric(pr$fit), as.numeric(fit$fitted), tolerance = 1e-8)
})

test_that("predict() returns a fit + se.fit when se.fit = TRUE", {
  d <- make_linear_data()
  fit <- krls(X = d$X, y = d$y, print.level = 0)
  pr <- predict(fit, newdata = d$X, se.fit = TRUE)
  expect_named(pr, c("fit", "se.fit", "vcov.fit", "newdata", "newdataK"))
  expect_equal(length(pr$fit), nrow(d$X))
  expect_equal(length(pr$se.fit), nrow(d$X))
  expect_true(all(pr$se.fit >= 0))
})

test_that("predict() rejects ncol mismatch", {
  d <- make_linear_data()
  fit <- krls(X = d$X, y = d$y, derivative = FALSE, vcov = FALSE,
              print.level = 0)
  expect_error(predict(fit, newdata = d$X[, 1, drop = FALSE]),
               "ncol")
})

test_that("predict() rejects non-krls input (clean error, not warning)", {
  expect_error(predict.krls(list(a = 1)),
               "krls")
})

test_that("predict() rejects se.fit = TRUE on a fit without vcov", {
  d <- make_linear_data()
  fit <- krls(X = d$X, y = d$y, derivative = FALSE, vcov = FALSE,
              print.level = 0)
  expect_error(predict(fit, newdata = d$X, se.fit = TRUE),
               "vcov")
})
