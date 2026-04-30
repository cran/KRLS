test_that("summary.krls() produces a coefficients table", {
  d <- make_linear_data()
  fit <- krls(X = d$X, y = d$y, print.level = 0)
  out <- capture.output(s <- summary(fit))

  expect_s3_class(s, "summary.krls")
  expect_true(any(grepl("R2", out)))
  expect_true(any(grepl("Average Marginal Effects", out)))
  expect_equal(colnames(s$coefficients),
               c("Est", "Std. Error", "t value", "Pr(>|t|)"))
  expect_equal(nrow(s$coefficients), ncol(d$X))
})

test_that("summary.krls() returns invisible NULL when no derivatives stored", {
  d <- make_linear_data()
  fit <- krls(X = d$X, y = d$y, derivative = FALSE, vcov = FALSE,
              print.level = 0)
  out <- capture.output(s <- summary(fit))
  expect_null(s)
  expect_true(any(grepl("derivative = TRUE", out)))
})

test_that("summary.krls() rejects non-krls input", {
  expect_error(summary.krls(list(a = 1)),
               "krls")
})
