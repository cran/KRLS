## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 7,
  fig.height = 4.5,
  dpi = 96
)
set.seed(2026)

## ----setup--------------------------------------------------------------------
library(KRLS)

## ----dgp----------------------------------------------------------------------
make_data <- function(n) {
  X <- cbind(rnorm(n), rnorm(n), rnorm(n), rbinom(n, 1, 0.4))
  colnames(X) <- c("x1", "x2", "x3", "x4")
  f <- sin(X[, 1]) + 0.3 * X[, 2] + 0.05 * X[, 3]^2 + 0.5 * X[, 4]
  list(X = X, y = f + rnorm(n, sd = 0.3), truth = f)
}

## ----bench, cache = FALSE-----------------------------------------------------
sizes <- c(500, 1000)
res <- vector("list", length(sizes))

for (i in seq_along(sizes)) {
  n <- sizes[i]
  tr <- make_data(n)
  te <- make_data(200)
  # Use the package default landmark count, min(500, n).
  m  <- min(500L, n)

  t_e <- system.time({
    fit_e <- krls(tr$X, tr$y, approx = "none",
                  derivative = FALSE, vcov = FALSE,
                  print.level = 0)
  })["elapsed"]
  t_n <- system.time({
    fit_n <- krls(tr$X, tr$y, approx = "nystrom",
                  derivative = FALSE, print.level = 0)
  })["elapsed"]

  pred_e <- predict(fit_e, newdata = te$X)$fit
  pred_n <- predict(fit_n, newdata = te$X)$fit

  res[[i]] <- data.frame(
    n          = n,
    m          = m,
    time_exact = round(as.numeric(t_e), 2),
    time_nys   = round(as.numeric(t_n), 3),
    speedup    = round(as.numeric(t_e) / as.numeric(t_n), 1),
    rmse_exact = round(sqrt(mean((pred_e - te$truth)^2)), 4),
    rmse_nys   = round(sqrt(mean((pred_n - te$truth)^2)), 4)
  )
}
do.call(rbind, res)

## ----big----------------------------------------------------------------------
n  <- 2000
tr <- make_data(n)
te <- make_data(200)
# Default is min(500, n) -- here 500 landmarks at n = 2000.
m  <- min(500L, n)

t_n <- system.time({
  fit_big <- krls(tr$X, tr$y, approx = "nystrom",
                  derivative = FALSE, print.level = 0)
})["elapsed"]

pred_big <- predict(fit_big, newdata = te$X)$fit
sprintf("n = %d, m = %d, time = %.3fs, R2 = %.3f, RMSE = %.3f",
        n, m, t_n, fit_big$R2,
        sqrt(mean((pred_big - te$truth)^2)))

## ----reuse--------------------------------------------------------------------
tr <- make_data(400)
fit <- krls(tr$X, tr$y, approx = "nystrom", nystrom_m = 30,
            derivative = FALSE, print.level = 0)

Z <- get_landmarks(fit)              # original X-scale matrix
y_perturbed <- tr$y + rnorm(400, sd = 0.05)

fit2 <- krls(tr$X, y_perturbed, approx = "nystrom",
             landmarks = Z, derivative = FALSE, print.level = 0)

# Same landmarks under the hood
identical(unname(fit$landmarks), unname(fit2$landmarks))

## ----gcv----------------------------------------------------------------------
tr <- make_data(500)
fit_loo <- krls(tr$X, tr$y, approx = "nystrom", nystrom_m = 25,
                derivative = FALSE, lambda_method = "loo",
                print.level = 0)
fit_gcv <- krls(tr$X, tr$y, approx = "nystrom", nystrom_m = 25,
                derivative = FALSE, lambda_method = "gcv",
                print.level = 0)
data.frame(
  method = c("loo", "gcv"),
  lambda = c(fit_loo$lambda, fit_gcv$lambda),
  R2     = c(fit_loo$R2, fit_gcv$R2)
)

