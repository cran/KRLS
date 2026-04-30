## Shared fixtures used across the test suite.

# Reproducible toy with a linear additive truth.
make_linear_data <- function(N = 80, seed = 20260429L) {
  set.seed(seed)
  x1 <- rnorm(N)
  x2 <- rbinom(N, 1, 0.3)
  y  <- x1 + 0.5 * x2 + rnorm(N, 0, 0.2)
  list(X = cbind(x1, x2), y = y)
}

# Reproducible toy with x1^3 nonlinearity.
make_nonlinear_data <- function(N = 80, seed = 20260429L) {
  set.seed(seed)
  x1 <- rnorm(N)
  x2 <- rbinom(N, 1, 0.3)
  y  <- x1^3 + 0.5 * x2 + rnorm(N, 0, 0.2)
  list(X = cbind(x1, x2), y = y)
}
