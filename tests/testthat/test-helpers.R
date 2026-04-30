test_that("gausskernel() returns a symmetric matrix with 1's on the diagonal", {
  set.seed(20260429L)
  X <- matrix(rnorm(60), 20, 3)
  K <- gausskernel(X = X, sigma = 3)
  expect_equal(dim(K), c(20L, 20L))
  expect_equal(unname(diag(K)), rep(1, 20), tolerance = 1e-12)
  expect_equal(K, t(K), tolerance = 1e-12)
})

test_that("solveforc() and looloss() work on a small example", {
  set.seed(20260429L)
  N  <- 50
  X  <- cbind(rnorm(N), rnorm(N))
  y  <- as.matrix(X[, 1] + rnorm(N, 0, 0.1))
  K  <- gausskernel(X, sigma = ncol(X))
  E  <- eigen(K, symmetric = TRUE)
  out <- solveforc(y = y, Eigenobject = E, lambda = 0.1)
  expect_equal(length(out$coeffs), N)
  expect_equal(length(out$Le), 1L)
  expect_equal(as.numeric(looloss(y = y, Eigenobject = E, lambda = 0.1)),
               as.numeric(out$Le), tolerance = 1e-12)
})
