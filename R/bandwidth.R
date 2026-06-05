## Bandwidth selection by maximizing off-diagonal var(K).
##
## Background: the Gaussian-kernel bandwidth sigma controls how fast
## kernel similarity decays with squared Euclidean distance:
##   K_ij = exp(-||x_i - x_j||^2 / sigma).
## At sigma -> 0 every off-diagonal K_ij -> 0 (var collapses to 0);
## at sigma -> Inf every off-diagonal K_ij -> 1 (var collapses to 0).
## Somewhere in between sits the sigma that makes the columns of K
## most informative -- the off-diagonals carry the most variation.
## This is the choice used in kbal::b_maxvarK and gpss:::getb_maxvar
## and the one Chad now recommends as the default for krls().
##
## The exact-path version evaluates Var(K[off-diagonal]) over the full
## n x n kernel. The Nystrom version evaluates Var(C) over the n x m
## cross-kernel between observations and landmarks -- following Chad's
## suggestion in the email thread that "a Nystrom version would just
## compute var on the C matrix instead of full K".

## .var_offdiag_K(): variance of the off-diagonal entries of a square
## symmetric matrix K, without copying. Off-diagonal because the
## diagonal is structurally 1 for Gaussian kernels and would otherwise
## dominate the variance estimate.
##
## Identity used: Var(off) = E[K_off^2] - (E[K_off])^2 where the means
## are taken over the n*(n-1) off-diagonal entries. We get these from
## sum(K^2) - n and sum(K) - n (subtracting the n-on-diagonal 1's).
.var_offdiag_K <- function(K) {
  n <- nrow(K)
  if (n < 2L) return(0)
  n_off <- n * (n - 1L)
  s2 <- sum(K * K) - n        # sum of off-diagonal squared entries
  s1 <- sum(K) - n            # sum of off-diagonal entries
  mean_off <- s1 / n_off
  s2 / n_off - mean_off * mean_off
}

## .var_C(): variance of the entries of an n x m cross-kernel matrix.
## Used in the Nystrom path where we never materialize the full K. No
## diagonal correction is needed because C is rectangular.
.var_C <- function(Cmat) {
  N <- length(Cmat)
  if (N < 2L) return(0)
  m <- mean(Cmat)
  sum((Cmat - m)^2) / N
}

## b_maxvarK(): exact-path bandwidth selector. X_proc is the
## already-preprocessed numeric matrix (continuous columns
## standardized, categoricals sqrt(0.5)-one-hot). Returns the sigma
## that maximizes Var(K_off) over the search interval.
##
## Search bounds default to (1e-6, max(c(2 * ncol(X_proc), 2000))) so
## the upper end scales with the dimension on wide-X problems where
## the pre-1.7 default sigma = ncol(X) sat further up the range.
b_maxvarK <- function(X_proc, search_lower = 1e-6,
                      search_upper = NULL, tol = .Machine$double.eps^0.25) {
  if (!is.matrix(X_proc)) X_proc <- as.matrix(X_proc)
  if (is.null(search_upper)) {
    search_upper <- max(2 * ncol(X_proc), 2000)
  }
  if (!(search_lower > 0 && search_upper > search_lower))
    stop("search_lower must be positive and < search_upper")
  obj <- function(sigma) {
    .var_offdiag_K(gausskernel(X_proc, sigma = sigma))
  }
  res <- stats::optimize(obj, interval = c(search_lower, search_upper),
                         maximum = TRUE, tol = tol)
  list(sigma = res$maximum, var_K = res$objective,
       search_lower = search_lower, search_upper = search_upper)
}

## b_maxvarK_nystrom(): Nystrom-path bandwidth selector. Operates on
## the cross-kernel C = K(X, Z) between observations X and landmarks
## Z, both already preprocessed. Cheaper than the exact-path version
## (O(n m) per sigma vs O(n^2)) and produces a sigma that's
## well-calibrated to the same low-rank feature map the fit will use.
b_maxvarK_nystrom <- function(X_proc, Z_proc, search_lower = 1e-6,
                              search_upper = NULL,
                              tol = .Machine$double.eps^0.25) {
  if (!is.matrix(X_proc)) X_proc <- as.matrix(X_proc)
  if (!is.matrix(Z_proc)) Z_proc <- as.matrix(Z_proc)
  if (ncol(X_proc) != ncol(Z_proc))
    stop("X_proc and Z_proc must have the same number of columns")
  if (is.null(search_upper)) {
    search_upper <- max(2 * ncol(X_proc), 2000)
  }
  if (!(search_lower > 0 && search_upper > search_lower))
    stop("search_lower must be positive and < search_upper")
  obj <- function(sigma) {
    Cmat <- .nystrom_cross_kernel(X_proc, Z_proc, sigma)
    .var_C(Cmat)
  }
  res <- stats::optimize(obj, interval = c(search_lower, search_upper),
                         maximum = TRUE, tol = tol)
  list(sigma = res$maximum, var_C = res$objective,
       search_lower = search_lower, search_upper = search_upper)
}

## Plain n x m Gaussian cross-kernel exp(-||x_i - z_j||^2 / sigma).
## Used by the Nystrom bandwidth selector; the main Nystrom fit uses
## the equivalent helper in R/nystrom.R, but exposing this local
## version keeps b_maxvarK_nystrom self-contained for testing.
.nystrom_cross_kernel <- function(X, Z, sigma) {
  # ||x_i - z_j||^2 = ||x_i||^2 + ||z_j||^2 - 2 x_i' z_j
  Xrow2 <- rowSums(X * X)
  Zrow2 <- rowSums(Z * Z)
  D2 <- outer(Xrow2, Zrow2, "+") - 2 * tcrossprod(X, Z)
  D2[D2 < 0] <- 0  # numerical floor for ties
  exp(-D2 / sigma)
}
