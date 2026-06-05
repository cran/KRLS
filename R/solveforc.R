solveforc <- function(y = NULL, Eigenobject = NULL, lambda = NULL, eigtrunc = NULL) {

  nn <- nrow(y)

  if (is.null(eigtrunc)) {
    Ginv <- tcrossprod(multdiag(X = Eigenobject$vectors,
                                d = 1 / (Eigenobject$values + lambda)),
                       Eigenobject$vectors)
  } else {
    # eigentruncation: keep only eigenvectors at least 'eigtrunc' times as large as the largest.
    lastkeeper <- max(which(Eigenobject$values >= eigtrunc * Eigenobject$values[1]))
    # We need at least one.
    lastkeeper <- max(1, lastkeeper)

    Ginv <- tcrossprod(multdiag(X = Eigenobject$vectors[, 1:lastkeeper, drop = FALSE],
                                d = 1 / (Eigenobject$values[1:lastkeeper] + lambda)),
                       Eigenobject$vectors[, 1:lastkeeper, drop = FALSE])
  }

  coeffs <- tcrossprod(Ginv, t(y))
  Le     <- crossprod(coeffs / diag(Ginv))

  list(coeffs = coeffs, Le = Le)
}


