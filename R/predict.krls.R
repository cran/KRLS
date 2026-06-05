predict.krls <-
function(object, newdata, se.fit = FALSE, ...) {

  if (!inherits(object, "krls"))
    stop("object is not of class 'krls'")

  is_nystrom <- !is.null(object$landmarks)

  if (isTRUE(se.fit) && is.null(object$vcov.c))
    stop("refit with krls(..., vcov = TRUE) to compute standard errors")

  has_prep <- !is.null(object$prep)

  if (has_prep) {
    # Modern path: object$X_proc is already in kernel space (continuous
    # standardized, categorical sqrt(0.5)-one-hot), and object$prep
    # holds the transformation to apply to newdata.
    X_anchor <- object$X_proc
    prep <- object$prep

    # Dual-dispatch on newdata. fdskrls() builds its X1/X0 by
    # modifying object$X (the raw input form, which may be a
    # data.frame) and passes that to predict() -- those cases
    # naturally route through .apply_preprocess_X(). The "already
    # processed" branch supports callers that have done the encoding
    # themselves (and matches the legacy contract where newdata could
    # be a fully numeric matrix in kernel-input units).
    nd_cols <- if (is.data.frame(newdata)) ncol(newdata) else
      ncol(as.matrix(newdata))
    if (nd_cols == prep$d_orig) {
      newdata_proc <- .apply_preprocess_X(newdata, prep)
    } else if (nd_cols == ncol(X_anchor)) {
      newdata_proc <- as.matrix(newdata)
      storage.mode(newdata_proc) <- "double"
    } else {
      stop("ncol(newdata) = ", nd_cols,
           " matches neither the raw input (", prep$d_orig,
           ") nor the processed (", ncol(X_anchor),
           ") column count of the training X")
    }
  } else {
    # Legacy path: object$X is the raw input matrix, no preprocessing
    # metadata stored. Reproduce the original predict.krls() logic by
    # re-standardizing object$X and newdata together.
    newdata <- as.matrix(newdata)
    if (ncol(object$X) != ncol(newdata))
      stop("ncol(newdata) differs from ncol(X) from fitted krls object")
    Xmeans <- colMeans(object$X)
    Xsd    <- apply(object$X, 2, sd)
    X_anchor <- scale(object$X, center = Xmeans, scale = Xsd)
    attr(X_anchor, "scaled:center") <- NULL
    attr(X_anchor, "scaled:scale")  <- NULL
    newdata_proc <- scale(newdata, center = Xmeans, scale = Xsd)
    attr(newdata_proc, "scaled:center") <- NULL
    attr(newdata_proc, "scaled:scale")  <- NULL
  }

  nn <- nrow(newdata_proc)
  if (is_nystrom) {
    anchors  <- object$landmarks       # already in standardized/processed space
    newdataK <- .gauss_cross_kernel(newdata_proc, anchors, object$sigma)
  } else {
    anchors  <- X_anchor
    newdataK <- matrix(
      gausskernel(rbind(newdata_proc, anchors),
                  sigma = object$sigma)[1:nn, (nn + 1):(nn + nrow(anchors))],
      nrow = nn, byrow = FALSE)
  }

  yfitted <- newdataK %*% object$coeffs

  if (se.fit) {
    vcov.c.raw  <- object$vcov.c * as.vector(1 / var(object$y))
    vcov.fitted <- tcrossprod(newdataK %*% vcov.c.raw, newdataK)
    vcov.fit    <- (apply(object$y, 2, sd)^2) * vcov.fitted
    se.fit      <- matrix(sqrt(diag(vcov.fit)), ncol = 1)
  } else {
    vcov.fit <- se.fit <- NULL
  }

  yfitted <- (yfitted * apply(object$y, 2, sd)) + mean(object$y)

  list(fit      = yfitted,
       se.fit   = se.fit,
       vcov.fit = vcov.fit,
       newdata  = newdata_proc,
       newdataK = newdataK)
}
