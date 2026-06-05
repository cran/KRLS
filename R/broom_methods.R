# tidy() / glance() / augment() methods for krls objects.
# Registered against the generics in the `generics` package
# (which is what `broom` re-exports). Discoverable via library(broom).

tidy.krls <-
function(x, conf.int = TRUE, conf.level = 0.95, ...)
  {
    if (!inherits(x, "krls"))
      stop("x is not of class 'krls'")
    if (is.null(x$derivatives))
      stop("\n tidy() requires marginal effects; refit with krls(..., derivative = TRUE)\n")

    est <- as.numeric(x$avgderivatives)
    n   <- nrow(x$X); d <- ncol(x$X)
    df  <- max(n - d, 1)

    # Standard errors / inference are only available when derivative
    # variances were computed. If they are absent (e.g. vcov=FALSE),
    # keep stable inference columns filled with NA.
    has_se <- !is.null(x$var.avgderivatives)
    if (has_se) {
      se   <- sqrt(as.numeric(x$var.avgderivatives))
      tval <- est / se
      pval <- 2 * stats::pt(abs(tval), df, lower.tail = FALSE)
    } else {
      se   <- rep(NA_real_, length(est))
      tval <- rep(NA_real_, length(est))
      pval <- rep(NA_real_, length(est))
    }

    out <- data.frame(
      term      = colnames(x$X),
      estimate  = est,
      std.error = se,
      statistic = tval,
      p.value   = pval,
      stringsAsFactors = FALSE
    )

    if (isTRUE(conf.int)) {
      if (has_se) {
        crit <- stats::qt(1 - (1 - conf.level) / 2, df)
        out$conf.low  <- est - crit * se
        out$conf.high <- est + crit * se
      } else {
        out$conf.low  <- NA_real_
        out$conf.high <- NA_real_
      }
    }

    # Pointwise distribution summaries (the unique KRLS contribution)
    qd <- apply(x$derivatives, 2, stats::quantile,
                probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
    out$q25    <- qd[1, ]
    out$median <- qd[2, ]
    out$q75    <- qd[3, ]

    if (!is.null(x$binaryindicator)) {
      out$binary <- as.logical(x$binaryindicator)
    }

    rownames(out) <- NULL
    out
  }

glance.krls <-
function(x, ...)
  {
    if (!inherits(x, "krls"))
      stop("x is not of class 'krls'")
    n <- nrow(x$X); d <- ncol(x$X)
    # Effective degrees of freedom = tr(S(lambda)).
    # Under approx = "nystrom" S(lambda) = U_p diag(Sigma2 / (Sigma2 + lambda))
    # U_p', and the cached Sigma2 from the SVD of Phi gives the trace
    # directly (the m-dimensional Nystrom spectrum). Under the exact path
    # we use the eigendecomposition of K stored on the fit object.
    eff_df <- if (!is.null(x$Sigma2) && !is.null(x$lambda)) {
      sum(x$Sigma2 / (x$Sigma2 + x$lambda))
    } else if (!is.null(x$K) && !is.null(x$lambda)) {
      ev <- tryCatch(eigen(x$K, symmetric = TRUE, only.values = TRUE)$values,
                     error = function(e) NA_real_)
      if (length(ev) == 1 && is.na(ev)) NA_real_
      else sum(ev / (ev + x$lambda))
    } else {
      NA_real_
    }
    looe_mean <- if (!is.null(x$Looe)) mean(x$Looe^2) else NA_real_

    is_nys     <- !is.null(x$approx) && x$approx == "nystrom"
    approx_lbl <- if (is.null(x$approx)) "none" else x$approx
    inference_lbl <- if (is.null(x$inference)) {
      if (!is.null(x$vcov.c)) "exact" else "none"
    } else x$inference

    data.frame(
      nobs         = n,
      n_predictors = d,
      r.squared    = x$R2,
      loo_mse      = looe_mean,
      lambda       = x$lambda,
      sigma        = x$sigma,
      eff_df       = eff_df,
      approx       = approx_lbl,
      nystrom_m    = if (is_nys) as.integer(x$nystrom_m) else NA_integer_,
      inference    = inference_lbl,
      stringsAsFactors = FALSE
    )
  }

augment.krls <-
function(x, data = NULL, ...)
  {
    if (!inherits(x, "krls"))
      stop("x is not of class 'krls'")
    n <- nrow(x$X)

    # In-sample fitted values (from the krls object); column-prefixed so
    # they don't collide with user data columns.
    .fitted   <- as.numeric(x$fitted)
    .resid    <- as.numeric(x$y) - .fitted
    derivs    <- as.matrix(x$derivatives)
    colnames(derivs) <- paste0(".dy_d_", colnames(x$X))

    if (is.null(data)) {
      out <- data.frame(.fitted = .fitted, .resid = .resid,
                        as.data.frame(x$X),
                        stringsAsFactors = FALSE)
    } else {
      if (NROW(data) != n)
        stop(sprintf(
          "\n augment(): supplied data has %d rows but the fit was built on %d \n",
          NROW(data), n))
      out <- data.frame(data, stringsAsFactors = FALSE)
      out$.fitted <- .fitted
      out$.resid  <- .resid
    }
    out <- cbind(out, derivs)
    rownames(out) <- NULL
    out
  }
