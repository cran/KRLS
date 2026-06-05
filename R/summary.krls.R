summary.krls <-
function(object, probs=c(.25,.5,.75),...)
      {

        if (!inherits(object, "krls")) {
          stop("object is not of class 'krls'")
        }

        cat("* *********************** *\n")
        cat("Model Summary:\n\n")
        cat("R2:",object$R2,"\n\n")

        d <- ncol(object$X)
        n <- nrow(object$X)

        # Nystrom-specific block: surface the approximation choices and
        # any diagnostic concerns (e.g. many landmark eigenvalues at the
        # relative-ridge floor, which signals near-collinear landmarks).
        if (!is.null(object$approx) && object$approx == "nystrom") {
          method_lbl <- switch(if (is.null(object$landmark_method)) "unknown" else object$landmark_method,
            random       = "random",
            kmeans       = "k-means",
            user_indices = "user-supplied indices",
            user_matrix  = "user-supplied matrix",
            object$landmark_method)
          cat("Approximation: Nystrom with m =", object$nystrom_m,
              "landmarks (", method_lbl, ")\n")
          cat("Inference:", switch(if (is.null(object$inference)) "none" else object$inference,
            conditional_nystrom = "conditional approximate",
            none                = "point estimates only (vcov = FALSE)",
            object$inference), "\n")
          if (!is.null(object$D_max_raw) && !is.null(object$D_min_raw) &&
              is.finite(object$D_max_raw) && object$D_max_raw > 0) {
            cond_str <- if (object$D_min_raw <= 0) {
              "rank-deficient"
            } else {
              sprintf("%.2g", object$D_max_raw / object$D_min_raw)
            }
            cat("Landmark kernel: condition =", cond_str, ",",
                object$floored_count, "of", object$nystrom_m,
                "eigenvalues at relative-ridge floor\n")
            # Loud-but-not-fatal warning if a large fraction were floored.
            if (!is.null(object$floored_count) &&
                object$floored_count > object$nystrom_m / 2) {
              cat("  Note: more than half of the landmark-kernel ",
                  "eigenvalues are at the relative-ridge floor; ",
                  "consider larger sigma, fewer landmarks, or more ",
                  "distinct landmark coordinates.\n", sep = "")
            }
          }
          cat("\n")
        }

        coefficients <- matrix(NA,d,0)
        rownames(coefficients) <- colnames(object$X)


        if(is.null(object$derivatives)){
          cat("\n")
          cat("recompute krls object with krls(...,derivative = TRUE) to get summary of marginal effects\n")
          return(invisible(NULL))
        }

        # average marginal effects
        est     <- t(object$avgderivatives)
        has_se  <- !is.null(object$var.avgderivatives)
        if (has_se) {
          se   <- sqrt(t(object$var.avgderivatives))
          tval <- est/se
          avgcoefficients <- cbind(est, se, tval,
                                   2 * pt(abs(tval), n - d, lower.tail = FALSE))
          colnames(avgcoefficients) <- c("Est", "Std. Error", "t value", "Pr(>|t|)")
        } else {
          # No variances available (e.g. approx = "nystrom"); print point
          # estimates only.
          avgcoefficients <- est
          colnames(avgcoefficients) <- "Est"
        }

       # add stars for binary
        if(sum(object$binaryindicator)>0){
          rownames(avgcoefficients)[object$binaryindicator] <- paste(rownames(avgcoefficients)[object$binaryindicator],"*",sep="")
        }

        cat("Average Marginal Effects:\n")
        print(avgcoefficients,...)
        if(sum(object$binaryindicator)>0){
        cat("\n(*) average dy/dx is for discrete change of dummy variable from min to max (i.e. usually 0 to 1))\n\n")
        }

        # quantiles of derivatives
        qderiv <- apply(object$derivatives,2,quantile,probs=probs)
        if(sum(object$binaryindicator)>0){
          colnames(qderiv)[object$binaryindicator] <- paste(colnames(qderiv)[object$binaryindicator],"*",sep="")
        }
        qderiv <- t(qderiv)

        cat("\n")
        cat("Quartiles of Marginal Effects:\n")
        print(qderiv,...)

        if(sum(object$binaryindicator)>0){
         cat("\n(*) quantiles of dy/dx is for discrete change of dummy variable from min to max (i.e. usually 0 to 1))\n\n")
        }

      ans <- list(
                 coefficients=avgcoefficients,
                 qcoefficients=qderiv)
      class(ans) <- "summary.krls"
      return(invisible(ans))
}

