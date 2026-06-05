lambdasearch <-
  function(L=NULL,
           U=NULL,
           y=NULL,
           Eigenobject=NULL,
           tol=NULL,
           noisy=FALSE,
           eigtrunc=NULL,
           lambda_method=c("loo","gcv")){

    lambda_method <- match.arg(lambda_method)

    n <- nrow(y)
    if(is.null(tol)){
      tol <- 10^-3 * n
    } else {
      stopifnot(is.vector(tol),
                length(tol)==1,
                is.numeric(tol),
                tol>0)
    }

  # get upper bound starting value
    if(is.null(U)){
    U <- n
    while(sum(Eigenobject$values / (Eigenobject$values + U)) < 1){
      U <- U-1
     }
    } else {
      stopifnot(is.vector(U),
                length(U)==1,
                is.numeric(U),
                U>0)
    }

  # get lower bound starting value
    if(is.null(L)){
      q <- which.min(abs(Eigenobject$values - (max(Eigenobject$values)/1000)))

      #L <- 0
      L = .Machine$double.eps  #CJH: to avoid Inf in next statement

      while(sum(Eigenobject$values / (Eigenobject$values + L)) > q){
        L <- L+.05
      }
    }  else {
      stopifnot(is.vector(L),
                length(L)==1,
                is.numeric(L),
                L>=0)
    }

    # Loss-function dispatch. LOO objective is the historical default; GCV
    # is the closed-form generalized-CV alternative -- comparable in cost
    # since both reuse the eigendecomposition.
    loss_fn <- if (lambda_method == "loo") {
      function(lam) as.numeric(looloss(lambda = lam, y = y,
                                       Eigenobject = Eigenobject,
                                       eigtrunc = eigtrunc))
    } else {
      function(lam) .gcv_loss_exact(lam, Eigenobject, y, eigtrunc)
    }

    # create new search values
    X1 <- L + (.381966)*(U-L)
    X2 <- U - (.381966)*(U-L)

    S1 <- loss_fn(X1)
    S2 <- loss_fn(X2)

    if(noisy){cat("L:",L,"X1:",X1,"X2:",X2,"U:",U,"S1:",S1,"S2:",S2,"\n") }

    while(abs(S1-S2)>tol){ # terminate if difference between S1 and S2 less than tolerance

     # update steps and use caching
      if(S1 < S2){
        U  <- X2
        X2 <- X1
        X1 <- L + (.381966)*(U-L)
        S2 <- S1
        S1 <- loss_fn(X1)

       } else { #S2 < S1
        L  <- X1
        X1 <- X2
        X2 <- U - (.381966)*(U-L)
        S1 <- S2
        S2 <- loss_fn(X2)
       }

      if(noisy){cat("L:",L,"X1:",X1,"X2:",X2,"U:",U,"S1:",S1,"S2:",S2,"\n") }
    }
    out <- if (S1 < S2) X1 else X2
    if (noisy) cat("Lambda:", out, "\n")
    return(invisible(out))
}

## GCV objective for the exact path. GCV(lambda) = RSS(lambda) /
## (1 - tr(S(lambda))/n)^2 with S(lambda) the kernel-ridge hat matrix.
## Computed directly from the eigendecomposition of K without
## materializing Ginv or S explicitly.
.gcv_loss_exact <- function(lambda, Eigenobject, y, eigtrunc = NULL) {
  vals <- Eigenobject$values
  vecs <- Eigenobject$vectors
  if (!is.null(eigtrunc)) {
    last <- max(which(vals >= eigtrunc * vals[1]))
    last <- max(1L, last)
    vals <- vals[1:last]
    vecs <- vecs[, 1:last, drop = FALSE]
  }
  yv     <- as.numeric(y)
  n      <- length(yv)
  weights <- vals / (vals + lambda)
  yhat   <- as.numeric(vecs %*% (weights * crossprod(vecs, yv)))
  RSS    <- sum((yv - yhat)^2)
  tr_S   <- sum(weights)
  denom  <- 1 - tr_S / n
  if (denom <= .Machine$double.eps) return(.Machine$double.xmax)
  RSS / (denom^2)
}
