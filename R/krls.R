krls <-
function(     X=NULL,
              y=NULL,
              whichkernel="gaussian",
              lambda=NULL,
              sigma=NULL,
              derivative=TRUE,
              binary=TRUE,
              vcov=TRUE,
              print.level=1,
              L=NULL,
              U=NULL,
              tol=NULL){
              
      
      # checks
      y <- as.matrix(y)
      X <- as.matrix(X)
      
      if (is.numeric(X)==FALSE){
       stop("X must be numeric")
      }
      if (is.numeric(y)==FALSE){
       stop("y must be numeric")
      }
       if (sum(is.na(X))>0){
       stop("X contains missing data")
      }
      if (sum(is.na(y))>0){
       stop("y contains missing data")
      }
      if (var(y)==0){
       stop("y does not vary")
      }
      n <- nrow(X)
      d <- ncol(X)
      if (n!=nrow(y)){
       stop("nrow(X) not equal to number of elements in y")
      }
      
      stopifnot(
                is.logical(derivative),
                is.logical(vcov),
                is.logical(binary)
                )
      
      if(derivative==TRUE){
        if(vcov==FALSE){
        stop("derivative==TRUE requires vcov=TRUE")
        }
      }
      
      # default sigma to dim of X 
      if(is.null(sigma)) { sigma <- d
      } else {
        stopifnot(is.vector(sigma),
                  length(sigma)==1,
                  is.numeric(sigma),
                  sigma>0)        
      }
      
      # column names
      if(is.null(colnames(X))){
      colnames(X) <- paste("x",1:d,sep="")
      }
         
      # scale
      X.init <- X
      X.init.sd <- apply(X.init,2,sd)
      y.init <- y
      y.init.sd <- apply(y.init,2,sd)
      y.init.mean <- mean(y.init)
      X <- scale(X,center=TRUE,scale=X.init.sd)    
      y <- scale(y,center=y.init.mean,scale=y.init.sd)
     
      # kernel matrix
      K <- NULL
      if(whichkernel=="gaussian"){ K <- gausskernel(X,sigma=sigma)}
      if(whichkernel=="linear"){K <- tcrossprod(X)}
      if(whichkernel=="poly2"){K <- (tcrossprod(X)+1)^2}
      if(whichkernel=="poly3"){K <- (tcrossprod(X)+1)^3}
      if(whichkernel=="poly4"){K <- (tcrossprod(X)+1)^4}
      if(is.null(K)){stop("No valid Kernel specified")}
      # eigenvalue decomposition
      Eigenobject <- eigen(K,symmetric=TRUE)
      # default lamda is chosen by leave one out optimization 
       if(is.null(lambda)) {
 
 #  Old Way to get lambda       
 # first try with max eigenvalue (increase interval in case of corner solution at upper bound)
 #   lowerb  <-  .Machine$double.eps
 #   upperb <- max(Eigenobject$values)
 #     if(print.level>1) { cat("Using Leave one out validation to determine lamnda. Search Interval: 0 to",round(upperb,3), "\n")}  
 #     lambda <- optimize(looloss,interval=c(lowerb,upperb),y=y,Eigenobject=Eigenobject)$minimum
 #     if(lambda >= (upperb - .5)){
 #       if(print.level>1) { cat("Increasing search window for Lambda that minimizes Loo-Loss \n")}  
 #      lambda <- optimize(looloss,interval=c(upperb,2*upperb),y=y,Eigenobject=Eigenobject)$minimum
 #      }
      # run golden section search
       if(print.level>2){noisy <- TRUE
       } else {
         noisy <- FALSE
       }   
       lambda<- lambdasearch(L=L,U=U,y=y,Eigenobject=Eigenobject,noisy=noisy)
         
       if(print.level>1) { cat("Lambda that minimizes Loo-Loss is:",round(lambda,5),"\n")}    
       
       } else {  # check user specified lamnbda
         stopifnot(is.vector(lambda),
                   length(lambda)==1,
                   is.numeric(lambda),
                   lambda>0)  
      }
      # solve given LOO optimal or user specified lambda
      out <-  solveforc(y=y,Eigenobject=Eigenobject,lambda=lambda)
      # fitted values
      yfitted <- K%*%out$coeffs

      
      ## var-covar matrix for c
      if(vcov==TRUE){      
        # sigma squared
        sigmasq <- as.vector((1/n) * crossprod(y-yfitted))
        Gdiag   <-  Eigenobject$values+lambda
        Ginvsq  <-  tcrossprod(Eigenobject$vectors %*% diag(Gdiag^-2,n,n),Eigenobject$vectors)      
        vcovmatc <- sigmasq*Ginvsq
        # var-covar for y hats
        vcovmatyhat <- crossprod(K,vcovmatc%*%K) 

      } else { 
        vcov.c      <- NULL
        vcov.fitted <- NULL
      }
                
      # compute derivatives
      derivmat<-NULL
      avgderiv <- derivmat <- varavgderivmat <- NULL
      if(derivative==TRUE){
        if(whichkernel!="gaussian"){
          stop("derivatives are only available when whichkernel='gaussian' is specified")
        }
      rows <- cbind(rep(1:nrow(X), each = nrow(X)), 1:nrow(X))
      distances <- X[rows[,1],] - X[ rows[,2],]    # d by n*n matrix of pairwise distances  
      derivmat <- matrix(NA,n,d)
      varavgderivmat <- matrix(NA,1,d)
      colnames(derivmat)       <- colnames(X)
      colnames(varavgderivmat) <- colnames(X)
      for(k in 1:d){       
       if(d==1){
             distk <-  matrix(distances,n,n,byrow=TRUE)
         } else {
             distk <-  matrix(distances[,k],n,n,byrow=TRUE) 
        }
         L <-  distk*K
         # pointwise derivatives
         derivmat[,k] <- (-2/sigma)*L%*%out$coeff
         # variance for average derivative
         varavgderivmat[1,k] <- (1/n^2)*sum((-2/sigma)^2 * crossprod(L,vcovmatc%*%L))     
      }
       # avg derivatives
       avgderiv <- matrix(colMeans(derivmat),nrow=1)
       colnames(avgderiv) <- colnames(X)
      # get back to scale
      derivmat <- scale(y.init.sd*derivmat,center=FALSE,scale=X.init.sd)
      attr(derivmat,"scaled:scale")<- NULL
      avgderiv <- scale(as.matrix(y.init.sd*avgderiv),center=FALSE,scale=X.init.sd)
      attr(avgderiv,"scaled:scale")<- NULL
      varavgderivmat <- (y.init.sd/X.init.sd)^2*varavgderivmat
      attr(varavgderivmat,"scaled:scale")<- NULL
      }
         
      # get back to scale
      yfitted     <- yfitted*y.init.sd+y.init.mean
      vcov.c      <- (y.init.sd^2)*vcovmatc
      vcov.fitted <- (y.init.sd^2)*vcovmatyhat
      Looe        <- out$Le*y.init.sd
      # R square
      R2 <- 1-(var(y.init-yfitted)/(y.init.sd^2))
            
      # indicator for binary predictors
      binaryindicator=matrix(FALSE,1,d)
      colnames(binaryindicator) <- colnames(X)
      
      # return
   z <- list(K=K,
             coeffs=out$coeffs,
             Looe=Looe,
             fitted=yfitted,
             X=X.init,
             y=y.init,
             sigma=sigma,
             lambda=lambda,
             R2 = R2,
             derivatives=derivmat,
             avgderivatives=avgderiv, 
             var.avgderivatives=varavgderivmat,
             vcov.c=vcov.c,
             vcov.fitted=vcov.fitted,
             binaryindicator=binaryindicator
            )          
  class(z) <- "krls"
    
  # add first differences if requested
    if(derivative==TRUE && binary==TRUE){
      z <- fdskrls(z)
    }
  
  # printing
      if(print.level>0 && derivative==TRUE){
      output <- setNames(as.vector(z$avgderivatives), colnames(z$avgderivatives))
      cat("\n Average Marginal Effects:\n \n")
      print(output)
      
      cat("\n Quartiles of Marginal Effects:\n \n")
      print(apply(z$derivatives,2,quantile,probs=c(.25,.5,.75)))
      }
      
  return(z)
            
}

