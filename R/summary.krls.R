summary.krls <-
function(object,binary=TRUE,probs=seq(.25,.75,.25),...)
      {
            
        if( class(object)!= "krls" ){
        warning("Object not of class 'krls'")
        UseMethod("summary")
        return(invisible(NULL))
        }
        
        d <- ncol(object$X)
        n <- nrow(object$X)
        
        coefficients <- matrix(NA,d,0)
        if(is.null(colnames(object$X))){
           rownames(coefficients) <- paste("x",1:d,sep="") 
          } else {
           rownames(coefficients) <- colnames(object$X) 
          }      
       
        if(is.null(object$derivatives)){
          cat("* *********************** *\n")
          cat("Model Summary:\n\n")
          cat("R2:",object$R2,"\n\n")
          cat("\n")
          cat("recompute krls object with krls(...,derivative = TRUE) to get summary of marginal effects\n")
          qcoefficients <- 
          dydx <- NULL
        } else {
        # store marginal effects to be returned
        dydx <-   object$derivatives
          
        # average marginal effects  
        est     <- t(object$avgderivatives)
        se     <- sqrt(t(object$var.avgderivatives))
        tval   <- est/se
        avgcoefficients <- cbind(est, se, tval, 2 * pt(abs(tval),n-d, lower.tail = FALSE))
        colnames(avgcoefficients) <- c("Est", "Std. Error", "t value", "Pr(>|t|)")
        
        coefficients <- cbind(avgcoefficients,coefficients)
 
        cat("* *********************** *\n")
        cat("Model Summary:\n\n")
        cat("R2:",object$R2,"\n\n")
        
        
        # first differences for binary variables
        binaryindicator <- NULL
        if(binary){
        lengthunique    <- function(x){length(unique(x))}
        # vector with positions of binary variables
        binaryindicator <-which(apply(object$X,2,lengthunique)==2)
        
        if(length(binaryindicator)==0){
          # no binary vars in X
         } else {
          # compute marginal differences from min to max 
          est <- se <- matrix(NA,nrow=length(binaryindicator),ncol=1)
          qest <- matrix(NA,nrow=length(binaryindicator),ncol=length(probs))
          diffsstore <- matrix(NA,nrow=n,ncol=length(binaryindicator))
          for(i in 1:length(binaryindicator)){
            X1 <- X0 <- object$X
            # test data with D=Max
            X1[,binaryindicator[i]] <- max(X1[,binaryindicator[i]])
            # test data with D=Min
            X0[,binaryindicator[i]] <- min(X0[,binaryindicator[i]])
            Xall      <- rbind(X1,X0)
            # contrast vector
            h         <- matrix(rep(c(1/n,-(1/n)),each=n),ncol=1)
            # fitted values
            pout      <- predict(object,newdata=Xall,se=TRUE)
            # diffs in means
            est[i,1] <- t(h)%*%pout$fit        
            # SE (multiply by sqrt2 to correct for using data twice )
            se[i,1] <- as.vector(sqrt(t(h)%*%pout$vcov.fit%*%h))*sqrt(2)
            # quantiles
            diffs <- pout$fit[1:n]-pout$fit[(n+1):(2*n)]
            qest[i,] <-  quantile(diffs,probs=probs)
            # store FD estimates for post hoc regressions
            diffsstore[,i] <- diffs 
          }
        tval   <- est/se
        fds    <- cbind(est, se, tval, 2 * pt(abs(tval),n-d, lower.tail = FALSE))
        coefficients[binaryindicator,1:4] <- fds
        rownames(coefficients)[binaryindicator] <- paste(rownames(coefficients)[binaryindicator],"*",sep="")
        dydx[,binaryindicator] <- diffsstore   
        }  
        } # close binary = TRUE
        
        cat("Average Marginal Effects:\n")
        print(coefficients)
        if(binary && is.null(binaryindicator)==FALSE){
        cat("\n(*) average dy/dx is for discrete change of dummy variable from min to max (i.e. usually 0 to 1))\n\n")
        }
        
        # quantiles derivatives
          qcoefficients <-  t(apply(object$derivatives,2,quantile,probs=probs))
           if(is.null(colnames(object$X))){
              rownames(qcoefficients)<- paste("x",1:d,sep="") 
             } else {
              rownames(qcoefficients) <- colnames(object$X)
             }
        # quantile of first differences if applicable    
        if(binary && is.null(binaryindicator)==FALSE){ 
         qcoefficients[binaryindicator,] <- qest
         rownames(qcoefficients)[binaryindicator] <- paste(rownames(qcoefficients)[binaryindicator],"*",sep="")
        }
          #
        cat("\n")
        cat("Quartiles of Marginal Effects:\n")
        print(qcoefficients)
        
        if(binary && is.null(binaryindicator)==FALSE){ 
         cat("\n(*) quantiles of dy/dx is for discrete change of dummy variable from min to max (i.e. usually 0 to 1))\n\n")
        }
        
        }
        
             
      ans <- list(
                coefficients=coefficients,
                qcoefficients=qcoefficients,
                dydx=dydx)
      class(ans) <- "summary.krls"  
      return(invisible(ans))
}

