solveforc <-
function(y=NULL,Eigenobject=NULL,lambda=NULL){
   nn <- nrow(y)
   Ginv    <- tcrossprod(tcrossprod(Eigenobject$vectors,diag(c(1/(Eigenobject$values+lambda)),nn,nn)),Eigenobject$vectors)
   coeffs  <- tcrossprod(Ginv,t(y))
   Le      <- as.vector(crossprod(coeffs/diag(Ginv)))
   return(list(coeffs=coeffs,
               Le=Le)
               )
  }

