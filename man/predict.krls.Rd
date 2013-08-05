\name{predict.krls}
\alias{predict.krls}
\title{
Predict method for Kernel-based Regularized Least Squares (KRLS) Model Fits
}
\description{
Predicted values and standard errors based on krls model object.
}
\usage{

 \method{predict}{krls}(object, newdata, se.fit = FALSE , \dots)
}
\arguments{
  \item{object}{Fitted krls model, i.e. an object of class krls}
  \item{newdata}{An n by k data frame or matrix with variables values at which to predict the outcome. Variable columns in newdata have to match the corresponding covariates 
used in the krls model object.
}

\item{se.fit}{logical flag if standard errors should be computed for pointwise predictions.}

\item{\dots}{additional arguments affecting the predictions produced.}

}
\details{

Function produces predicted values, obtained by evaluating the fitted krls function with 
the newdata (ie. the test points). The prediction at a new test point x_i is based on f(x_i)= sum_j=1^n c_j K_{x_j}(x_i) where K is the kernel matrix and thus K_{x_j}
is a vector whose j-th entry is K(x_j,x_i) (e.g. the distance between the test point x_i and the training point x_j). The training points are passed to the function
with the krls object.



}
\value{

  \item{fit}{n by 1 vector of fitted values for n test points.}
  \item{se.fit}{n by 1 vector of standard errors for the fitted values for n test points (NULL unless se.fit=TRUE is specified).}
  \item{vcov.fit}{n by n variance-covariance matrix for the fitted values for n test points (NULL unless se.fit=TRUE is specified).}
  \item{newdata}{n by k data matrix for test points}
  \item{newdataK}{n by m data matrix for pairwise Gauss Kernel distances between test points and m training points from krls model fit.}

}

\author{
Jens Hainmueller (MIT) and Chad Hazlett (MIT)
}

\seealso{
\code{\link{krls}}
}
\examples{


## 2D example:
# predictor data
X <- matrix(seq(-3,3,.1))
# true function
Ytrue <- sin(X)
# add noise 
Y     <- sin(X) + rnorm(length(X),sd=.3)
# approximate function using KRLS
out <- krls(y=Y,X=X,derivative=TRUE)
# get fitted values and ses
fit <- predict.krls(out,newdata=X,se.fit=TRUE)
# results
par(mfrow=c(2,1))
plot(y=Ytrue,x=X,type="l",col="red",ylim=c(-1.2,1.2),lwd=2,main="f(x)")
points(y=fit$fit,X,col="blue",pch=19)
arrows(y1=fit$fit+1.96*fit$se.fit,
       y0=fit$fit-1.96*fit$se.fit,
       x1=X,x0=X,col="blue",length=0)
legend("bottomright",legend=c("true f(x)=sin(x)","KRLS fitted f(x)"),
       lty=c(1,NA),pch=c(NA,19),lwd=c(2,NA),col=c("red","blue"),cex=.8)

plot(y=cos(X),x=X,type="l",col="red",ylim=c(-1.2,1.2),lwd=2,main="df(x)/dx")
points(y=out$derivatives,X,col="blue",pch=19)

legend("bottomright",legend=c("true df(x)/dx=cos(x)","KRLS fitted df(x)/dx"),
       lty=c(1,NA),pch=c(NA,19),lwd=c(2,NA),col=c("red","blue"),,cex=.8)

## 3D example
# plot true function
par(mfrow=c(1,2))
f<-function(x1,x2){ sin(x1)*cos(x2)}
x1 <- x2 <-seq(0,2*pi,.2)
z   <-outer(x1,x2,f)
persp(x1, x2, z,theta=30,main="true f(x1,x2)=sin(x1)cos(x2)")
# approximate function with KRLS
# data and outcomes
X <- cbind(sample(x1,200,replace=TRUE),sample(x2,200,replace=TRUE))
y   <- f(X[,1],X[,2])+ runif(nrow(X))
# fit surface
krlsout <- krls(X=X,y=y)
# plot fitted surface
ff  <- function(x1i,x2i,krlsout){predict(object=krlsout,newdata=cbind(x1i,x2i))$fit}
z   <- outer(x1,x2,ff,krlsout=krlsout)
persp(x1, x2, z,theta=30,main="KRLS fitted f(x1,x2)")


}
