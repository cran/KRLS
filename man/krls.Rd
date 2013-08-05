\name{krls}
\alias{krls}
\title{
Kernel-based Regularized Least Squares (KRLS)
}
\description{
Function that implements Kernel-based Regularized Least Squares (KRLS), a machine learning method that fits 
flexible solution surfaces of the form y=f(X) that arise in regression or classification problems 
without relying on linearity or other assumptions that use the columns of the predictor matrix X directly as basis functions (such as additivity). KRLS finds the best fitting function by minimizing a Tikhonov regularization problem 
with a square loss using Gaussian Kernels as radial basis functions (see details). 

}
\usage{
krls(X = NULL, y = NULL, whichkernel = "gauss", lambda = NULL,
sigma = NULL, derivative = TRUE, print.level = 0)
}

\arguments{
  \item{X}{
A n by k data numeric matrix that contains the values of k predictor variables for i=1,...,n observations. The matrix may not contain missing values. 
}
  \item{y}{
A n by 1 data numeric matrix or vector that contains the values of response variable for all observations. This vector may not contain missing values. 
}
  \item{whichkernel}{
Which kernel should be used. Must be one of: "gauss","linear","poly1","poly2","poly3", or "poly4" (see details). Default is Gaussian Kernel.
}
  \item{lambda}{
A positive scalar that specifies the parameter for the regularizer. It governs the tradeoff between 
model fit and complexity. By default lambda is chosen by minimizing the sum of the squared leave-one-out errors.
}
  \item{sigma}{
Optional parameter to specify the width of the Gaussian kernel (see \code{\link{gausskernel}} for details). 
By default sigma is set equal to k (the number of dimensions) which ensures a reasonable scaling of the distances on average. Future releases 
may include an option for choosing sigma based on cross-validation.
}
  \item{derivative}{
Logical that specifies whether derivatives should be computed. Currently, derivatives are only implemented for the Gaussian Kernel.
}
  \item{print.level}{
Positive integer that determines the level of printing.
}
}
\details{
Kernel-based Regularized Least Squares arises
as a Tikhonov minimization problem with a square loss, ie. we search over a 
space of possible functions H and choose the best function f according to the rule: 

argmin_{f in H} sum_i^n (y_i - f(x_i))^2 + lambda ||f||_H^2

where (y_i - f(x_i))^2 a loss function that computes how `wrong' the function 
is at each observation i and lambda ||f||_H^2 is the regularizer that measures the 
complexity of the function according to the L_2 norm |f||^2 = integral f(x)^2 dx. Note that
lambda is the parameter that governs the complexity of the tradeoff between model fit and complexity. 

The representer theorem states that under fairly general conditions, the function that minimizes the regularized loss 
within the hypothesis space established by the choice of a (positive semidefinite) kernel function k is of the form:
f(x_j)= sum_i^n c_i k(x_i,x_j) where the kernel function k(x_i,x_j) measures the distance 
between two observations x_i and x_j, K is the kernel matrix with all 
pairwise distances K_ij=k(x_i,x_j), and c is the vector of choice coefficients for each observation i such that y=Kc. Accordingly, the krls function solves
the following minimization problem:

argmin_{f in H} sum_i^n (y - Kc)'(y-Kc)+ lambda c'Kc

which is convex in c and solved by c=(K +lambda I)^-1 y, a linear solution that 
provides a fit that will be potentially highly non-linear in terms of the predictors.

The function currently implements KRLS for the following Kernels:

"gaussian": k(x_i,x_j)=exp(-|| x_i - x_j ||^2 / sigma^2) where ||x_i - x_j|| is the Euclidean distance.

"linear": k(x_i,x_j)=x_i'x_j

"poly1", "poly2", "poly3", "poly4" which are polynomials based on  k(x_i,x_j)=(x_i'x_j +1)^p where p is the order.


For the Gaussian kernel the sigma width parameter is set to the number of dimensions by default. Unless otherwise specified by the user, the 
lambda parameter is chosen by minimization of the leave-one-out error. 

The function also computes the variance-covariance matrix for the choice coefficients c and fitted values y=Kc based on 
a variance estimator developed in Hainmueller and Hazlett (2011).

For the Gaussian kernel, the function also computes the pointwise partial derivatives of the fitted function wrt to each predictor in X if \code{derivative=TRUE} is specified. 
Average derivatives are also computed with variances.

}
\value{
A list object of class "krls" with the following elements:
  \item{K}{n by n matrix of pairwise kernel distances between observations.}
  \item{coeffs }{n by 1 vector of choice coefficients c.}
  \item{Le}{scalar with sum of squared leave-one-out error.}
  \item{fitted}{n by 1 vector of fitted values.}
  \item{X}{original n by k data matrix.}
  \item{y}{original n by 1 outcome matrix.} 
  \item{sigma}{sigma value used for Gaussian Kernel.}   
  \item{lambda}{lambda value used (user specified or based on LOO minimization).}  
  \item{R2}{scalar with R squared value}
  \item{vcov.c}{n by n variance covariance matrix for choice coefficients.}
  \item{vcov.fitted}{n by n variance covariance matrix for fitted values.}  
  \item{derivatives}{n by k matrix of pointwise derivatives based on gaussian kernel (NULL unless \code{derivative=TRUE} is specified.}
  \item{avgderivatives}{1 by k matrix of average derivative based on gaussian kernel (NULL unless \code{derivative=TRUE} is specified.}
   \item{var.avgderivatives}{1 by k matrix of variances for average derivative based on gaussian kernel (NULL unless \code{derivative=TRUE} is specified.}  

  }
\references{

Hainmueller, J. and Hazlett, C. 2011. Kernel Regularized Least Squares: Moving Beyond Linearity and Additivity Without Sacrifcing Interpretability. Working Paper MIT.

Rifkin, R. 2002. Everything Old is New Again: A fresh look at historical approaches in machine learning. Thesis, MIT. September, 2002.

Evgeniou, T., Pontil, M., and Poggio, T. (2000). Regularization networks and support vector
machines. Advanced In Computational Mathematics, 13(1):1-50.

Schoelkopf, B., Herbrich, R. and Smola, A.J. (2001) A generalized representer theorem. In 14th Annual
Conference on Computational Learning Theory, pages 416-426.

Kimeldorf, G.S. Wahba, G. 1971. Some results on Tchebycheffian spline functions. Journal of Mathematical Analysis and Applications, 33:82-95.

}
\author{
Jens Hainmueller (MIT) and Chad Hazlett (MIT)
}
\note{
The function requires the storage of a n by n kernel matrix and therefore exceed the memory limits for very large datasets (eg. more than 10,000 datapoints).
}


\seealso{
\code{\link{predict.krls}} for fitted values and predictions.
}
\examples{

# Linear example
# set up data
N <- 200
x1 <- rnorm(N)
x2 <- rbinom(N,size=1,prob=.2)
y <- x1 + .5*x2 + rnorm(N,0,.15)
X <- cbind(x1,x2)
# fit model
krlsout <- krls(X=X,y=y)
# summarize marginal effects and contribution of each variable
summary(krlsout)
# plot marginal effects and conditional expectation plots
plot(krlsout)


# non-linear example
# set up data
N <- 200
x1 <- rnorm(N)
x2 <- rbinom(N,size=1,prob=.2)
y <- x1^3 + .5*x2 + rnorm(N,0,.15)
X <- cbind(x1,x2)

# fit model
krlsout <- krls(X=X,y=y)
# summarize marginal effects and contribution of each variable
summary(krlsout)
# plot marginal effects and conditional expectation plots
plot(krlsout)

## 2D example:
# predictor data
X <- matrix(seq(-3,3,.1))
# true function
Ytrue <- sin(X)
# add noise 
Y     <- sin(X) + rnorm(length(X),sd=.3)
# approximate function using KRLS
out <- krls(y=Y,X=X)
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
