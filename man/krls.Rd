\name{krls}
\alias{krls}
\title{
Kernel-based Regularized Least Squares (KRLS)
}
\description{
Function that implements Kernel-based Regularized Least Squares (KRLS), a machine learning method described in Hainmueller and Hazlett (2013) that fits flexible solution surfaces of the form y=f(X) that arise in regression or classification problems without relying on linearity, additivity, or other assumptions that use the columns of the predictor matrix X directly as basis functions. KRLS finds the best fitting function by minimizing a Tikhonov regularization problem with a squared loss using Gaussian Kernels as radial basis functions (see details). 

}
\usage{
krls(X = NULL, y = NULL, whichkernel = "gaussian", lambda = NULL,
sigma = NULL, derivative = TRUE, binary= TRUE, vcov=TRUE, 
print.level = 1)
}

\arguments{
  \item{X}{
\var{N} by \var{D} data numeric matrix that contains the values of \var{D} predictor variables for \eqn{i=1,\ldots,N} observations. The matrix may not contain missing values. 
}
  \item{y}{
\var{N} by \var{1} data numeric matrix or vector that contains the values of the response variable for all observations. This vector may not contain missing values. 
}
  \item{whichkernel}{
String vector that specifies which kernel should be used. Must be one of \code{gaussian}, \code{linear}, \code{poly1}, \code{poly2}, \code{poly3}, or \code{poly4} (see details). Default is \code{gaussian}.
}
  \item{lambda}{
A positive scalar that specifies the \eqn{\lambda}{lambda} parameter for the regularizer (see details). It governs the tradeoff between model fit and complexity. By default, this parameter is chosen by minimizing the sum of the squared leave-one-out errors.
}
  \item{sigma}{
A positive scalar that specifies the bandwidth of the Gaussian kernel (see \code{\link{gausskernel}} for details). 
By default, the bandwidth is set equal to \var{D} (the number of dimensions) which ensures a reasonable scaling of average distance between observations in the standardized data that is used for the fitting. 
}
  \item{derivative}{
Logical that specifies whether pointwise partial derivatives should be computed. Currently, derivatives are only implemented for the Gaussian Kernel.
}
  \item{binary}{
Logical that specifies whether first-differences instead of pointwise partial derivatives should be computed for binary predictors. Ignored unless \code{derivative=TRUE}. 
}
  \item{vcov}{
Logical that specifies whether variance-covariance matrix for the choice coefficients \var{c} and fitted values should be computed. Note that \code{derivative=TRUE} requires that \code{vcov=TRUE}. 
}

\item{print.level}{
Positive integer that determines the level of printing. Set to 0 for no printing and 2 for more printing.
}
}
\details{

\code{krls} implements the Kernel-based Regularized Least Squares (KRLS) estimator as described in Hainmueller and Hazlett (2013). Please consult this reference for any details.

Kernel-based Regularized Least Squares (KRLS) arises as a Tikhonov minimization problem with a squared loss. Assume we have data of the from \eqn{y_i,\,x_i}{y_i, x_i} where \var{i} indexes observations, \eqn{y_i \in R}{y_i in R} is the outcome and \eqn{x_i \in R^D}{x_i in R^D} is a D-dimensional vector of predictor values. Then KRLS searches over a space of functions \eqn{H} and chooses the best fitting function \eqn{f} according to the rule: 

\deqn{argmin_{f \in H} \sum_i^N (y_i - f(x_i))^2 + \lambda ||f||_{H^2}}{%
      argmin_{f in H} sum_i^N (y_i - f(x_i))^2 + lambda || f ||_H^2} 

where \eqn{(y_i - f(x_i))^2} is a loss function that computes how `wrong' the function 
is at each observation \var{i} and \eqn{|| f ||_{H^2}}{|| f ||_H^2} is the regularizer that measures the complexity of the function according to the \eqn{L_2} norm \eqn{||f||^2 = \int f(x)^2 dx}{||f||^2 = int f(x)^2 dx}. \eqn{\lambda}{lambda} is the scalar regularization parameter that governs the tradeoff between model fit and complexity. By default, \eqn{\lambda}{lambda} is chosen by minimizing the sum of the squared leave-one-out errors, but it can also be specified by the user in the \code{lambda} argument to implement other approaches. 

Under fairly general conditions, the function that minimizes the regularized loss 
within the hypothesis space established by the choice of a (positive semidefinite) kernel function \eqn{k(x_i,x_j)} is of the form

\deqn{f(x_j)= \sum_i^N c_i k(x_i,x_j)}{%
      f(x_j)= sum_i^N c_i k(x_i,x_j)} 

where the kernel function \eqn{k(x_i,x_j)} measures the distance 
between two observations \eqn{x_i} and \eqn{x_j} and \eqn{c_i} is the choice coefficient for each observation \eqn{i}. Let \eqn{K} be the \eqn{N} by \eqn{N} kernel matrix with all pairwise distances \eqn{K_ij=k(x_i,x_j)} and \eqn{c} be the  \eqn{N} by \eqn{1} vector of choice coefficients for all observations then in matrix notation the space is \eqn{y=Kc}. 

Accordingly, the \code{krls} function solves the following minimization problem

\deqn{argmin_{f \in H} \sum_i^n (y - Kc)'(y-Kc)+ \lambda c'Kc}{%
      argmin_{f in H} sum_i^n (y - Kc)'(y-Kc)+ lambda c'Kc}

which is convex in \eqn{c} and solved by \eqn{c=(K +\lambda I)^-1 y}{c=(K +lambda I)^-1 y} where \eqn{I} is the identity matrix. Note that this linear solution provides a flexible fitted response surface that typically reduces misspecification bias because it can learn a wide range of nonlinear and or nonadditive functions of the predictors. 

If \code{vcov=TRUE} is specified, \code{krls} also computes the variance-covariance matrix for the choice coefficients \eqn{c} and fitted values \eqn{y=Kc} based on a variance estimator developed in Hainmueller and Hazlett (2013). Note that both matrices are \var{N} by \var{N} and therefore this results in increased memory and computing time.

By default, \code{krls} uses the Gaussian Kernel (\code{whichkernel = "gaussian"}) given by 

\deqn{k(x_i,x_j)=exp(\frac{-|| x_i - x_j ||^2}{\sigma^2})}{%
      k(x_i,x_j)=exp(-|| x_i - x_j ||^2 / sigma^2)}

where \eqn{||x_i - x_j||} is the Euclidean distance. The kernel bandwidth \eqn{\sigma^2}{sigma^2} is set to \eqn{D}, the number of dimensions, by default, but the user can also specify other values using the \code{sigma} argument to implement other approaches. 


If \code{derivative=TRUE} is specified, \code{krls} also computes the pointwise partial derivatives of the fitted function wrt to each predictor using the estimators developed in Hainmueller and Hazlett (2013). These can be used to examine the marginal effects of each predictor and how they marginal effects vary across the covariate space. Average derivatives are also computed with variances. Note that the \code{derivative=TRUE} option results in increased computing time and is only supported for the Gaussian kernel, i.e. when \code{whichkernel = "gaussian"}. Also \code{derivative=TRUE} requires that \code{vcov=TRUE}. 

If \code{binary=TRUE} is also specified, the function will identify binary predictors and return first differences for these predictors instead of partial derivatives. First differences are computed going from the minimum to the maximum value of each binary predictor. Note that first differences are more appropriate to summarize the effects for binary predictors (see Hainmueller and Hazlett (2013) for details). 

A few other kernels are also implemented, but derivatives are currently not supported for these: "linear": \eqn{k(x_i,x_j)=x_i'x_j}, "poly1", "poly2", "poly3", "poly4" are polynomial kernels based on  \eqn{k(x_i,x_j)=(x_i'x_j +1)^p} where \eqn{p} is the order.

}
\value{
A list object of class \code{krls} with the following elements:
  \item{K}{\var{N} by \var{N} matrix of pairwise kernel distances between observations.}
  \item{coeffs }{\var{N} by 1 vector of choice coefficients \var{c}.}
  \item{Le}{scalar with sum of squared leave-one-out error.}
  \item{fitted}{\var{N} by 1 vector of fitted values.}
  \item{X}{original \var{N} by \var{D} predictor data matrix.}
  \item{y}{original \var{N} by 1 matrix of values of the outcome variable.} 
  \item{sigma}{scalar with value of bandwidth, \eqn{\sigma^2}{sigma^2}, used for the Gaussian kernel.}   
  \item{lambda}{scalar with value of regularization parameter, \eqn{\lambda}{lambda}, used (user specified or based on leave-one-out cross-validation).}  
  \item{R2}{scalar with value of R-squared}
  \item{vcov.c}{\var{N} by \var{N} variance covariance matrix for choice coefficients (\code{NULL} unless \code{vcov=TRUE} is specified).}
  \item{vcov.fitted}{\var{N} by \var{N} variance covariance matrix for fitted values (\code{NULL} unless \code{vcov=TRUE} is specified).}  
  \item{derivatives}{\var{N} by \var{D} matrix of pointwise partial derivatives based on the Gaussian kernel (\code{NULL} unless \code{derivative=TRUE} is specified. If \code{binary=TRUE} is specified, first differences are returned for binary predictors.}
  \item{avgderivatives}{1 by \var{D} matrix of average derivative based on the Gaussian kernel (\code{NULL} unless \code{derivative=TRUE} is specified. If \code{binary=TRUE} is specified, average first differences are returned for binary predictors.}
   \item{var.avgderivatives}{1 by \var{D} matrix of variances for average derivative based on gaussian kernel (\code{NULL} unless \code{derivative=TRUE} is specified. If \code{binary=TRUE} is specified, variances for average first differences are returned for binary predictors.}  
   \item{binaryindicator}{1 by \var{D} matrix that indicates for each predictor if it is treated as binary or not (evaluates to FALSE unless \code{binary=TRUE} is specified and a predictor is recognized binary.}  

  }
\references{

Hainmueller, J. and Hazlett, C. 2013. Kernel Regularized Least Squares: Moving Beyond Linearity and Additivity Without Sacrifcing Interpretability. Working Paper MIT \url{http://dx.doi.org/10.2139/ssrn.2046206}.

Rifkin, R. 2002. Everything Old is New Again: A fresh look at historical approaches in machine learning. Thesis, MIT. September, 2002.

Evgeniou, T., Pontil, M., and Poggio, T. (2000). Regularization networks and support vector
machines. Advances In Computational Mathematics, 13(1):1-50.

Schoelkopf, B., Herbrich, R. and Smola, A.J. (2001) A generalized representer theorem. In 14th Annual
Conference on Computational Learning Theory, pages 416-426.

Kimeldorf, G.S. Wahba, G. 1971. Some results on Tchebycheffian spline functions. Journal of Mathematical Analysis and Applications, 33:82-95.

}
\author{
Jens Hainmueller (MIT) and Chad Hazlett (MIT)
}
\note{
The function requires the storage of a \var{N} by \var{N} kernel matrix and can therefore exceed the memory limits for very large datasets.

Setting \code{derivative=FALSE} and \code{vcov=FALSE} is useful to reduce computing time if pointwise partial derivatives and or variance covariance matrices are not needed.

}


\seealso{
\code{\link{predict.krls}} for fitted values and predictions. \code{\link{summary.krls}} for summary of the fit. \code{\link{plot.krls}} for plots of the fit. 
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
fit <- predict(out,newdata=X,se.fit=TRUE)
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

 \keyword{multivariate, smooth, kernels, machine learning, regression, classification}
