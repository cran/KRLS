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
  \item{object}{Fitted KRLS model, i.e. an object of class \code{krls}}
  \item{newdata}{A data frame or matrix with variables values at which to predict the outcome. Number and order of columns in \code{newdata} have to match the corresponding predictors used in the fitted krls model given in \code{object}.
}

\item{se.fit}{logical flag if standard errors should be computed for pointwise predictions.}

\item{\dots}{additional arguments affecting the predictions produced.}

}
\details{

Function produces predicted values, obtained by evaluating the fitted krls function with 
the newdata (ie. the test points). The prediction at a new test point \eqn{x_i} is based on 

\deqn{f(x_i)= \sum_j=1^n c_j K_{x_j}(x_i)}{%
      f(x_i)= sum_j=1^n c_j K_{x_j}(x_i)} 

where \eqn{K} is the kernel matrix and thus \eqn{K_{x_j}}
is a vector whose j-th entry is \eqn{K(x_j,x_i)} (e.g. the distance between the test point \eqn{x_i} and the training point \eqn{x_j}). The training points are passed to the function with the krls fit in \code{object}.

When data are missing in \code{newdata} during prediction, the value of each \eqn{k(x_i,x_j)} is computed by using an adjusted Euclidean distance in the kernel definition. Assume \eqn{x} is \var{D}-dimensional but a given pair of observations \eqn{x_i} and \eqn{x_j} have only \eqn{D' < D} non-missing dimensions in common. The adjusted Euclidean distance computes the sum of squared differences over the \eqn{D'} non-missing dimensions, rescales this sum by \eqn{D/D'}, and takes the square root.   The result corresponds to an assumption that conditional on the observed data, the missing values would not have contributed new information predictive of the outcome.
}
\value{

  \item{fit}{\var{M} by 1 vector of fitted values for \var{M} test points.}
  \item{se.fit}{\var{M} by 1 vector of standard errors for the fitted values for \var{M} test points (\code{NULL} unless \code{se.fit=TRUE} is specified).}
  \item{vcov.fit}{\var{M} by \var{M} variance-covariance matrix for the fitted values for \var{M} test points (\code{NULL} unless \code{se.fit=TRUE} is specified).}
  \item{newdata}{\var{M} by \var{D} data matrix of of \var{M} test points with \var{D} predictors.}
  \item{newdataK}{\var{M} by \var{N} data matrix for pairwise Gauss Kernel distances between \var{M} test points and \var{N} training points from krls model fit in \code{object}.}

}

\author{
Jens Hainmueller (Stanford) and Chad Hazlett (MIT)
}

\seealso{
\code{\link{krls}}
}
\examples{

# make up data
X <- seq(-3,3,.1)
Y <- sin(X) + rnorm(length(X),.1)

# fit krls
krlsout <- krls(y=Y,X=X)

# get in-sample prediction 
predin <- predict(krlsout,newdata=X,se.fit=TRUE)

# get out-of-sample prediction
X2 <- runif(5)
predout <- predict(krlsout,newdata=X2,se.fit=TRUE)

# plot true function and predictions
plot(y=sin(X),x=X,type="l",col="red",ylim=c(-1.8,1.8),lwd=2,ylab="f(X)")
points(y=predin$fit,x=X,col="blue",pch=19)
arrows(y1=predin$fit+2*predin$se.fit,
       y0=predin$fit-2*predin$se.fit,
       x1=X,x0=X,col="blue",length=0)
       
points(y=predout$fit,x=X2,col="green",pch=17)
arrows(y1=predout$fit+2*predout $se.fit,
       y0=predout$fit-2*predout $se.fit,
       x1=X2,x0=X2,col="green",length=0)

legend("bottomright",
       legend=c("true f(x)=sin(X)",
                "KRLS fitted in-sample",
                "KRLS fitted out-of-sample"),
       lty=c(1,NA,NA),pch=c(NA,19,17),
       lwd=c(2,NA,NA),
       col=c("red","blue","green"),
       cex=.8)

}
