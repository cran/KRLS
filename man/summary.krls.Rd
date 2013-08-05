\name{summary.krls}
\alias{summary.krls}
\title{
Summary method for Kernel-based Regularized Least Squares (KRLS) Model Fits
}
\description{
Summarizes average partial derivatives (i.e. marginal effect) and the distribution of the partial derivatives for each predictor. For binary predictors, the first difference is computed instead of the average partial derivatives.
}
\usage{

 \method{summary}{krls}(object, binary=TRUE,probs=seq(.25,.75,.25),\dots)
}
\arguments{
  \item{object}{Fitted krls model, i.e. an object of class krls}
  
\item{binary}{logical flag if first differences, instead of average partial derivatives, should be computed for binary predictors.}

\item{probs}{numeric vector with numbers between 0 and 1 that specify the quantiles of the pointwise marginal effects for the summary.}

\item{\dots}{additional arguments to be passed to lower level functions}

}
\details{

Notice that the partial derivatives can only be summarized if the krls object was computed with \code{krls(,derivatives=TRUE)}.

}
\value{

  \item{coefficients}{matrix with average partial derivates and or first differences (point estimates, standart errors, t-values, p-values).}
  \item{qcoefficients}{matrix with 1st, 2nd, and 3rd quatriles of distribution of pointwise marinal effects.}
\item{dydx}{matrix with pointwise partial derivatives where, if binary is set to TRUE, for the binary varibales the pointwise first diferences are reported instead of the partial derivatives.}

}

\author{
Jens Hainmueller (MIT) and Chad Hazlett (MIT)
}

\seealso{
\code{\link{krls}}
}
\examples{
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

}
