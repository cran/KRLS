\name{summary.krls}
\alias{summary.krls}
\title{
Summary method for Kernel-based Regularized Least Squares (KRLS) Model Fits
}
\description{
Summarizes average partial derivatives (i.e. marginal effects) and the distribution of the partial derivatives for each predictor. For binary predictors, the marginal effects are the first differences if \code{krls(,derivatives=TRUE,binary=TRUE)} was specified.
}
\usage{

 \method{summary}{krls}(object, probs=c(.25,.5,.75),\dots)
}
\arguments{
  \item{object}{Fitted krls model, i.e. an object of class krls}
  
\item{probs}{numeric vector with numbers between 0 and 1 that specify the quantiles of the pointwise marginal effects for the summary (see the \code{\link{quantile}} function for details).}

\item{\dots}{additional arguments to be passed to lower level functions}

}
\details{

Notice that the partial derivatives can only be summarized if the krls object was computed with \code{krls(,derivatives=TRUE)}.

}
\value{

  \item{coefficients}{matrix with average partial derivates and or first differences (point estimates, standart errors, t-values, p-values).}
  \item{qcoefficients}{matrix with 1st, 2nd, and 3rd quatriles of distribution of pointwise marinal effects.}

}

\author{
Jens Hainmueller (Stanford) and Chad Hazlett (MIT)
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
