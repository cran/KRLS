\name{plot.krls}
\alias{plot.krls}
\title{
Plot method for Kernel-based Regularized Least Squares (KRLS) Model Fits
}
\description{
Produces two types of plots. The first type of plot shows histograms for the pointwise partial derivatives to examine the heterogeneity in the marginal effects of each predictor (\code{which}==1). The second type of plot shows estimates of the conditional expectation functions of \eqn{E[Y|X]} for each predictor (\code{which==2}). For each plot, the predictor of interest varies from its 1st to its 3rd quartile values, while the other predictors are kept at the means (or other values specified in \code{setx}). For binary varibales the \eqn{E[Y|X]} are predicted at the max and the min value of the predictor (instead of the range from the 1st to the 3rd quantile).}
\usage{
 \method{plot}{krls}(x,which=c(1:2),
 main="distributions of pointwise marginal effects",
 setx="mean",ask = prod(par("mfcol")) < nplots,nvalues=50,probs=c(.25,.75),\dots)
}

\arguments{
  \item{x}{An object of class "\code{krls}" that results from call to \code{\link{krls}}.}
  \item{which}{if a subset of the plots is required, specify a subset of the numbers \code{1:2}.}
  \item{main}{main title for histograms of pointwise partial derivatives.}
  \item{setx}{either one of \code{mean} or \code{median} to hold other predictors at their mean or median values for the conditional expectation plots. Alternativley the user can specific a numeric vector with predictor values at which the other predictors should be fixed for the conditional expectation plots. If specifed in this way there must be one value per predictor and the order of the values much match the order of the predictor used in the predictor matrix of the krls fit passed in \code{x}.}
  \item{ask}{logical; if \code{TRUE}, the user is asked before each plot, see \code{\link{par}} (\code{ask=.}).}
  \item{nvalues}{scalar that specifies the number of values at which conditional expectations should be plotted.}
  \item{probs}{vector with numbers between 0 and 1 that specify the quantiles that determine the range for of the predictor values for which the conditional expectation should be plotted. By default we vary each predictor from the 1st quartile to the 3rd quartile value.}
  \item{\dots}{additional arguments to be passed to lower level functions}
}
\details{

Notice that the historgrams for the partial derivatives can only be plotted if the KRLS object was computed with \code{krls(,derivatives=TRUE)}.

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
