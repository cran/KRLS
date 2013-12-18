\name{solveforc}
\alias{solveforc}
\title{
Solve for Choice Coefficients in KRLS 
}
\description{
Internal function that computes choice coefficients for KRLS given a fixed value for lambda (the parameter that governs the tradeoff between model fit and complexity in KRLS).
This function is called internally by \code{\link{krls}}. It would normally not be called by the user directly.
}
\usage{
solveforc(y = NULL, Eigenobject = NULL, lambda = NULL)
}
\arguments{
  \item{y}{
n by 1 matrix of outcomes.
}
  \item{Eigenobject}{
Object from call to \code{\link{eigen}} that contains spectral decomposition of the n by n Kernel matrix.
}
  \item{lambda}{
Positive scalar value for lamnbda parameter.
}
}
\details{
Function relies on fast eigenvalue decomposition method described in  method Rifkin and Lippert (2007).
}
\value{
  \item{coeffs }{n by 1 one matrix of choice coefficients for KRLS model.}
  \item{Le }{n by 1 matrix of errors from leave-one-out validation.}
}
\references{
Rifkin, Ryan M. and Lippert, Ross A. (2007). Notes on Regularized Least Squares. MIT-CSAIL-TR-2007-025. CBCL-268
}
\author{
Jens Hainmueller (Stanford) and Chad Hazlett (MIT)
}

\seealso{
\code{\link{krls}}
}
