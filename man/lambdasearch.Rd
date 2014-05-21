\name{lambdasearch}
\alias{lambdasearch}
\title{
Leave-one-out optimization to find \eqn{\lambda}{lambda}
}
\description{
Function conducts leave-one-out optimization to find \eqn{\lambda}{lambda} using a golden search search with caching. This function is called internally by \code{\link{krls}}. It would normally not be called by the user directly.

}
\usage{
lambdasearch(L=NULL,
           U=NULL,
           y=NULL,
           Eigenobject=NULL,
           tol=NULL,
           noisy=FALSE,
           eigtrunc=NULL)
}
\arguments{
  \item{L}{Non-negative scalar that determines the lower bound of the search window. Default is \code{NULL} which means that the lower bound is found using an algorithm (see details).}
\item{U}{Positive scalar that determines the upper bound of the search window. Default is \code{NULL} which means that the upper bound is found using an algorithm (see details).}
\item{y}{\var{N} by \var{1} matrix of outcomes.}
\item{Eigenobject}{List that contains the eigenvalues and eigenvectors of the kernel matrix \var{K}.}
\item{tol}{Positive scalar that determines the tolerance used in the optimization routine used to find \eqn{\lambda}{lambda}. Default is \code{NULL} which means that convergence is achieved when the difference in the sum of squared leave-one-out errors between the \code{i} and the \code{i+1} iteration is less than \code{N * 10^-3}.}
\item{noisy}{If \code{TRUE}, the function will print details of the golden section search.}
\item{eigtrunc}{
Positive scalar value that determines truncation of eigenvalues for lamnda search window. See \code{\link{krls}} for details. Default is \code{NULL} which means no truncation.
}
}
\details{

By default, upper bound is found as follows:
Set j to n, decrease by one until the following is longer true:
sum(EigenValues / (EigenValues + j)) < 1.

By default, upper bound is found as follows:
Get the position, q, of the eigenvalue that is closest to max(Eigenvalue)/1000.
Set j to 0, increase in steps of 0.05 until the below is longer true: 
sum(EigenValues / (EigenValues + j)) > q.

}
\value{
A scalar that contains the \eqn{\lambda}{lambda} that minimizes the sum of squared leave-one-out errors.
}

\author{
Jens Hainmueller (Stanford) and Chad Hazlett (MIT)
}

\seealso{
\code{\link{krls}}
}
