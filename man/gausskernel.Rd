\name{gausskernel}
\alias{gausskernel}
\title{
Gaussian Kernel Distance Computation
}
\description{
Given a \var{N} by \var{D} numeric data matrix, this function computes the \var{N} by \var{N} distance matrix with the pairwise distances between the rows of the data matrix as measured by a Gaussian Kernel. 
}
\usage{
gausskernel(X = NULL, sigma = NULL)
}

\arguments{
  \item{X}{\var{N} by \var{N} numeric data matrix.
}
  \item{sigma}{Positive scalar that specifies the bandwidth of the Gaussian kernel (see details).
}
}
\details{

Given two \var{D} dimensional vectors \eqn{x_i} and \eqn{x_j}. The Gaussian kernel is defined as

\deqn{k(x_i,x_j)=exp(\frac{-|| x_i - x_j ||^2}{\sigma^2})}{%
      k(x_i,x_j)=exp(-|| x_i - x_j ||^2 / sigma^2)}

where \eqn{||x_i - x_j||} is the Euclidean distance given by

\deqn{||x_i - x_j||=((x_i1-x_j1)^2 + (x_i2-x_j2)^2 + ... + (x_iD-x_jD)^2)^.5}

and \eqn{\sigma^2}{sigma^2} is the bandwidth of the kernel.

Note that the Gaussian kernel is a measure of similarity between \eqn{x_i} and \eqn{x_j}. It evalues to 1 if the \eqn{x_i} and \eqn{x_j} are identical, and approaches 0 as \eqn{x_i} and \eqn{x_j} move further apart. 

The function relies on the \code{\link{dist}} function in the stats package for an initial estimate of the euclidean distance.

}
\value{
An \var{N} by \var{N} numeric distance matrix that contains the pairwise distances between the rows in \var{X}.
}

\author{
Jens Hainmueller (MIT) and Chad Hazlett (MIT)
}

\seealso{
\code{\link{dist}} function in the stats package.
}
\examples{
X <- matrix(rnorm(6),ncol=2)
gausskernel(X=X,sigma=1)
}