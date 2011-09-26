\name{gausskernel}
\alias{gausskernel}
\title{
Gaussian Kernel Distance Computation
}
\description{
This function computes and returns a distance matrix 
computed by using a Gaussian Kernel distance measure to compute the distances between the rows of a data matrix.  
}
\usage{
gausskernel(X = NULL, sigma = NULL)
}

\arguments{
  \item{X}{ an n by k numeric matrix.
}
  \item{sigma}{ a scalar value that specifies the width of the Gaussian kernel sigma^2 (see details).
}
}
\details{

Given two k dimensional vectors x_i and x_j. The Gaussian kernel is defined as: k(x_i,x_j)= exp(-|| x_i - x_j ||^2/sigma^2). 
Where ||x_i - x_j|| is the Euclidean distance between x_i and x_j ||x_i - x_j||=((x_i1-x_j1)^2 + (x_i2-x_j2)^2 + ... + (x_ik-x_jk)^2)^.5. 
Note that the Gaussian kernel is a measure of similarity between the two input vectors 
(it evalues to 1 if the two input arguments are identical, and approaches 0 as the two vectors get very far apart.). 

The function relies on the \code{\link{dist}} function in the stats package for an initial estimate of the euclidean distance.

}
\value{
An n by n numeric distance matrix that contains the pairwise distances between teh rows of X.
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