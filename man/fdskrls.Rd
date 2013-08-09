\name{fdskrls}
\alias{fdskrls}
\title{
Compute first differences with KRLS
}
\description{
Internal function that is called by  \code{\link{krls}} to computes first differences for binary predictors in the X matrix. It would normally not be called by the user directly.}
\usage{
fdskrls(object,\dots)
}
\arguments{
  \item{object}{Object from call to \code{\link{krls}}.
}
\item{\dots}{additional arguments to be passed to lower level functions}
}

\value{
A object of class \code{krls} where the derivatives, average derivatives, and the varinaces of the average derivatives are 
replaced with the first differences for binary predictors. The binaryindicator is also updated and set to TRUE for binary predictors.
}
\author{
Jens Hainmueller (MIT) and Chad Hazlett (MIT)
}

\seealso{
\code{\link{krls}}
}
