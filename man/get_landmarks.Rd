\name{get_landmarks}
\alias{get_landmarks}
\title{
Extract the landmark coordinates from a Nystrom KRLS fit
}
\description{
Returns the landmark coordinate matrix used by an \code{approx = "nystrom"}
fit. Landmarks are stored internally in the standardized X-space the model
operates in; by default this accessor un-standardizes them so the returned
matrix is in the original \code{X} units and can be passed back through
\code{krls(..., landmarks = ...)} on a comparable dataset without a
standardize-twice bug.
}
\usage{
get_landmarks(fit, scale = c("original", "standardized"))
}
\arguments{
  \item{fit}{
A fitted \code{krls} object built with \code{approx = "nystrom"}.
}
  \item{scale}{
String selecting the coordinate system of the returned matrix.
\code{"original"} (the default) un-standardizes the stored landmarks using
the training data's column means and standard deviations.
\code{"standardized"} returns the internally stored matrix unchanged.
}
}
\value{
An \var{m} by \var{D} numeric matrix of landmark coordinates with column
names inherited from \code{fit$X}.
}
\author{
Jens Hainmueller (Stanford) and Chad Hazlett (UCLA)
}
\seealso{
\code{\link{krls}}
}
\examples{
\dontrun{
fit <- krls(X, y, approx = "nystrom", nystrom_m = 30)
Z   <- get_landmarks(fit)        # original-scale m x D matrix
# Reuse the same landmarks on a sensitivity-check refit with a new y:
fit2 <- krls(X, y_alt, approx = "nystrom", landmarks = Z)
}
}
