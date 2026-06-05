## Accessor for the landmark coordinates on a Nystrom fit. Internally
## the fit stores landmarks in standardized X-space; this function
## un-standardizes them on request so users can pass them back through
## `krls(..., landmarks = ...)` without hitting the standardize-twice
## bug.

get_landmarks <-
function(fit, scale = c("original", "standardized"))
  {
    if (!inherits(fit, "krls"))
      stop("fit is not of class 'krls'")
    if (is.null(fit$landmarks))
      stop("fit was not built with approx = 'nystrom'; ",
           "no landmarks to return")
    scale <- match.arg(scale)
    if (scale == "standardized")
      return(fit$landmarks)
    # Un-standardize using the training X's column means / SDs (the
    # same centers/scales used at fit time).
    Xmeans <- colMeans(fit$X)
    Xsds   <- apply(fit$X, 2, sd)
    out    <- sweep(sweep(fit$landmarks, 2, Xsds, `*`), 2, Xmeans, `+`)
    colnames(out) <- colnames(fit$X)
    out
  }
