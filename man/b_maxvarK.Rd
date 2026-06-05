\name{b_maxvarK}
\alias{b_maxvarK}
\alias{b_maxvarK_nystrom}
\title{
Gaussian kernel bandwidth by maximum off-diagonal var(K)
}
\description{
Pick the Gaussian kernel bandwidth \eqn{\sigma}{sigma} that maximizes
the variance of the off-diagonal entries of the kernel matrix \var{K}.
The idea is that this choice makes the columns of \var{K} most
informative: at very small or very large \eqn{\sigma}{sigma}, the
off-diagonal entries collapse toward 0 or 1 respectively and carry
little discriminating signal.

This is the bandwidth-selection convention used in
\code{kbal::b_maxvarK} and the GPSS package; in KRLS 1.7+ it is also
the default when \code{sigma = NULL} in \code{\link{krls}}.

\code{b_maxvarK_nystrom} is the Nystrom-aware variant: it evaluates
the variance of the entries of the \var{n} by \var{m} cross-kernel
\eqn{C = K(X, Z)}{C = K(X, Z)} between observations and landmarks
instead of forming the full \var{n} by \var{n} kernel.
}
\usage{
b_maxvarK(X_proc, search_lower = 1e-6, search_upper = NULL,
          tol = .Machine$double.eps^0.25)

b_maxvarK_nystrom(X_proc, Z_proc, search_lower = 1e-6,
                  search_upper = NULL,
                  tol = .Machine$double.eps^0.25)
}
\arguments{
  \item{X_proc}{
A numeric matrix in the same kernel-ready form that \code{\link{krls}}
uses internally: continuous columns standardized to \code{sd = 1},
categorical columns one-hot encoded and multiplied by
\code{sqrt(0.5)}. Use \code{krls}'s internal preprocessing if you want
to recover the same matrix it uses.
}
  \item{Z_proc}{
Landmark matrix for the Nystrom variant. Same column structure as
\code{X_proc}.
}
  \item{search_lower, search_upper}{
Lower and upper bounds of the bandwidth search interval passed to
\code{\link[stats]{optimize}}. \code{search_upper} defaults to
\code{max(2 * ncol(X_proc), 2000)}.
}
  \item{tol}{
Tolerance forwarded to \code{\link[stats]{optimize}}. Default is
\code{.Machine$double.eps^0.25} (the \code{optimize} default).
}
}
\value{
A list with components
  \item{sigma}{the selected bandwidth.}
  \item{var_K (or var_C)}{the value of the off-diagonal variance at
the selected bandwidth.}
  \item{search_lower, search_upper}{the bounds used for the search.}
}
\references{
Hazlett, C. (2020). Kernel Balancing: A Flexible Non-Parametric
Weighting Procedure for Estimating Causal Effects. \emph{Statistica
Sinica}.
}
\seealso{
\code{\link{krls}}, \code{\link{gausskernel}}.
}
\keyword{methods}
