\name{krls}
\alias{krls}
\title{
Kernel-based Regularized Least Squares (KRLS)
}
\description{

Function implements Kernel-Based Regularized Least Squares (KRLS), a machine learning method described in Hainmueller and Hazlett (2014) that allows users to solve regression and classification problems without manual specification search and strong functional form assumptions. KRLS finds the best fitting function by minimizing a Tikhonov regularization problem with a squared loss, using Gaussian Kernels as radial basis functions. KRLS reduces misspecification bias since it learns the functional form from the data. Yet, it nevertheless allows for interpretability and inference in ways similar to ordinary regression models. In particular, KRLS provides closed-form estimates for the predicted values, variances, and the pointwise partial derivatives that characterize the marginal effects of each independent variable at each data point in the covariate space. The distribution of pointwise marginal effects can be used to examine effect heterogeneity and or interactions. 
}
\usage{
krls(X = NULL, y = NULL, whichkernel = "gaussian", lambda = NULL,
sigma = NULL, derivative = TRUE, binary= TRUE, vcov=TRUE,
print.level = 1,L=NULL,U=NULL,tol=NULL,eigtrunc=NULL,data=NULL,
approx = c("auto", "none", "nystrom"), nystrom_m = NULL,
landmarks = NULL, landmark_method = c("random", "kmeans"),
nystrom_eps = sqrt(.Machine$double.eps),
landmark_seed = NULL,
lambda_method = c("loo", "gcv"),
cat_columns = NULL)
}

\arguments{
  \item{X}{
For the matrix interface (the original): an \var{N} by \var{D} numeric data matrix that contains the values of \var{D} predictor variables for \eqn{i=1,\ldots,N} observations. The matrix may not contain missing values or constants. Note that no intercept is required since the function operates on demeaned data and subtracting the mean of \var{y} is equivalent to including an (unpenalized) intercept into the model.

For the formula interface (added in 1.2-0): a two-sided formula of the form \code{y ~ x1 + x2 + ...} with the response on the left-hand side and the covariates on the right. \code{data} must be supplied. The intercept column is dropped automatically. Both interfaces produce the same fit; the formula form is purely an ergonomic wrapper.
}
  \item{y}{
\var{N} by \var{1} data numeric matrix or vector that contains the values of the response variable for all observations. This vector may not contain missing values. Ignored when \code{X} is a formula (the response is taken from the formula).
}
  \item{data}{
For the formula interface: a \code{data.frame} containing the variables referenced in the formula. Required when \code{X} is a formula; ignored otherwise.
}
  \item{whichkernel}{
String vector that specifies which kernel should be used. Must be one of \code{gaussian}, \code{linear}, \code{poly1}, \code{poly2}, \code{poly3}, or \code{poly4} (see details). Default is \code{gaussian}.
}
  \item{lambda}{
A positive scalar that specifies the \eqn{\lambda}{lambda} parameter for the regularizer (see details). It governs the tradeoff between model fit and complexity. By default, this parameter is chosen by minimizing the sum of the squared leave-one-out errors.
}
  \item{sigma}{
A positive scalar that specifies the bandwidth of the Gaussian kernel (see \code{\link{gausskernel}} for details).

\strong{Default changed in 1.7-0.} When \code{sigma = NULL}, \code{krls()} now selects the bandwidth that maximizes the off-diagonal variance of the kernel matrix \var{K} (see \code{\link{b_maxvarK}}), matching the convention used by \code{kbal::b_maxvarK} and \code{gpss}. The pre-1.7 default was \code{sigma = ncol(X_processed)}; pin that value explicitly to reproduce older fits. A one-line \code{message()} is emitted on the first call where the new default kicks in.
}
  \item{derivative}{
Logical that specifies whether pointwise partial derivatives should be computed. Currently, derivatives are only implemented for the Gaussian Kernel.
}
  \item{binary}{
Logical that specifies whether first-differences instead of pointwise partial derivatives should be computed for binary predictors. Ignored unless \code{derivative=TRUE}. 
}
\item{vcov}{
Logical that specifies whether variance-covariance matrix for the choice coefficients \var{c} and fitted values should be computed. Note that \code{derivative=TRUE} requires that \code{vcov=TRUE} for the exact path. Under \code{approx = "nystrom"}, \code{vcov=TRUE} computes conditional approximate variances for the landmark coefficients and average derivatives without storing the full fitted-value covariance matrix.
}

\item{print.level}{
Positive integer that determines the level of printing. Set to 0 for no printing and 2 for more printing.
}

\item{L}{
Non-negative scalar that determines the lower bound of the search window for the leave-one-out optimization to find \eqn{\lambda}{lambda}. Default is \code{NULL} which means that the lower bound is found by using an algorithm outlined in \code{\link{lambdasearch}}.
}

\item{U}{
Positive scalar that determines the upper bound of the search window for the leave-one-out optimization to find \eqn{\lambda}{lambda}. Default is \code{NULL} which means that the upper bound is found by using an algorithm outlined in \code{\link{lambdasearch}}.
}

\item{tol}{
Positive scalar that determines the tolerance used in the optimization routine used to find \eqn{\lambda}{lambda}. Default is \code{NULL} which means that convergence is achieved when the difference in the sum of squared leave-one-out errors between the \var{i} and the \var{i+1} iteration is less than \var{N * 10^-3}.
}

\item{eigtrunc}{
Positive scalar that determines how much eignvalues should be trunacted for finding the upper bound of the search window in the algorithm outlined in \code{\link{lambdasearch}}. If \code{eigtrunc} is set to \var{10^-6} this means that we keep only eigenvalues that are \var{10^-6} as large as the first. Default is \code{eigtrunc=NULL} which means no truncation is used.
}

\item{approx}{
String selecting the approximation regime.

\code{"auto"} (the default) uses the exact \var{N} by \var{N} kernel
when the sample size does not exceed the would-be landmark count
(\code{nystrom_m}, default \code{min(500, N)}), and switches to the
Nystrom approximation otherwise. A one-line message is emitted when
the switch happens, so the choice is visible at the call site.
Conditions that are unsupported under Nystrom (non-Gaussian kernel,
non-null \code{eigtrunc}, \code{derivative = TRUE} with
\code{vcov = FALSE}) fall through to the exact path silently.

\code{"none"} forces the exact path regardless of sample size; this
is the historical KRLS behavior and gives identical fits to KRLS
\eqn{\le}{<=} 1.5-2.

\code{"nystrom"} forces the low-rank Nystrom approximation: an
explicit \var{m}-dimensional feature map replaces the full kernel and
the ridge problem is solved in \var{m}-space. Time becomes
\var{O(N m^2 + m^3)} and memory \var{O(N m)}, so KRLS is feasible at
\var{N} well beyond the exact path's ~5000-row ceiling. When
\code{vcov = TRUE}, coefficient, prediction, and average-derivative
standard errors are computed as conditional approximate Nystrom
standard errors -- conditional on the selected landmarks, fixed
\eqn{\lambda}{lambda}, and the low-rank feature approximation. The
full \var{N} by \var{N} fitted-value covariance matrix is not stored.
Only the Gaussian kernel is supported.
}

\item{nystrom_m}{
Integer number of landmark points to use when \code{approx = "nystrom"}.
Ignored otherwise. Default is \code{min(500, N)}, which keeps every
row as a landmark for small samples and caps at 500 for larger ones.
Larger \code{nystrom_m} approaches exact KRLS at a quadratic cost;
smaller trades fit quality for runtime/memory. The internal
approximation is \eqn{K \approx C\,W_{reg}^{-1}\,C^\top}{K = C
W_reg^(-1) C'} where \eqn{C = K(X, Z)}{C = K(X, Z)},
\eqn{W = K(Z, Z)}{W = K(Z, Z)}, and \eqn{W_{reg}}{W_reg} is \eqn{W}
with eigenvalues floored at \code{nystrom_eps * max(eigenvalues of W)}.
}

\item{landmarks}{
Optional landmark specification for \code{approx = "nystrom"}. One of:
\code{NULL} (auto-select \code{nystrom_m} landmarks via
\code{landmark_method}); an integer vector of row indices into \code{X};
or an \var{m} by \var{D} numeric matrix of landmark coordinates in the
original \code{X} units. If indices are supplied, the resolved standardized
landmark matrix and the original index vector are both stored on the fit
object; if a matrix is supplied, only the matrix is stored.
}

\item{landmark_method}{
String selecting the auto-selection rule for landmarks when
\code{landmarks = NULL}. \code{"random"} (the default) samples uniformly
without replacement and is bit-reproducible under \code{set.seed()}.
\code{"kmeans"} runs \code{stats::kmeans()} on the standardized
predictors and uses cluster centers as landmarks -- typically more
robust on imbalanced predictor distributions, but only approximately
reproducible across R versions because the algorithm's initialization
is implementation-dependent.
}

\item{nystrom_eps}{
Positive scalar controlling the relative-ridge stabilization of the
landmark kernel \var{W}. Eigenvalues of \var{W} below
\code{nystrom_eps * max(eigenvalues of W)} are floored at that threshold
before the feature map is built; this stabilizes the
\eqn{W^{-1/2}}{W^(-1/2)} operation when landmarks are near-collinear in
feature space. Documented as numerical stabilization, not a modeling
parameter. Default is \code{sqrt(.Machine$double.eps)}.
}

\item{landmark_seed}{
Optional integer seed used \emph{only} for landmark selection when
\code{landmarks = NULL} and the chosen \code{landmark_method} is
stochastic (random sampling or k-means). If supplied, the global
\code{.Random.seed} is saved before landmark selection and restored
on exit, so the caller's downstream RNG state is unaffected. This
makes Nystrom fits bit-reproducible without forcing the user to
remember a \code{set.seed()} call in their own code and without the
silent-seed-reset failure mode that can plague simulations. Default
is \code{NULL}, in which case landmark selection draws from
whichever RNG state the caller has in effect.
}

\item{cat_columns}{
Optional categorical-column declaration (new in 1.7-0). One of:
\code{NULL} (the default; no one-hot expansion, but a single warning
is emitted if columns of \code{X} look categorical -- factor /
character / logical / numeric with few unique values -- and were not
declared); a character vector of column names in \code{X}; or an
integer vector of column indices.

Declared columns are one-hot encoded with all levels (no reference
cell) and multiplied by \code{sqrt(0.5)}, matching the convention
used in \code{kbal::one_hot} and \code{gpss}. The factor of
\code{sqrt(0.5)} compensates for the doubled squared-Euclidean
distance from carrying both the "in level k" and "not in level k"
indicators, so two observations differing in one category contribute
1 (not 2) to \eqn{||x_i - x_j||^2}{||x_i - x_j||^2} -- equivalent to
a Hamming distance. Continuous columns are still standardized to
\code{sd = 1}.

There is no autodetection. Pass \code{cat_columns = integer(0)} (or
\code{character(0)}) to acknowledge that you have looked at \var{X}
and nothing is categorical; this silences the heuristic warning.
}

\item{lambda_method}{
String selecting the objective minimized to choose
\eqn{\lambda}{lambda} when \code{lambda} is not supplied.
\code{"loo"} (the default) uses the closed-form leave-one-out sum of
squared errors -- the historical KRLS criterion. \code{"gcv"} uses
generalized cross-validation,
\eqn{RSS(\lambda) / (1 - \mathrm{tr}(S(\lambda))/n)^2}{RSS(lambda) /
(1 - tr(S(lambda))/n)^2}, computed in closed form from the kernel
eigendecomposition (exact path) or the cached SVD of \eqn{\Phi}{Phi}
(Nystrom path). Both methods are evaluated at the same per-iteration
cost; GCV is often preferred when one or two observations have
disproportionate leverage. Ignored when \code{lambda} is supplied
directly.
}

}
\details{

\code{krls} implements the Kernel-based Regularized Least Squares (KRLS) estimator as described in Hainmueller and Hazlett (2014). Please consult this reference for any details.

Kernel-based Regularized Least Squares (KRLS) arises as a Tikhonov minimization problem with a squared loss. Assume we have data of the from \eqn{y_i,\,x_i}{y_i, x_i} where \var{i} indexes observations, \eqn{y_i \in R}{y_i in R} is the outcome and \eqn{x_i \in R^D}{x_i in R^D} is a \var{D}-dimensional vector of predictor values. Then KRLS searches over a space of functions \eqn{H} and chooses the best fitting function \eqn{f} according to the rule: 

\deqn{argmin_{f \in H} \sum_i^N (y_i - f(x_i))^2 + \lambda ||f||_{H^2}}{%
      argmin_{f in H} sum_i^N (y_i - f(x_i))^2 + lambda || f ||_H^2} 

where \eqn{(y_i - f(x_i))^2} is a loss function that computes how `wrong' the function 
is at each observation \var{i} and \eqn{|| f ||_{H^2}}{|| f ||_H^2} is the regularizer that measures the complexity of the function according to the \eqn{L_2} norm \eqn{||f||^2 = \int f(x)^2 dx}{||f||^2 = int f(x)^2 dx}. \eqn{\lambda}{lambda} is the scalar regularization parameter that governs the tradeoff between model fit and complexity. By default, \eqn{\lambda}{lambda} is chosen by minimizing the sum of the squared leave-one-out errors, but it can also be specified by the user in the \code{lambda} argument to implement other approaches. 

Under fairly general conditions, the function that minimizes the regularized loss 
within the hypothesis space established by the choice of a (positive semidefinite) kernel function \eqn{k(x_i,x_j)} is of the form

\deqn{f(x_j)= \sum_i^N c_i k(x_i,x_j)}{%
      f(x_j)= sum_i^N c_i k(x_i,x_j)} 

where the kernel function \eqn{k(x_i,x_j)} measures the distance 
between two observations \eqn{x_i} and \eqn{x_j} and \eqn{c_i} is the choice coefficient for each observation \eqn{i}. Let \eqn{K} be the \eqn{N} by \eqn{N} kernel matrix with all pairwise distances \eqn{K_ij=k(x_i,x_j)} and \eqn{c} be the  \eqn{N} by \eqn{1} vector of choice coefficients for all observations then in matrix notation the space is \eqn{y=Kc}. 

Accordingly, the \code{krls} function solves the following minimization problem

\deqn{argmin_{f \in H} \sum_i^n (y - Kc)'(y-Kc)+ \lambda c'Kc}{%
      argmin_{f in H} sum_i^n (y - Kc)'(y-Kc)+ lambda c'Kc}

which is convex in \eqn{c} and solved by \eqn{c=(K +\lambda I)^-1 y}{c=(K +lambda I)^-1 y} where \eqn{I} is the identity matrix. Note that this linear solution provides a flexible fitted response surface that typically reduces misspecification bias because it can learn a wide range of nonlinear and or nonadditive functions of the predictors. 

If \code{vcov=TRUE} is specified, exact KRLS also computes the variance-covariance matrix for the choice coefficients \eqn{c} and fitted values \eqn{y=Kc} based on a variance estimator developed in Hainmueller and Hazlett (2014). Note that both exact matrices are \var{N} by \var{N} and therefore this results in increased memory and computing time. Under \code{approx = "nystrom"}, \code{vcov=TRUE} computes conditional approximate variances in the \var{m}-landmark feature space and does not store the full fitted-value covariance matrix.

By default, \code{krls} uses the Gaussian Kernel (\code{whichkernel = "gaussian"}) given by 

\deqn{k(x_i,x_j)=exp(\frac{-|| x_i - x_j ||^2}{\sigma^2})}{%
      k(x_i,x_j)=exp(-|| x_i - x_j ||^2 / sigma^2)}

where \eqn{||x_i - x_j||} is the Euclidean distance. The kernel bandwidth \eqn{\sigma^2}{sigma^2} is set to \eqn{D}, the number of dimensions, by default, but the user can also specify other values using the \code{sigma} argument to implement other approaches. 


If \code{derivative=TRUE} is specified, \code{krls} also computes the pointwise partial derivatives of the fitted function wrt to each predictor using the estimators developed in Hainmueller and Hazlett (2014). These can be used to examine the marginal effects of each predictor and how the marginal effects vary across the covariate space. Average derivatives are also computed with variances when \code{vcov=TRUE}. Note that the \code{derivative=TRUE} option results in increased computing time and is only supported for the Gaussian kernel, i.e. when \code{whichkernel = "gaussian"}. For the exact path, \code{derivative=TRUE} requires \code{vcov=TRUE}; under \code{approx = "nystrom"}, \code{derivative=TRUE, vcov=FALSE} returns derivative point estimates without standard errors.

If \code{binary=TRUE} is also specified, the function will identify binary predictors and return first differences for these predictors instead of partial derivatives. First differences are computed going from the minimum to the maximum value of each binary predictor. Note that first differences are more appropriate to summarize the effects for binary predictors (see Hainmueller and Hazlett (2014) for details). 

A few other kernels are also implemented, but derivatives are currently not supported for these: "linear": \eqn{k(x_i,x_j)=x_i'x_j}, "poly1", "poly2", "poly3", "poly4" are polynomial kernels based on  \eqn{k(x_i,x_j)=(x_i'x_j +1)^p} where \eqn{p} is the order.

}
\value{
A list object of class \code{krls} with the following elements:
  \item{K}{\var{N} by \var{N} matrix of pairwise kernel distances between observations.}
  \item{coeffs }{\var{N} by 1 vector of choice coefficients \var{c}.}
  \item{Le}{scalar with sum of squared leave-one-out errors.}
  \item{fitted}{\var{N} by 1 vector of fitted values.}
  \item{X}{original \var{N} by \var{D} predictor data matrix.}
  \item{y}{original \var{N} by 1 matrix of values of the outcome variable.} 
  \item{sigma}{scalar with value of bandwidth, \eqn{\sigma^2}{sigma^2}, used for the Gaussian kernel.}   
  \item{lambda}{scalar with value of regularization parameter, \eqn{\lambda}{lambda}, used (user specified or based on leave-one-out cross-validation).}  
  \item{R2}{scalar with value of R-squared}
  \item{vcov.c}{Variance covariance matrix for choice coefficients (\var{N} by \var{N} for exact KRLS; \var{m} by \var{m} conditional approximate covariance under \code{approx = "nystrom"}; \code{NULL} unless \code{vcov=TRUE} is specified).}
  \item{vcov.fitted}{\var{N} by \var{N} variance covariance matrix for fitted values under exact KRLS (\code{NULL} unless \code{vcov=TRUE} is specified, and always \code{NULL} under \code{approx = "nystrom"} to preserve the memory benefit).}
  \item{derivatives}{\var{N} by \var{D} matrix of pointwise partial derivatives based on the Gaussian kernel (\code{NULL} unless \code{derivative=TRUE} is specified. If \code{binary=TRUE} is specified, first differences are returned for binary predictors.}
  \item{avgderivatives}{1 by \var{D} matrix of average derivative based on the Gaussian kernel (\code{NULL} unless \code{derivative=TRUE} is specified. If \code{binary=TRUE} is specified, average first differences are returned for binary predictors.}
   \item{var.avgderivatives}{1 by \var{D} matrix of variances for average derivative based on gaussian kernel (\code{NULL} unless \code{derivative=TRUE} is specified. If \code{binary=TRUE} is specified, variances for average first differences are returned for binary predictors.}  
   \item{binaryindicator}{1 by \var{D} matrix that indicates for each predictor if it is treated as binary or not (evaluates to FALSE unless \code{binary=TRUE} is specified and a predictor is recognized binary.}  

  }
\references{

Jeremy Ferwerda, Jens Hainmueller, Chad J. Hazlett (2017). Kernel-Based Regularized Least Squares in {R} ({KRLS}) and {Stata} ({krls}).
 Journal of Statistical Software, 79(3), 1-26. doi:10.18637/jss.v079.i03
  
Hainmueller, J. and Hazlett, C. (2014). Kernel Regularized Least Squares: Reducing Misspecification Bias with a Flexible and Interpretable Machine Learning Approach. Political Analysis, 22(2)

Rifkin, R. 2002. Everything Old is New Again: A fresh look at historical approaches in machine learning. Thesis, MIT. September, 2002.

Evgeniou, T., Pontil, M., and Poggio, T. (2000). Regularization networks and support vector machines. Advances In Computational Mathematics, 13(1):1-50.

Schoelkopf, B., Herbrich, R. and Smola, A.J. (2001) A generalized representer theorem. In 14th Annual Conference on Computational Learning Theory, pages 416-426.

Kimeldorf, G.S. Wahba, G. 1971. Some results on Tchebycheffian spline functions. Journal of Mathematical Analysis and Applications, 33:82-95.

}
\author{
Jens Hainmueller (Stanford) and Chad Hazlett (UCLA)
}
\note{
The function requires the storage of a \var{N} by \var{N} kernel matrix and can therefore exceed the memory limits for very large datasets.

Setting \code{derivative=FALSE} and \code{vcov=FALSE} is useful to reduce computing time if pointwise partial derivatives and or variance covariance matrices are not needed.

}


\seealso{
\code{\link{predict.krls}} for fitted values and predictions. \code{\link{summary.krls}} for summary of the fit. \code{\link{plot.krls}} for plots of the fit. 
}
\examples{

# Linear example
# set up data
N <- 200
x1 <- rnorm(N)
x2 <- rbinom(N,size=1,prob=.2)
y <- x1 + .5*x2 + rnorm(N,0,.15)
X <- cbind(x1,x2)
# fit model
krlsout <- krls(X=X,y=y)
# summarize marginal effects and contribution of each variable
summary(krlsout)
# plot marginal effects and conditional expectation plots
plot(krlsout)


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

## 2D example:
# predictor data
X <- matrix(seq(-3,3,.1))
# true function
Ytrue <- sin(X)
# add noise 
Y     <- sin(X) + rnorm(length(X),sd=.3)
# approximate function using KRLS
out <- krls(y=Y,X=X)
# get fitted values and ses
fit <- predict(out,newdata=X,se.fit=TRUE)
# results
par(mfrow=c(2,1))
plot(y=Ytrue,x=X,type="l",col="red",ylim=c(-1.2,1.2),lwd=2,main="f(x)")
points(y=fit$fit,X,col="blue",pch=19)
arrows(y1=fit$fit+1.96*fit$se.fit,
       y0=fit$fit-1.96*fit$se.fit,
       x1=X,x0=X,col="blue",length=0)
legend("bottomright",legend=c("true f(x)=sin(x)","KRLS fitted f(x)"),
       lty=c(1,NA),pch=c(NA,19),lwd=c(2,NA),col=c("red","blue"),cex=.8)

plot(y=cos(X),x=X,type="l",col="red",ylim=c(-1.2,1.2),lwd=2,main="df(x)/dx")
points(y=out$derivatives,X,col="blue",pch=19)

legend("bottomright",legend=c("true df(x)/dx=cos(x)","KRLS fitted df(x)/dx"),
       lty=c(1,NA),pch=c(NA,19),lwd=c(2,NA),col=c("red","blue"),,cex=.8)

## 3D example
# plot true function
par(mfrow=c(1,2))
f<-function(x1,x2){ sin(x1)*cos(x2)}
x1 <- x2 <-seq(0,2*pi,.2)
z   <-outer(x1,x2,f)
persp(x1, x2, z,theta=30,main="true f(x1,x2)=sin(x1)cos(x2)")
# approximate function with KRLS
# data and outcomes
X <- cbind(sample(x1,200,replace=TRUE),sample(x2,200,replace=TRUE))
y   <- f(X[,1],X[,2])+ runif(nrow(X))
# fit surface
krlsout <- krls(X=X,y=y)
# plot fitted surface
ff  <- function(x1i,x2i,krlsout){predict(object=krlsout,newdata=cbind(x1i,x2i))$fit}
z   <- outer(x1,x2,ff,krlsout=krlsout)
persp(x1, x2, z,theta=30,main="KRLS fitted f(x1,x2)")

}

 \keyword{multivariate}
\keyword{smooth}
\keyword{regression}
\keyword{classification}
