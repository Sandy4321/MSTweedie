\name{cv.MSTweedie}
\alias{cv.MSTweedie}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-Validation on the Multi-sourse sparse Tweedie model
}
\description{
This function performs k-fold cross-validation (CV) of the multi-source sparse Tweedie model (MSTweedie) partly based on the \code{\link{glmnet::cv}} function. The Tweedie model deviance is the statistical criterion for model selection.
}
\usage{
cv.MSTweedie(x, y, w, source, rho, nlambda = 100, lambda,
      lambda.min = 0.001, nfolds = 10, kktstop = F, foldid,
      adaptive = 0, x.normalize = TRUE, reg = c("L2", "Linf"),
      eps = 5e-04, sr = TRUE, maxit = 10000,
      pf = rep(1, nvars), alpha = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data frame containing the predictors, the responses (identifying the sources either by different columns in the simultaneous case or via an additionnal index column) and, optionnaly, the observation weigths.
}
  \item{y}{
Either (1) a single integer identifying the column of \code{x} containing the response (requires \code{source} to be specified), (2) a vector of integers indentifying which columns of \code{x} are the responses (simultaneous case).
}
  \item{w}{
A single integer identifying the column of \code{x} containing the observation weights. If this argument is missing, equal weight is assumed.
}
  \item{source}{
When \code{y} is a single integer, this arguments identifies the column of \code{x} which indexes the different sources. Disregard is \code{y} is a vector or list of vectors.
}
  \item{rho}{
Power used for the mean-variance relation of the Tweedie distribution. Possible range is [1,2], default is 1.5.
}
  \item{nlambda}{
The length of the regularization path. Disregarded if \code{lambda} is specified, default if 100.
}
  \item{lambda.min}{
The fraction of the first regularization parameter (which is computed to be the smallest such that no predictors are included) defining the last regularization parameter. Disregarded if \code{lambda} is specified; possible range is (0,1), default is 1e-3.
}
  \item{lambda}{
(Optional) User specified sequence of regularization parameter with positive values. When omitted, the sequence is computed starting from the smallest value excluding all predictors from the model and decreasing to a fraction \code{lambda.min} of that starting value by logarithmic decreaments.
}
  \item{nfolds}{
Number of folds in the k-folds CV. Default is 10.
}
  \item{kktstop}{
Logical flag for using the KKT conditions to stop the fit before the end of the regularization parameter sequence. Default is \code{FALSE}. Only applies to the preliminary fit.
}
  \item{foldid}{
(Optional) An list of vector of values between 1 and \code{nfolds} identifying what fold each observation is in. If supplied, \code{nfolds} can be missing. If missing, the folds are constructed randomly.
}
  \item{adaptive}{
Exponent of the adaptive penalty weights; suggested value is 1 (See reference for details.) When the argument is 0, no adaptation is performed.
}
  \item{x.normalize}{
Logical flag for stadardization of the predictors prior to fitting the model. If \code{TRUE}, each predictors in each source is centered to zero and scaled to variance 1. After the fit of the model, the coefficients are returned on the original scale. Default is \code{FALSE}.
}
  \item{reg}{
Either \code{"Linf"} for using \eqn{L_\infty}-regularization in the fit or \code{"L2"} for the \eqn{L_2}-regularization. Default is \code{"Linf"}.
}
  \item{eps}{
Convergence threshold. Default is 5e-4.
}
  \item{sr}{
Logical flag for using the strong rule in the fit. Default is \code{TRUE}.
}
  \item{maxit}{
Maximum number of inner-loop iterations. Default is 10,000.
}
  \item{pf}{
Penalty weights in the penalty term by feature. Mostly used internally when the Adaptive Lasso is used in cross-validation. Expects a vector of length \code{nvars}, default is 1.
}
  \item{alpha}{
Parameter controlling the balance between across-feature and within-feature sparsity in the penalty term
\deqn{(1-\alpha)||\beta||_q +\alpha||\beta||_1.} Possible range is [0,1], default is 0.
}
  \item{\dots}{
   Further arguments to be passed to \code{MSTweedie}.
}
}
\details{
The function runs the \code{MSTweedie} function on the first dataset to get the sequence of regularization paramter and the coefficient estimates. Then, it performs CV along the solution path and computes the out-of-sample Tweedie deviance based on the predicted responses. For each value of the regularization parameter, the error is averaged over the number of folds used and the standard error is computed.
}
\value{
\item{lambda}{A vector containing the sequence of regularization parameters.}
\item{cvm}{A vector containing the Tweedie deviance mean across the folds along the solution path.}
\item{cvsd}{A vector containing the Tweedie deviance standard error across the folds along the solution path.}
\item{cvupper}{\code{cvm+cvsd}.}
\item{cvlo}{\code{cvm-cvsd}.}
\item{reg}{The type of regularization used in the algorithm.}
\item{alpha}{The value of the argument \code{alpha} used.}
\item{name}{\code{"Tweedie Deviance"}, the loss function used.}
\item{MSTweedie.fit}{The \code{MSTweedie} object (fitted on the whole dataset.)}
\item{time}{Computing time.}
\item{lambda.min}{Regularization parameter at minimum CV error.}
\item{lambda.1se}{Largest regularization parameter within one standard error of the minimum.}
}
\references{
Fontaine, S., Yang, Y., Fan, B., Qian, W. and Gu, Y. (2018). "A Unified Approach to Sparse Tweedie Model
with Big Data Applications to Multi-Source
Insurance Claim Data Analysis," to be submitted.

Friedman, J., Hastie, T., Simon, N., Qian, J. and Tibshirani, R. (2017). "\code{glmnet}: Lasso and Elastic-Net Regularized Generalized Linear Models." A vignette for \code{R} package \code{glmnet}. Available from \url{https://cran.r-project.org/web/packages/glmnet}.

\code{glmnet} package.
}

\author{
Simon Fontaine, Yi Yang, Bo Fan, Wei Qian and Yuwen Gu.

Maintainer: Simon Fontaine \email{fontaines@dms.umontreal.ca}
}

\examples{
#import package
library(MSTweedie)

#load data
data(AutoClaim)

# performs 10-folds CV with L1/Linf regularization
cv<-cv.MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')
}
\seealso{
\code{\link{MSTweedie}},
\code{\link{coef.cv.MSTweedie}},
\code{\link{plot.cv.MSTweedie}},
\code{\link{predict.cv.MSTweedie}}
}
