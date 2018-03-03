\name{MSTweedie}
\alias{MSTweedie}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Regularization path for the Multi-source sparse Tweedie model
}
\description{
This function fits the sparse Tweedie model on multi-source datasets along a sequence of regularization parameters lambda. The optimization is done by a Fortran95 routine.
}
\usage{
MSTweedie(x, y, w, source, rho = 1.5,
      nlambda = 100, lambda.min, lambda, x.normalize = T,
      eps, sr = T, kktstop = F, reg = c("L2", "Linf"),
      alpha = 0, dfmax = nvars + 1, pmax = min(dfmax * 1.2, nvars),
      pf = rep(1, nvars), maxit = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Either (1) a data frame containing the predictors, the responses (identifying the sources either by different columns in the simultaneous case or via an additionnal index column) and, optionnaly, the observation weigths or (2) a list of matrices containing only the predictors (mostly used internally for cross-validation.)
}
  \item{y}{
Either (1) a single integer identifying the column of \code{x} containing the response (requires \code{source} to be specified), (2) a vector of integers indentifying which columns of \code{x} are the responses (simultaneous case) or (3) a list of vector of responses (mostly used internally for cross-validation.)
}
  \item{w}{
(Optional) Either (1) a single integer identifying the column of \code{x} containing the observation weights or (2) a list of vector of weights (mostly used internally for cross-validation.) If this argument is missing, equal weight is assumed.
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
  \item{x.normalize}{
Logical flag for stadardization of the predictors prior to fitting the model. If \code{TRUE}, each predictors in each source is centered to zero and scaled to variance 1. After the fit of the model, the coefficients are returned on the original scale. Default is \code{FALSE}.
}
  \item{eps}{
Convergence threshold. Default is 1e-3.
}
  \item{sr}{
Logical flag for using the strong rule in the fit. Default is \code{TRUE}.
}
  \item{kktstop}{
Logical flag for using the KKT conditions to stop the fit before the end of the regularization parameter sequence. Default is \code{FALSE}.
}
  \item{reg}{
Either \code{"Linf"} for using \eqn{L_\infty}-regularization in the fit or \code{"L2"} for the \eqn{L_2}-regularization. Default is \code{"Linf"}.
}
  \item{alpha}{
Parameter controlling the balance between across-feature and within-feature sparsity in the penalty term
\deqn{(1-\alpha)||\beta||_q +\alpha||\beta||_1.} Possible range is [0,1], default is 0.
}
  \item{dfmax}{
Maximum number of variables included in the model at a single time. Default is \code{nvars+1}.
}
  \item{pmax}{
Limits the number of features ever to be nonzero. The difference with \code{dfmax}, is that if, a variable eventually exits the model, it will still be counted here. Default is \code{min(dfmax*1.2,nvars)}.
}
  \item{pf}{
Penalty weights in the penalty term by feature. Mostly used intternaly when the Adaptive Lasso is used in cross-validation. Expects a vector of length \code{nvars}, default is 1.
}
  \item{maxit}{
Maximum number of inner-loop iterations. Default is 10,000.
}
}
\details{
The sequence of regularization parameters implies a sequence of models fitted by the IRLS-BSUM algorithm described in the reference. For each value of the parameter, this function yield a model optimizing the penalzed Tweedie log-likelihood of multi-source data. The type of sparsity can be controlled by the arguments \code{reg} and \code{alpha}.

The computation time is influence by the arguments \code{eps}, \code{nlambda}, \code{lambda.min} (or \code{lambda}) and \code{maxit}. Consider ajusting these parameters to speed up computation. Small values of regularization parameters are the often the longest to fit; the \code{kktstop} argument can stop the algorithm before the end if convergence is judged sufficient in term of KKT conditions.

To pass sources with missing features compared to other sources, simply add a column of zero instead.
}
\value{
An object with S3 class \code{MSTweedie} :
\item{beta0}{A \code{ntaks*nlambda} matrix of parameter estimates for the intercept.}
\item{beta}{A list of length \code{nlambda} containing \code{nvars*ntaks} matrix of parameter estimates for the features.}
\item{df}{The number of included variables along the regularization path.}
\item{lambda}{The sequence of regularization parameters.}
\item{npasses}{The number of inner-loop iterations.}
\item{idvars}{The index of the variables in order of inclusion in the model.}
\item{dim}{The dimesions of the model (\code{nvars,ntasks}).}
\item{call}{The original call that produce this object.}
\item{pf}{The penalty factors for the features.}
\item{eps}{The convergence threshold used in the algorithm.}
\item{kkt}{A \code{nvars*ntasks*nlambda} array containing the values of the KKT conditions.}
\item{norm}{A \code{nvars*nlambda} matrix containing the norm of the features along the regularization path.}
\item{reg}{The type of regularization used in the algorithm.}
\item{alpha}{The value of the argument \code{alpha} used.}
\item{y}{A list of length \code{ntasks} containing the vectors of the responses for each source.}
\item{x}{A list of length \code{ntasks} containing matrices of the features for each source.}
\item{w}{A list of length \code{ntasks} containing the vectors of the observation weights for each source.}
\item{rho}{The power of the mean-variance relation used in the algorithm.}
\item{M}{A \code{nvars*ntasks*nlambda} array containing flags for the KKT conditions.}
\item{time}{Computing time.}
}
\references{
Fontaine, S., Yang, Y., Fan, B., Qian, W. and Gu, Y. (2018). "A Unified Approach to Sparse Tweedie Model
with Big Data Applications to Multi-Source
Insurance Claim Data Analysis," to be submitted.
}
\author{
Simon Fontaine, Yi Yang, Bo Fan, Wei Qian and Yuwen Gu.

Maintainer: Simon Fontaine \email{fontaines@dms.umontreal.ca}
}
\examples{
# import package
library(MSTweedie)

# load data
data(AutoClaim)

# fit the MSTweedie model with L1/Linf regularization
# y=1 sets CLM_AMT5 as the response
# source=4 sets REVOLKED as the source index
fit <- MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')

}
\seealso{
\code{\link{MSTweedie}},
\code{\link{coef.MSTweedie}},
\code{\link{print.MSTweedie}},
\code{\link{plot.MSTweedie}},
\code{\link{kkt.check}},
\code{\link{predict.MSTweedie}}
}
