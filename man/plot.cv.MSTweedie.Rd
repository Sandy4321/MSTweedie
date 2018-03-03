\name{plot.cv.MSTweedie}
\alias{plot.cv.MSTweedie}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Cross-Validation curve of cv.MSTweedie object
}
\description{
This function produces a plot of the mean Tweedie deviance against the logarithm of the solution path of a \code{cv.MSTweedie} object. The mean is accompanied by a one-standard-error interval. The regularization parameter value minimizing the CV deviance and the largest value within one standard error are indicated by vertical lines.
}
\usage{
\method{plot}{cv.MSTweedie}(x, sign.lambda = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{cv.MSTweedie} object.
}
  \item{sign.lambda}{
Whether to plo from left-to-right or the opposite. The Default is 1 corresponding to right-to-left.
}
  \item{\dots}{
Additionnal parameters passed to \code{\link{plot}}.
}
}
\value{
A plot is produced
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

# performs 10-folds CV with L1/Linf regularization
cv <- cv.MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')

# plot CV deviance mean and std. err., and lambda.min, lambda.1se
plot.cv.MSTweedie(cv)
}
\seealso{
\code{\link{cv.MSTweedie}},
\code{\link{plot}}
}
