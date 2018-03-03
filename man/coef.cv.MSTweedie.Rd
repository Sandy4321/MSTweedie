\name{coef.cv.MSTweedie}
\alias{coef.cv.MSTweedie}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the estimated coefficients of a cv.MSTweedie object
}
\description{
This function is a wrapper function for the \code{coef.MSTweedie} when applied to a \code{cv.MSTweedie} object.
}
\usage{\method{coef}{cv.MSTweedie}(cv, s = c("lambda.1se", "lambda.min"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cv}{
\code{cv.MSTweedie} object
}
  \item{s}{
Either a vector of regularization parameters (must match those of \code{fit}), a vector of indices of regularization parameters. Default is the whole solution path or one of \code{"lambda.min"} or \code{"lambda.1se"}. Default is \code{"lambda.1se"}.
}
}
\details{
Returns the estimated coefficient in the \code{MSTweedie} object within the \code{cv.MSTweedie} object at the specified \code{s} values of the regularization parameter. When \code{"lambda.min"} or \code{"lambda.1se"} is supplied for \code{s}, the respective value of the regularization parameter is used.
}
\value{
A list of length \code{ntasks} of matrices of dimension \code{nvars*length(s)} containing the estimated coefficients at each values of \code{s}.
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
#import package
library(MSTweedie)

#load data
data(AutoClaim)

# performs 10-folds CV with L1/Linf regularization
cv <- cv.MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')

# extract coefficients at lambda.1se
coef.cv.MSTweedie(cv)
}
\seealso{
\code{\link{cv.MSTweedie}},
\code{\link{coef.MSTweedie}}
}
