\name{predict.cv.MSTweedie}
\alias{predict.cv.MSTweedie}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction based on a cv.MSTweedie object
}
\description{
This function is a wrapper function for the \code{predict.MSTweedie} when applied to a \code{cv.MSTweedie} object.
}
\usage{
\method{predict}{cv.MSTweedie}(object, newx, s = c("lambda.1se", "lambda.min"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
\code{cv.MSTweedie} object.
}
  \item{newx}{
A list of matrices all with \code{nvars} columns used to predict the responses. Each matrix is associated with a source from the original \code{object} MSTweedie object in the order of argument \code{tasks}. Default is the original dataset (\code{object$MSTweedie.fit$x}).
}
  \item{s}{
Either a vector of regularization parameters (must match those of \code{object}), a vector of indices of regularization parameters. Default is the whole solution path or one of \code{"lambda.min"} or \code{"lambda.1se"}. Default is \code{"lambda.1se"}.
}
}
\value{
A list of length \code{ntasks} of vectors corresponding to the prediction of either the reponses or the link function (see \code{type}).
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

# extract coefficients at lambda.1se
head(predict.cv.MSTweedie(cv)[[1]])
}
\seealso{
\code{\link{cv.MSTweedie}},
\code{\link{predict.MSTweedie}}
}
