\name{predict.MSTweedie}
\alias{predict.MSTweedie}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction based on a MSTweedie object
}
\description{
This function predicts the response (or the link function) given features at some regularization parameter.
}
\usage{\method{predict}{MSTweedie}(fit, newx, tasks, s, type = c("response", "link"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
\code{MSTweedie} object.
}
  \item{newx}{
A list of matrices all with \code{nvars} columns used to predict the responses. Each matrix is associated with a source from the original \code{fit} MSTweedie object in the order of argument \code{tasks}. Default is the original dataset (\code{fit$x}).
}
  \item{tasks}{
A vector of the same length as \code{newx} containing the index of the source associated with the respective matrix in \code{newx}. Default is all sources in original order.
}
  \item{s}{
Either an integer identifying the index of regularization parameter to use or the parameter itself (must match one of \code{fit}.)
}
  \item{type}{
Whether to return the \code{"response""} or the \code{"link""} (linear predictor).
}
}
\value{
A list of length \code{length(tasks)} of vectors corresponding to the prediction of either the reponses or the link function (see \code{type}).
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
fit <- MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')

# predict first source at 34th lambda
head(predict.MSTweedie(fit, s=34L)[[1]])
}
\seealso{
\code{\link{MSTweedie}},
\code{\link{coef.MSTweedie}}
}
