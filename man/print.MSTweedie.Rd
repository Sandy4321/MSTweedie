\name{print.MSTweedie}
\alias{print.MSTweedie}
\title{
Print the results of a MSTweedie object
}
\description{
This function prints the number of variables included in the model and the value of the regularization parameter along the solution path.
}
\usage{\method{print}{MSTweedie}(x, digits = max(3, getOption("digits") - 3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{MSTweedie} object.
}
  \item{digits}{
Number of significant digits in print display.
}
}
\value{
A two column matrix with \code{nlambda} rows containing the number of variables included in each model with the value of the regularization parameter.
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
fit <- MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')

# prints number of selected variables along solution path
print.MSTweedie(fit)
}
\seealso{
\code{\link{MSTweedie}}
}
