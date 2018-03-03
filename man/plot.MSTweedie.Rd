\name{plot.MSTweedie}
\alias{plot.MSTweedie}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the solution path of a MSTweedie object
}
\description{
This function produces coefficient profile plot(s) along the sequence of regularization parameters of a \code{MSTweedie} object.
}
\usage{
\method{plot}{MSTweedie}(x, log.lambda = TRUE, type.coef = c("coef", "norm"),
      lambda.min = NA, lambda.1se = NA, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{MSTweedie} object.
}
  \item{log.lambda}{
Logical flag: whether to plot against the logarithm of the regularization parameters (default) or against the reg. parameters.
}
  \item{type.coef}{
Whether to plot the coefficients in each source (\code{"coef"}) or to plot the norm across all sources (\code{"norm"}), where the norm is given by \code{x$reg}.
}
  \item{lambda.min}{
If supplied, this traces a vertical line at \code{lambda.min} (given by CV.)
}
  \item{lambda.1se}{
If supplied, this traces a vertical line at \code{lambda.1se} (given by CV.)
}
  \item{\dots}{
Additionnal parameters passed to \code{\link{plot}}.
}
}
\value{
Plot(s) are produced.
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

# plot solution path of the norm of the coefficients
plot.MSTweedie(fit, type.coef='norm')

# plot solution path of the the coefficients
par(mfrow=c(2,1))
plot.MSTweedie(fit, type.coef='coef')
}
\seealso{
\code{\link{MSTweedie}},
\code{\link{plot}},
\code{\link{coef.MSTweedie}}
}
