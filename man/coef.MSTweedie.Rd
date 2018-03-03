\name{coef.MSTweedie}
\alias{coef.MSTweedie}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the extimated coefficients of a MSTweedie object
}
\description{
This function extracts the extimated coefficients of a MSTweedie object at given regularization parameters and given sources.
}
\usage{\method{coef}{MSTweedie}(fit, s = NULL, tasks)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
\code{MSTweedie} object.
}
  \item{s}{
Either a vector of regularization parameters (must match those of \code{fit}) or a vector of indices of regularization parameters. Default is the whole solution path.
}
\item{tasks}{
Indices of the sources of which coefficients must be returned. Default is all sources.
}

}

\value{
A list of length \code{length(tasks)} of matrices of dimension \code{nvars*length(s)} containing the estimated coefficients at each values of \code{s}.
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

# L1/Linf regularization
fit <- MSTweedie(x=AutoClaim, y=1, source=4, reg='Linf')

# extract coefficients at 34th to 36th lambda
coef.MSTweedie(fit, s=34:36)
}
\seealso{
\code{\link{MSTweedie}},
\code{\link{plot.MSTweedie}},
\code{\link{predict.MSTweedie}}
}
