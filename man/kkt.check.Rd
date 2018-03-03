\name{kkt.check}
\alias{kkt.check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
KKT conditions plots for MSTweedie objects
}
\description{
This functions plots the values of the KKT conditions in the solution path of a \code{MSTweedie} object.
}
\usage{
kkt.check(fit, eps = fit$eps, cond = c(1, 2),
      from = 1, to = length(fit$lambda), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
\code{MSTweedie} object.
}
  \item{eps}{
Value of the threshold to add to the plots. Default is the convergence threshold of the \code{fit} object, i.e. \code{fit$eps}.
}
  \item{cond}{
A vector of integers between 1 and 4 indicating which conditions to plot. See "Details".
}
  \item{from}{
The index of the first regularization parameter to include. Default is 1.
}
  \item{to}{
The index of the last regularization parameter to include. Default is \code{length(fit$lambda)}.
}
  \item{\dots}{
Additionnal graphical parameters to pass to \code{plot}.
}
}
\details{
The reference contains detailed explanation of these KKT conditions.

If \code{code} contains 1, then the sequence of values for the non-zero KKT condition is plotted:
\deqn{||U_j||=\lambda v_j,\beta_j\neq 0.}

If \code{code} contains 2, then the sequence of values for the zero KKT condition is plotted:

\deqn{||U_j||\leq\lambda v_j,\beta_j= 0.}

If \code{code} contains 3, then the sequence of the detailed values for the non-zero KKT condition is plotted:

\deqn{-U_j/\lambda v_j=0,\beta_j\neq 0, k\not\in M(\beta_j).}

If \code{code} contains 4, then the sequence of the aggregate KKT condition is plotted:
\deqn{\frac{-\sum_{j=1}^pU_j^\top\beta_j}{\sum_{j=1}^pv_j||\beta_j||} = \lambda.}
}
\value{
Multiple plots are produced.
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

# plot the first two kkt conditions
par(mfrow=c(2,1))
kkt.check(fit)
}
\seealso{
\code{\link{MSTweedie}}
}
