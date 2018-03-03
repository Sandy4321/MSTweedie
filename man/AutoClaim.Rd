\name{AutoClaim}
\alias{AutoClaim}
\docType{data}
\title{
AutoClaim dataset
}
\description{
The motor insurance dataset is originially retrieved the \code{\link{cplm}} package. It contains insurance claim data as well as information on the policyholder. From the original dataset, only part of the variables are kept and some are transformed (see description below). Missing values are imputated via the \code{\link{rrcovNA::impSeq}} function.
}
\usage{data("AutoClaim")}
\format{
  A data frame with 10296 observations on the following 35 variables.
  \describe{
    \item{\code{CLM_AMT5}}{Aggregate claim loss of policy (in thousands)}
    \item{\code{KIDSDRIV}}{Number of child passengers}
    \item{\code{TRAVTIME}}{Commute time}
    \item{\code{CAR_USE}}{(1) Private or (2) commercial use}
    \item{\code{BLUEBOOK}}{(log) car value}
    \item{\code{NPOLICY}}{Number of policies}
    \item{\code{RED_CAR}}{Whether the color of the car is (2) car or (1) not}
    \item{\code{REVOLKED}}{Whether the policyholder's license was (2) revoked in the past or (1) not}
    \item{\code{MVR_PTS}}{Number of motor vehicule record points}
    \item{\code{HOMEKIDS}}{Number of children at home}
    \item{\code{GENDER}}{Gender of policyholder : (1) female or (2) male}
    \item{\code{MARRIED}}{Whether the policyholder is (2) married or (1) not}
    \item{\code{PARENT1}}{Whether (2) the policyholder grew up in a single-parent family or (1) not}
    \item{\code{AREA}}{(1) Rural or (2) urban area}
    \item{\code{CAR_TYPE_2-6}}{(0-1 dummy variables) Type of car :  (base) Panel Truck, (2) Pickup,(3) Sedan, (4) Sports Car, (5) SUV, (6) Van}
    \item{\code{JOBCLASS_2-9}}{(0-1 dummy variables) Job class of policyholder: (base) Unknown, (2) Blue Collar, (3) Clerical, (4) Doctor, (5) Home Maker, (6) Lawyer, (7) Manager, (8) Professionnal, (9) Student}
    \item{\code{MAX_EDUC_2-5}}{(0-1 dummy variables) Maximal level of education of policyholder: (base) less than High School, (2) Bachelors, (3) High School, (4) Masters, (5) PhD}
    \item{\code{AGE_CAT_2-5}}{(0-1 dummy variables) Age category of policyholder : (base) <30, (2) [30,40), (3) [40,50), (4) [50,60), (5) 60+].}
  }
}
\source{

\code{cplm} package.
}
\references{
Fontaine, S., Yang, Y., Fan, B., Qian, W. and Gu, Y. (2018). "A Unified Approach to Sparse Tweedie Model
with Big Data Applications to Multi-Source
Insurance Claim Data Analysis," to be submitted.

Zhang, Y. (2013). "\code{cplm}: Compound Poisson Linear Models." A vignette for \code{R} package \code{cplm}. Available from \url{http://cran.r-project.org/web/packages/cplm}.

Todorov, V. (2016). "\code{rrcovNA}: Scalable Robust Estimators with High Breakdown Point for Incomplete Data." A vignette for \code{R} package \code{rrcovNA}. Available from \url{https://cran.r-project.org/web/packages/rrcovNA}.
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

#display head of dataset
head(AutoClaim)

#classify the policies by REVOLKED and whether there was a claim or not
table(AutoClaim$REVOLKED, AutoClaim$CLM_AMT5 > 0)
}
\keyword{datasets}
