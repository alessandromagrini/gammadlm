\name{tsDiff}
\alias{tsDiff}
\title{Time series differencing}
\description{Application of differencing to a multivariate time series, eventually structured as a panel.}
\usage{tsDiff(var.names=NULL, panel=NULL, time=NULL, data, log=FALSE, ndiff=0)}
\arguments{
  \item{var.names}{Character vector including the name of the variables to be differenced. If \code{NULL} (the default), all the quantitative variables in the dataset provided to argument \code{data} will be differenced, with the exception of the variable indicated in \code{time}.}
  \item{panel}{Character containing the name of the units of observation in case of a panel structure. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), a single unit of observation is assumed.}
  \item{time}{Character containing the name of the time variable, which must be in numeric or date format. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), data are assumed to be temporally ordered.}
  \item{data}{Object of class \code{data.frame} containing the variables in \code{var.names}, \code{panel} and \code{time}.}
  \item{log}{Logical value indicating whether the logarithmic transformation should be applied to all variables in \code{var.names}, or a vector of logical values indicating whether the logarithmic transformation should be applied to each variable in \code{var.names}. In the second case, the length of \code{log} must be the same of \code{var.names}. Default is \code{FALSE}, meaning no logarithmic transformation.}
  \item{ndiff}{Non-negative integer value indicating the order of differencing for all variables in \code{var.names}, or a vector of non-negative integer values indicating the order of differencing for each variable in \code{var.names}. In the latter case, the length of \code{ndiff} must be the same of \code{var.names}. Default is 0, meaning no differencing.}
}
\value{The object provided to argument \code{data} where variables in  \code{var.names} have been log transformed and/or differenced.}
\note{The first order difference of logarithmic values provides the log returns, which approximate the proportional changes with respect to the previous time point.}
%\details{}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\link{adfTest}.}
\examples{
data(USstock)
mydata <- USstock[which(USstock$Date>="2020-04-01"),]

# setting ndiff=0 and log=T produces the log returns
mydataLR <- tsDiff(var.names=c("BTC","DJA","IXIC","GSPC"), time="Date",
  data=mydata, log=TRUE, ndiff=1)
summary(mydataLR)

# same result by omitting 'var.names'
mydataLR2 <- tsDiff(time="Date", data=mydata, log=TRUE, ndiff=1)
summary(mydataLR2)
}