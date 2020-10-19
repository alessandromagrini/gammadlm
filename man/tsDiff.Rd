\name{tsDiff}
\alias{tsDiff}
\title{Time series differencing}
\description{Application of differencing to several time series.}
\usage{tsDiff(var.names=NULL, time.name=NULL, data, ndiff=0, log=F)}
\arguments{
  \item{var.names}{Character vector including the name of the variables to be differenced. If \code{NULL} (the default), all the quantitative variables in the dataset provided to argument \code{data} will be differenced, with the exception of the variable indicated in \code{time.name}.}
  \item{time.name}{Character including the name of the time variable, which must be in numeric or date format. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), data are assumed to be temporally ordered.}
  \item{data}{Object of class \code{data.frame} containing the variables in \code{var.names} and \code{time.name}.}
  \item{ndiff}{Non-negative integer value indicating the order of differencing for all variables in \code{var.names}, or a vector of non-negative integer values indicating the order of differencing for each variable in \code{var.names}. In the second case, the length of \code{ndiff} must be the same of \code{var.names}. Default is 0, meaning no differencing.}
  \item{log}{Logical value indicating whether the logarithmic transformation should be applied to all variables in \code{var.names}, or a vector of logical values indicating whether the logarithmic transformation should be applied to each variable in \code{var.names}. In the second case, the length of \code{log} must be the same of \code{var.names}. Default is \code{FALSE}, meaning no logarithmic transformation.}
}
\value{The object provided to argument \code{data} where variables in \code{var.names} have been differenced and/or log transformed.}
%\details{}
\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\link{adfTest}.}
\examples{
data(USstock)
mydata <- USstock[which(USstock$Date>="2020-04-01"),]

# setting ndiff=0 and log=T produces the log return
mydataLR <- tsDiff(var.names=c("BTC","DJA","IXIC","GSPC"), time.name="Date",
  data=mydata, ndiff=1, log=TRUE)
summary(mydataLR)

# same result by omitting 'var.names'
mydataLR2 <- tsDiff(time.name="Date", data=mydata, ndiff=1, log=TRUE)
summary(mydataLR2)
}