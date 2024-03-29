\name{lagCoef}
\alias{lagCoef}
\title{Estimated dynamic multipliers}
\description{See the estimated dynamic multipliers for each explanatory variable with lags.}
\usage{lagCoef(x, conf=0.95, cumulative=FALSE, max.lag=NULL, max.quantile=0.999)}
\arguments{
  \item{x}{An object of class \code{gammadlm}.}
  \item{conf}{Numerical value indicating the level of confidence intervals. If \code{NULL}, confidence intervals are not provided. Default is 0.95.}
  \item{cumulative}{Logical value indicating whether cumulative multipliers should be returned. Default is \code{FALSE}.}
  \item{max.lag}{Non-negative integer value indicating the lag up to which multipliers should be returned. If \code{NULL} (the default), it is set accordingly to argument \code{max.quantile}.}
  \item{max.quantile}{Numerical value indicating the order of the quantile lag up to which multipliers should be returned. Default is 0.999 (99.9th percentile).}  
}
\value{A list with one component for each explanatory variable with lags.
Each component is an object of class \code{data.frame} containing estimation, asymptotic standard error and confidence interval for each lag up to \code{max.lag}.}
%\details{}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\link{gammadlm}; \link{plot.gammadlm}.}
\examples{
\dontrun{
data(btc)
mydata <- btc[which(btc$Date>="2020-04-01"),]
mydataLR <- preProcess(c("BTC","DJA","IXIC","GSPC"),
  time="Date", data=mydata, box.cox=0, ndiff=1)
m0 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR)

lagCoef(m0)  ## multipliers shown up to the 99.9th percentile lag
lagCoef(m0, max.lag=11)  ## multipliers shown up to lag 11
lagCoef(m0, cumulative=TRUE)  ## cumulative multipliers}
}