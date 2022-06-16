\name{lagCoef}
\alias{lagCoef}
\title{Estimated dynamic coefficients}
\description{See the estimated dynamic coefficients for each explanatory variable with lags.}
\usage{lagCoef(x, conf=0.95, cumulative=FALSE, max.lag=NULL, max.quantile=0.999)}
\arguments{
  \item{x}{An object of class \code{gammadlm}.}
  \item{conf}{Numerical value indicating the level of confidence intervals. If \code{NULL}, confidence intervals are not provided. Default is 0.95.}
  \item{cumulative}{Logical value indicating whether cumulative coefficients should be returned. Default is \code{FALSE}.}
  \item{max.lag}{Non-negative integer value indicating the lag up to which coefficients should be returned. If \code{NULL} (the default), it is set accordingly to argument \code{max.quantile}.}
  \item{max.quantile}{Numerical value indicating the order of the quantile lag up to which coefficients should be returned. Default is 0.999 (99.9th percentile).}  
}
\value{A list with one component for each explanatory variable with lags.
Each component is an object of class \code{data.frame} containing estimation, asymptotic standard error and confidence interval for each lag up to \code{max.lag}.}
%\details{}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\link{gammadlm}; \link{plot.gammadlm}.}
\examples{
data(btc)
mydata <- btc[which(btc$Date>="2020-04-01"),]
mydataLR <- preProcess(c("BTC","DJA","IXIC","GSPC"),
  time="Date", data=mydata, box.cox=0, ndiff=1)

dval <- list(DJA=0.85,IXIC=0.75,GSPC=0.55)
lval <- list(DJA=0.05, IXIC=0.35,GSPC=0.45)
m0 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR,
  control=list(delta.lim=dval, lambda.lim=lval))

lagCoef(m0)  ## coefficients shown up to the 99.9th percentile lag
lagCoef(m0, max.lag=11)  ## coefficients shown up to lag 11
lagCoef(m0, cumulative=TRUE)  ## cumulative coefficients
}