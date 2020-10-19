\name{lagCoef}
\alias{lagCoef}
\title{Estimated dynamic coefficients}
\description{Estimated dynamic coefficients for each explanatory variable.}
\usage{lagCoef(x, cumulative=FALSE, max.lag=NULL, max.quantile=0.999)}
\arguments{
  \item{x}{An object of class \code{gammadlm}.}
  \item{cumulative}{Logical value indicating whether cumulative coefficients should be returned. Default is \code{FALSE}.}
  \item{max.lag}{Non-negative integer value indicating the lag up to which coefficients should be returned. If \code{NULL} (the default), it is set accordingly to argument \code{max.quantile}.}
  \item{max.quantile}{Numerical value indicating the order of the quantile lag up to which coefficients should be returned. Default is 0.999 (99.9th percentile).}  
}
\value{A list with one component for each explanatory variable with lags.
Each component is an object of class \code{data.frame} with lags as observations and two columns containing estimation and asymptotic standard error.}
%\details{}
\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\link{gammadlm}; \link{plot.gammadlm}.}
\examples{
data(USstock)
mydata <- USstock[which(USstock$Date>="2020-04-01"),]
mydataLR <- tsDiff(time.name="Date", data=mydata, ndiff=1, log=TRUE)

dval <- list(DJA=0.85,IXIC=0.75,GSPC=0.55)
lval <- list(DJA=0.5, IXIC=0.35,GSPC=0.45)
mod <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR,
  control=list(delta.lim=dval, lambda.lim=lval))

lagCoef(mod)  ## coefficients shown up to the 99.9th percentile lag
lagCoef(mod, max.lag=11)  ## coefficients shown up to lag 11
lagCoef(mod, cumulative=TRUE)  ## cumulative coefficients
}