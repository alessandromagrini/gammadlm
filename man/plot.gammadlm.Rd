\name{plot.gammadlm}
\alias{plot.gammadlm}
\title{Graphics for the estimated lag distributions}
\description{Display the estimated lag distribution of each explanatory variable with lags.}
\usage{
\method{plot}{gammadlm}(x, x.names=NULL, conf=0.95, max.lag=NULL, max.quantile=0.999,
  xlim=NULL, ylim=NULL, add.legend=TRUE, cex.legend=1, digits=4, grid.length=100,
  main=NULL, ylab=NULL, xlab=NULL, ...)
}
\arguments{
  \item{x}{Object of class \code{gammadlm}.}
  \item{x.names}{Character vector including the name of the variables for which the lag distribution should be displayed. Unknown variables will be ignored. If \code{NULL} (the default), the lag distribution of all the variables with lags will be displayed.}
  \item{conf}{Numerical value indicating the level of confidence bands. Default is 0.95. If \code{NULL}, confidence bands will not be displayed.}
  \item{max.lag}{Non-negative integer value indicating the lag up to which each lag distribution should be displayed. If \code{NULL} (the default), it is set accordingly to argument \code{max.quantile}.}
  \item{max.quantile}{Numerical value indicating the order of the quantile up to which each lag distribution should be displayed. Default is 0.999 (99.9th percentile).} 
  \item{xlim}{Numerical vector of length 2 indicating the range of the x-axis, which is applied to all graphics (optional).}
  \item{ylim}{Numerical vector of length 2 indicating the range of the y-axis, which is applied to all graphics (optional).}
  \item{add.legend}{Logical value indicating whether a legend with numerical information should be added to the graphics. Default is \code{TRUE}.}
  \item{cex.legend}{Size of the legend. Default is 1.}
  \item{digits}{Integer non-negative value indicating the number of decimal places to be used in the legend. Default is 4. Ignored if \code{add.legend}=\code{FALSE}.}
  \item{grid.length}{Numerical value no less than 100 indicating the resolution of the interpolation. Default is 100.}
  \item{main}{Vector of characters including the title for each graphic. If \code{NULL} (the default), the name of the explanatory variables is used.}
  \item{ylab}{Text for y-axis, which is applied to all graphics (optional).}
  \item{xlab}{Text for x-axis, which is applied to all graphics (optional).}
  \item{...}{Further graphical parameters.}
}
%\value{}
%\details{}
%\note{}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\link{gammadlm}.}
\examples{
data(btc)
mydata <- btc[which(btc$Date>="2020-04-01"),]
mydataLR <- preProcess(time="Date", data=mydata, box.cox=0, ndiff=1)

dval <- list(DJA=0.85,IXIC=0.75,GSPC=0.55)
lval <- list(DJA=0.05, IXIC=0.35,GSPC=0.45)
m0 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR,
  control=list(delta.lim=dval, lambda.lim=lval))

plot(m0)  ## all the lag distributions
plot(m0, x.names=c("DJA","IXIC"))  ## just the ones of 'DJA' and 'IXIC'
}
