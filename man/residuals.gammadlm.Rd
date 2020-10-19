\name{residuals.gammadlm}
\alias{residuals.gammadlm}
\title{Residuals and graphical model diagnostics}
\description{Obtain residuals and graphical diagnostics for an estimated Gamma distributed-lag model.}
\usage{
\method{residuals}{gammadlm}(object, plot=FALSE, cex.lab=1, cex.axis=1, ...)
}
\arguments{
  \item{object}{Object of class \code{gammadlm}.}
  \item{plot}{Logical value indicating whether the graphical diagnostics should be produced.}
  \item{cex.lab}{Size of the axis text labels.}
  \item{cex.axis}{Size of the axis tick labels.}
  \item{...}{Further graphical parameters.}
}
\value{Residuals extracted from \code{object}, or \code{NULL} if argument \code{plot} is set to \code{TRUE}.}
%\details{}
%\note{}
\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\link{gammadlm}.}
\examples{
data(USstock)
mydata <- USstock[which(USstock$Date>="2020-04-01"),]
mydataLR <- tsDiff(time.name="Date", data=mydata, ndiff=1, log=TRUE)

dval <- list(DJA=0.85,IXIC=0.75,GSPC=0.55)
lval <- list(DJA=0.5, IXIC=0.35,GSPC=0.45)
mod <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR,
  control=list(delta.lim=dval, lambda.lim=lval))

residuals(mod)  ## returns the residuals
residuals(mod, plot=TRUE)  ## displays graphical model diagnostics
}
