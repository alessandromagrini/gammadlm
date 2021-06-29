\name{gammaWeights}
\alias{gammaWeights}
\alias{gammaQuantile}
\alias{gammaKernel}
\title{Functionalities for the Gamma lag distribution}
\description{Obtain weights, quantiles and kernel projection for the desired Gamma lag distribution.}
\usage{gammaWeights(k, par, offset=0, normalize=TRUE)
gammaQuantile(prob, par, offset=0)
gammaKernel(x, par, panel=NULL, offset=0, normalize=TRUE)
}
\arguments{
  \item{k}{Numerical vector indicating the lags for which the weights should be computed.}
  \item{prob}{Numerical vector indicating the order of the quantiles to be computed.}
  \item{x}{Numerical vector representing temporally ordered data for which the kernel projection should be returned.}
  \item{par}{Numerical vector of length 2 representing the shape parameters of the Gamma lag distribution.}
  \item{panel}{Character containing the name of the units of observation in case of a panel structure. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), a single unit of observation is assumed.}
  \item{offset}{Numerical value representing the offset of the Gamma lag distribution. Default is 0.}
  \item{normalize}{Logical value indicating whether the weights should be normalized to have sum 1. Default is \code{TRUE}.}
}
%\value{}
\details{Function \code{gammaWeights} provides the weights, function \code{gammaQuantile} computes the quantiles, and function \code{gammaKernel} returns the kernel projection.}
%\note{}
%\details{}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
%\seealso{\link{gammadlm}.}
\examples{
## examples for a Gamma lag distribution with delta=0.6 and lambda=0.3

# weights
gammaWeights(0:12, par=c(0.6,0.3))  ## at lags from 0 to 12
gammaWeights(10, par=c(0.6,0.3))    ## at lag 10

# quantiles
gammaQuantile(0.5, par=c(0.6,0.3))   ## median
gammaQuantile(0.95, par=c(0.6,0.3))  ## 95th percentile
gammaQuantile(0.99, par=c(0.6,0.3))  ## 99th percentile

# kernel projection
set.seed(100); xval <- rnorm(10)
gammaKernel(xval, par=c(0.6,0.3))

# kernel projection under a panel structure
set.seed(100); xval <- rnorm(20)
gr <- c(rep(0,10),rep(1,10))
gammaKernel(xval, par=c(0.6,0.3), panel=gr)
}
