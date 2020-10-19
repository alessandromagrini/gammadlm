\name{adfTest}
\alias{adfTest}
\title{Augmented Dickey-Fuller test}
\description{Augmented Dickey-Fuller test with automated selection of the lag length.}
\usage{adfTest(x, max.lag=NULL)}
\arguments{
  \item{x}{Numerical vector containing temporally ordered data.}
  \item{max.lag}{Non-negative integer representing the maximum lag length at which the test should be performed. If \code{NULL} (the default), it is taken as \code{trunc((length(x)-1)^(1/3))}}
}
\value{A list with three components:
\itemize{
\item{\code{statistic}: }{the test statistic};
\item{\code{lag.selected}: }{the lag length selected to perform the test};
\item{\code{p.value}: }{the p-value of the test}.
}
}
\details{The variable subjected to the test must be quantitative.

The test is performed starting from a specified lag length, which is sequentially decreased until the test statistic is lower than 1.6 in absolute value or the lag length is 0.

The null hypothesis is the presence of unit root, thus p-values higher than 0.05 indicate non-stationarity and the need of differencing the time series.

Missing values, if present, are deleted out and a warning is returned.}
\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
D. A. Dickey, W. A. Fuller (1981). Likelihood ratio statistics for autoregressive time series with a unit root. \emph{Econometrica}, 49(4), 1057-1072.
}
\seealso{\link{tsDiff}.}
\examples{
data(USstock)
mydata <- USstock[which(USstock$Date>="2020-04-01"),]
adfTest(mydata$BTC)  ## unit root
adfTest(diff(log(mydata$BTC)))  ## data in log return: ok
}