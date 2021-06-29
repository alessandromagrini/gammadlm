\name{adfTest}
\alias{adfTest}
\title{Augmented Dickey-Fuller test}
\description{Augmented Dickey-Fuller test with automated selection of the lag length for a multivariate time series, eventually structured as a panel.}
\usage{adfTest(var.names=NULL, panel=NULL, time=NULL, data, log=FALSE, ndiff=0, max.lag=NULL)}
\arguments{
  \item{var.names}{Character vector including the name of the variables to be differenced. If \code{NULL} (the default), all the quantitative variables in the dataset provided to argument \code{data} will be differenced, with the exception of the variable indicated in \code{time}.}
  \item{panel}{Character containing the name of the units of observation in case of a panel structure. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), a single unit of observation is assumed.}
  \item{time}{Character containing the name of the time variable, which must be in numeric or date format. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), data are assumed to be temporally ordered.}
  \item{data}{Object of class \code{data.frame} containing the variables in \code{var.names}, \code{panel} and \code{time}.}
  \item{log}{Logical value indicating whether the logarithmic transformation should be applied to all variables in \code{var.names}, or a vector of logical values indicating whether the logarithmic transformation should be applied to each variable in \code{var.names}. In the second case, the length of \code{log} must be the same of \code{var.names}. Default is \code{FALSE}, meaning no logarithmic transformation.}
  \item{ndiff}{Non-negative integer value indicating the order of differencing for all variables in \code{var.names}, or a vector of non-negative integer values indicating the order of differencing for each variable in \code{var.names}. In the latter case, the length of \code{ndiff} must be the same of \code{var.names}. Default is 0, meaning no differencing.}
  \item{max.lag}{Non-negative integer representing the maximum lag length at which the test should be performed. If \code{NULL} (the default), it is taken as \code{trunc((length(x)-1)^(1/3))}.}
}
\value{One list for each variable in \code{var.names}, each with three components:
\itemize{
\item{\code{statistic}: }{the test statistic};
\item{\code{lag.selected}: }{the lag length selected to perform the test};
\item{\code{p.value}: }{the p-value of the test, which is a single value if \code{panel} is \code{NULL} and a vector of values otherwise, one for each unit of observation plus another one indicating the combined p-value};
\item{\code{log}: }{logical value indicating whether the logaritmic transformation has been applied};
\item{\code{ndiff}: }{the order of differencing}.
}
}
\note{The first order difference of logarithmic values provides the log returns, which approximate the proportional changes with respect to the previous time point.}
\details{The variable subjected to the test must be quantitative.

The null hypothesis is the presence of unit root, thus p-values higher than 0.05 indicate the absence of weak stationarity and thus the need of differencing the time series.

The test is performed starting from a specified lag length, which is sequentially decreased until the test statistic is lower than 1.6 in absolute value or the lag length is 0 (Ng \& Perron, 2001).

In case of a panel structure, p-values are combined according to the method by Demetrescu \emph{et al.} (2006).

Missing values, if present, are deleted out and a warning is returned.}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
M. Demetrescu, U. Hassler, and A. Tarcolea (2006). Combining Significance of Correlated Statistics with Application to Panel Data. \emph{Oxford Bulletin of Economics and Statistics}, 68(5), 647-663.
  
D. A. Dickey, and W. A. Fuller (1981). Likelihood Ratio Statistics for Autoregressive Time Series with a Unit Root.

S. Ng, and W. P. Perron (2001). Lag Length  Selection and the Construction of Unit Root Tests with Good Size and Power. \emph{Econometrica}, 60: 1519-1554..
}
\seealso{\link{tsDiff}.}
\examples{
data(USstock)
mydata <- USstock[which(USstock$Date>="2020-04-01"),]

# adf test for time series in level: there are some unit roots
adfTest(var.names=c("BTC","DJA","IXIC","GSPC"), time="Date",
  data=mydata, log=TRUE)
  
# adf test for time series in log return: ok
adfTest(var.names=c("BTC","DJA","IXIC","GSPC"), time="Date",
  data=mydata, log=TRUE, ndiff=1)
}