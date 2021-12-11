\name{adfTest}
\alias{adfTest}
\title{Augmented Dickey-Fuller test}
\description{Augmented Dickey-Fuller test with automated selection of the lag length for a multivariate time series, eventually structured as a panel.}
\usage{adfTest(var.names=NULL, unit=NULL, time=NULL, data, box.cox=1, ndiff=0, max.lag=NULL)}
\arguments{
  \item{var.names}{Character vector including the name of the variables to be differenced. If \code{NULL} (the default), all the quantitative variables in the dataset provided to argument \code{data} will be differenced, with the exception of the variable indicated in \code{time}.}
  \item{unit}{Character containing the name of the units of observation in case of a panel structure. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), a single unit of observation is assumed.}
  \item{time}{Character containing the name of the time variable, which must be in numeric or date format. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), data are assumed to be temporally ordered.}
  \item{data}{Object of class \code{data.frame} containing the variables in \code{var.names}, \code{unit} and \code{time}.}
  \item{box.cox}{Named vector indicating the parameter of the Box-Cox transformation for variables in \code{x.names}. If \code{box.cox} is of length one and has no names, the same parameter is used for all variables in \code{x.names}. If a variable in \code{x.names} has no name in \code{box.cox}, value 1 is assumed, meaning no transformation. Value 0 of the parameter equates to the logarithmic transformation. Default value of \code{box.cox} is 1, meaning no transformation for all variables in \code{x.names}.}
  \item{ndiff}{Named vector with non-negative integer values indicating the number of differences for variables in \code{x.names}. If \code{ndiff} is of length one and has no names, the same number of differences is used for all variables in \code{x.names}. If a variable in \code{x.names} has no name in \code{ndiff}, value 0 is assumed, meaning no differencing. Default value of \code{ndiff} is 0, meaning no differencing for all variables in \code{x.names}.}
  \item{max.lag}{Non-negative integer representing the maximum lag length at which the test should be performed. If \code{NULL} (the default), it is taken as the squared root of the length of the time series.
%\code{trunc((length(x)-1)^(1/3))}
.}
}
\value{One list for each variable in \code{var.names}, each with three components:
\itemize{
\item{\code{statistic}: }{the test statistic};
\item{\code{lag.selected}: }{the lag length selected to perform the test};
\item{\code{p.value}: }{the p-value of the test, which is a single value if \code{unit} is \code{NULL} and a vector of values otherwise, one for each unit of observation plus another one indicating the combined p-value};
\item{\code{box.cox}: }{parameter of the Box-Cox transformation};
\item{\code{ndiff}: }{order of differencing}.
}
}
\note{The first order difference of logarithmic values (\code{box.cox}=0 and \code{ndiff}=1) provides the log returns, which approximate the proportional changes with respect to the previous time point.

If the parameter of the Box-Cox transformation is set to 0 for a non-negative variable, the transformation will not be applied. i.e., the parameter is set to 1, and a warning is returned.}
\details{The variable subjected to the test must be quantitative.

The null hypothesis is the presence of unit roots, thus p-values higher than 0.05 indicate the absence of weak stationarity and thus the need of differencing the time series.

The lag length at which the test is performed is selected through backward elimination starting from the lag length specified in argument \code{max.lag}.

In case of a panel structure, p-values are combined according to the method by Demetrescu \emph{et al.} (2006).

Missing values, if present, are deleted out and a warning is returned.}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
M. Demetrescu, U. Hassler, and A. Tarcolea (2006). Combining significance of correlated statistics with application to panel data. \emph{Oxford Bulletin of Economics and Statistics}, 68(5), 647-663.
  
D. A. Dickey, and W. A. Fuller (1981). Likelihood ratio statistics for autoregressive time series with a unit root. \emph{Econometrica}, 49(4): 1057-1072.

%S. Ng, and W. P. Perron (2001). Unit Root Tests in ARMA Models with Data-Dependent Methods for the Selection of the Truncation Lag. \emph{Journal of the American Statistical Association}, 90(429): 268-281.
}
\seealso{\link{tsDiff}.}
\examples{
data(USstock)
mydata <- USstock[which(USstock$Date>="2020-04-01"),]

# adf test for time series in level ('box.cox'=1 by default):
#   there are some unit roots
adfTest(var.names=c("BTC","DJA","IXIC","GSPC"), time="Date",
  data=mydata)
  
# adf test for time series in log return: ok
adfTest(var.names=c("BTC","DJA","IXIC","GSPC"), time="Date",
  data=mydata, box.cox=0, ndiff=1)
}