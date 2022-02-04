\name{unirootTest}
\alias{unirootTest}
\title{Unit root tests}
\description{Augmented Dickey-Fuller (ADF) and Kwiatkowski-Phillips-Schmidt-Shin (KPSS) tests for a multivariate time series, eventually structured as a panel.}
\usage{unirootTest(var.names=NULL, unit=NULL, time=NULL, data, box.cox=1, ndiff=0, max.lag=NULL)}
\arguments{
  \item{var.names}{Character vector including the name of the variables to be differenced. If \code{NULL} (the default), all the quantitative variables in the dataset provided to argument \code{data} will be differenced, with the exception of the variable indicated in \code{time}.}
  \item{unit}{Character containing the name of the units of observation in case of a panel structure. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), a single unit of observation is assumed.}
  \item{time}{Character containing the name of the time variable, which must be in numeric or date format. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), data are assumed to be temporally ordered.}
  \item{data}{Object of class \code{data.frame} containing the variables in \code{var.names}, \code{unit} and \code{time}.}
  \item{box.cox}{Named vector indicating the parameter of the Box-Cox transformation for variables in \code{x.names}. If \code{box.cox} is of length one and has no names, the same parameter is used for all variables in \code{x.names}. If a variable in \code{x.names} has no name in \code{box.cox}, value 1 is assumed, meaning no transformation. Value 0 of the parameter equates to the logarithmic transformation. Default value of \code{box.cox} is 1, meaning no transformation for all variables in \code{x.names}.}
  \item{ndiff}{Named vector with non-negative integer values indicating the number of differences for variables in \code{x.names}. If \code{ndiff} is of length one and has no names, the same number of differences is used for all variables in \code{x.names}. If a variable in \code{x.names} has no name in \code{ndiff}, value 0 is assumed, meaning no differencing. Default value of \code{ndiff} is 0, meaning no differencing for all variables in \code{x.names}.}
  \item{max.lag}{Non-negative integer representing the maximum lag length at which to perform the tests. If \code{NULL} (the default), it is taken as the squared root of the length of the time series.
%\code{trunc((length(x)-1)^(1/3))}
.}
}
\value{One list for each variable in \code{var.names}, each with three components:
\itemize{
\item{\code{statistic}: }{test statistic for both tests};
\item{\code{lag.selected}: }{lag length selected for both tests};
\item{\code{p.value}: }{p-value of both tests, which is a single value if \code{unit} is \code{NULL}, otherwise one value for each unit of observation plus another one indicating the combined p-value};
\item{\code{box.cox}: }{parameter of the Box-Cox transformation};
\item{\code{ndiff}: }{order of differencing}.
}
}
\note{The first order difference of logarithmic values (\code{box.cox}=0 and \code{ndiff}=1) provides the log returns, which approximate the proportional changes with respect to the previous time point.

If the parameter of the Box-Cox transformation is set to 0 for a non-negative variable, the transformation will not be applied. i.e., the parameter is set to 1, and a warning is returned.}
\details{The variable subjected to the test must be quantitative.

The null hypothesis of the ADF test is the presence of unit roots, while the null hypothesis of the KPSS test is the absence of unit roots.
Therefore, p-value higher than 0.05 for the ADF test or p-value lower than 0.05 for the KPSS test suggest the presence of unit roots and the need of further differencing.

The lag length at which to perform the tests is selected through backward elimination starting from the lag length specified in argument \code{max.lag}.

In case of a panel structure, p-values are combined according to the method by Demetrescu \emph{et al.} (2006).

Missing values internal to the time series are imputed through linear interpolation, otherwise they are deleted out.}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
M. Demetrescu, U. Hassler, and A. Tarcolea (2006). Combining significance of correlated statistics with application to panel data. \emph{Oxford Bulletin of Economics and Statistics}, 68(5), 647-663.
  
D. A. Dickey, and W. A. Fuller (1981). Likelihood ratio statistics for autoregressive time series with a unit root. \emph{Econometrica}, 49(4): 1057-1072.

D. Kwiatkowski, P. C. B. Phillips, P. Schmidt, and Y. Shin (1992). Testing the null hypothesis of stationarity against the alternative of a unit root. \emph{Journal of Econometrics}, 54(1-3): 159-178.

%S. Ng, and W. P. Perron (2001). Unit Root Tests in ARMA Models with Data-Dependent Methods for the Selection of the Truncation Lag. \emph{Journal of the American Statistical Association}, 90(429): 268-281.
}
\seealso{\link{preProcess}.}
\examples{
data(btc)
mydata <- btc[which(btc$Date>="2020-04-01"),]

# tests on time series in level ('box.cox'=1 by default):
#   there are some unit roots
unirootTest(var.names=c("BTC","DJA","IXIC","GSPC"), time="Date",
  data=mydata)
  
# tests on time series in log return: ok
unirootTest(var.names=c("BTC","DJA","IXIC","GSPC"), time="Date",
  data=mydata, box.cox=0, ndiff=1)
}