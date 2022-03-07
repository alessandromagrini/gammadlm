\name{preProcess}
\alias{preProcess}
\title{Preprocessing of time series data}
\description{Application of the Box-Cox transformation and/or differencing to a multivariate time series, eventually structured as a balanced or unbalanced panel.}
\usage{preProcess(var.names=NULL, unit=NULL, time=NULL, data, box.cox=1, ndiff=0,
  imputation=TRUE, em.control=list(nlags=NULL, tol=1e-4, maxit=1000, quiet=FALSE))}
\arguments{
  \item{var.names}{Character vector including the name of the variables to be differenced. If \code{NULL} (the default), all the quantitative variables in the dataset provided to argument \code{data} will be differenced, with the exception of the variable indicated in \code{time}.}
  \item{unit}{Character containing the name of the units of observation in case of a panel structure. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), a single unit of observation is assumed.}
  \item{time}{Character containing the name of the time variable, which must be in numeric or date format. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), data are assumed to be temporally ordered.}
  \item{data}{Object of class \code{data.frame} containing the variables in \code{var.names}, \code{unit} and \code{time}.}
  \item{box.cox}{Named vector with non-negative real values indicating the parameters of the Box-Cox transformation (Box & Cox, 1964) for variables in \code{x.names}. If \code{box.cox} is of length one and has no names, the same parameter is used for all variables in \code{x.names}. If a variable in \code{x.names} has no name in \code{box.cox}, value 1 is assumed, meaning no transformation. Value 0 of the parameter equates to the logarithmic transformation. Default value of \code{box.cox} is 1, meaning no transformation for all variables in \code{x.names}.}
  \item{ndiff}{Named vector with non-negative integer values indicating the number of differences for variables in \code{x.names}. If \code{ndiff} is of length one and has no names, the same number of differences is used for all variables in \code{x.names}. If a variable in \code{x.names} has no name in \code{ndiff}, value 0 is assumed, meaning no differencing. Default value of \code{ndiff} is 0, meaning no differencing for all variables in \code{x.names}.}
  \item{imputation}{Logical value indicating whether imputation of missing values should be performed (see 'Details'). Default is \code{TRUE}.}
  \item{em.control}{List including control options for the EM algorithm (see 'Details') with the following components:
  \itemize{
  \item{\code{nlags}: }{Non-negative integer value indicating the number of lags to consider. If \code{NULL} (the default), it is taken as \code{trunc((n-1)^(1/3))}, where \code{n} is the sample size, otherwise it is right-truncated at \code{trunc(n*2/3)}.}
  \item{\code{tol}: }{Positive number indicating the tolerance. Default is 0.001.}
  \item{\code{maxit}: }{Positive integer indicating the number of iterations. Default is 1000.}
  \item{\code{quiet}: }{Logical value indicating whether prompt messages should be suppressed. Default is \code{FALSE}.}
  }
  Ignored if \code{imputation} is \code{FALSE}.}
}
\value{The object provided to argument \code{data} where variables in  \code{var.names} have been transformed and/or differenced, and the following attributes are added:

%- \code{inits}: if \code{unit}=\code{NULL}, named vector indicating the first observed value for each variable in \code{x.names}, otherwise a named matrix indicating the first observed value for each unit of observation (by row) and each variable in \code{x.names} (by column);
%
- \code{box.cox}: named vector indicating the parameter of the Box-Cox transformation for each variable in \code{x.names};

- \code{ndiff}: named vector indicating the number of differences for each variable in \code{x.names}.
}
\note{The first order difference of logarithmic values (\code{box.cox}=0 and \code{ndiff}=1) provides the log returns, which approximate the proportional changes with respect to the previous time point.

If a variable contains negative values, the Box-Cox transformation will be not applied and a warning is returned.

If the number of differencing exceeds \code{n-5}, where \code{n} is the sample size, differencing will be not applied and a warning is returned.}
\details{
Imputation of missing values is performed after the Box-Cox transformation (Box & Cox, 1964) and differencing.
A vector autoregressive model is assumed where the expectation of missing values is computed through the Expectation-Maximization (EM, Dempster et al., 1977) algorithm.
}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
G. E. P. Box, and D. R. Cox (1964). An analysis of transformations. \emph{Journal of the Royal Statistical Society}, Series B (Methodological), 26(2): 211-252. DOI: 10.1111/j.2517-6161.1964.tb00553.x

A. P. Dempster, N. M. Laird, and D. B. Rubin (1977). Maximum likelihood from incomplete data via the EM algorithm. \emph{Journal of the Royal Statistical Society}n Series B (Methodological), 39(1): 1-38. DOI: 10.1111/j.2517-6161.1977.tb01600.x
}
\seealso{\link{unirootTest}.}
\examples{
data(btc)
mydata <- btc[which(btc$Date>="2020-04-01"),]

# first order differencing ('box.cox'=1 by default)
mydataD <- preProcess(var.names=c("BTC","DJA","IXIC","GSPC"), time="Date",
  data=mydata, ndiff=1)
summary(mydataD)
plot(ts(mydataD[,c("BTC","DJA","IXIC","GSPC")]), main="")

# setting box.cox=0 and ndiff=1 produces the log returns
mydataLR <- preProcess(var.names=c("BTC","DJA","IXIC","GSPC"), time="Date",
  data=mydata, box.cox=0, ndiff=1)
summary(mydataLR)
plot(ts(mydataLR[,c("BTC","DJA","IXIC","GSPC")]), main="")

# same result by omitting 'var.names'
mydataLR2 <- preProcess(time="Date", data=mydata, box.cox=0, ndiff=1)
summary(mydataLR2)
plot(ts(mydataLR2[,c("BTC","DJA","IXIC","GSPC")]), main="")
}