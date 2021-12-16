\name{gammadlm}
\alias{gammadlm}
\title{Estimation of a Gamma distributed-lag model}
\description{Maximum likelihood estimation of a Gamma distributed-lag model with multiple explanatory variables using the hill climbing algorithm. A panel structure can be taken into account.}
\usage{gammadlm(y.name, x.names, z.names=NULL, unit=NULL, time=NULL, data,
  offset=rep(0,length(x.names)), control=list(nstart=50, grid.by=0.05,
  delta.lim=NULL, lambda.lim=NULL, peak.lim=NULL, length.lim=NULL), quiet=FALSE)}
\arguments{
  \item{y.name}{Character including the name of the response variable, that must be a quantitative variable. If a vector with length greater than 1 is provided, only the first element is considered.}  
  \item{x.names}{Character vector of length 1 or greater including the names of the explanatory variables with lags, that must be quantitative variables. If the name of the response variable is indicated in \code{x.names}, then the autoregressive lag distribution will be estimated.}
  \item{z.names}{Character vector including the names of the explanatory variables without lags (optional). They may be either quantitative or qualitative variables. If \code{NULL} (the default), no explanatory variable without lags is included in the model.}
  \item{unit}{Character containing the name of the units of observation in case of a panel structure. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), a single unit of observation is assumed.}
  \item{time}{Character containing the name of the time variable, which must be in numeric or date format. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), data are assumed to be temporally ordered.}
  \item{data}{Object of class \code{data.frame} containing the variables in the model. Variables in \code{y.name}, \code{x.names}, \code{z.names}, \code{unit} and \code{time} cannot contain missing values.}
  %\item{par}{} %%%
  \item{offset}{Numerical vector indicating, in order, the offset for each variable in \code{x.names}. Default is 0. If the name of the response variable is indicated in \code{x.names} and its offset is less than 1, then it will be set automatically to 1.}
  \item{control}{A list including control options for the hill climbing algorithm.
   \itemize{
    \item{\code{nstart}: }{positive integer value indicating the number of restarts. If equal to 1, then all the shape parameters are initialized to value 0. Default is 50.}
    %\item{\code{max.try}: }{}
    %\item{\code{max.start}: }{} 
    \item{\code{grid.by}: }{positive value no greater than 0.1, indicating the increment in grid search. Default is 0.05.}
    \item{\code{delta.lim}: }{a named list with one component for each variable in \code{x.names}, that must be either a numerical vector of length 2 indicating the minimum and the maximum value of \eqn{\delta}, or a numerical value indicating the exact value of \eqn{\delta}. If there is no component in \code{delta.lim} for a certain variable in \code{x.names}, then the theoretical range [0,1) is assumed for \eqn{\delta}.
%If, for an explanatory variable, minimum and maximum values coincide, then \eqn{\delta} parameter is kept fixed at their values.
}
    \item{\code{lambda.lim}: }{the same as \code{delta.lim}, but it is about \eqn{\lambda} parameters.}
    \item{\code{peak.lim}: }{the same as \code{delta.lim}, but it is about the peak of the lag distributions.}
    \item{\code{length.lim}: }{the same as \code{delta.lim}, but it is about the 99.9th percentile of the lag distributions.}
    }}
  \item{quiet}{Logical value indicating whether prompt messages should be displayed during the execution. Default is \code{TRUE}.}
  }
\value{An object of class \code{lm} and \code{gammadlm}, including all the components of an object of class \code{lm} plus the following components:
\itemize{
  \item{\code{offset}: }{vector including the offset of the lag distributions;}
  \item{\code{par}: }{matrix including the shape parameters for each variable in \code{x.names} (by column);}
  \item{\code{variables}: }{list including the names of the variables provided to arguments \code{y.name}, \code{x.names} and \code{z.names}.}
  \item{\code{unit.id}: }{list including the row names of the observations for each unit.}
  \item{\code{data}: }{data.frame including the data used for parameter estimation.}
  \item{\code{local.max}: }{list including all the models fitted at each restart.}
  }
}
\details{
All S3 methods for class \code{lm} are also available for class \code{gammadlm}.
Furthermore, method \link{lagCoef} can be used to see the estimated dynamic coefficients, and method \code{plot} to display the estimated lag distributions.
%and method \code{residuals} to obtain graphical model diagnostics.
}
\note{
Weak stationarity of all time series (expected value and autocorrelation function independent of time) is a basic assumption of the model.
However, it is sufficient that all time series do not contain unit roots to get valid results (Granger \& Newbold, 1974).
Before calling the function \link{gammadlm}, the user is strongly recommended to check the absence of unit roots in each time series through the function \link{unirootTest}, and to apply differencing through the function \link{tsDiff} to each time series with unit roots.

Function \link{gammadlm} automatically checks the absence of unit roots in the residuals and returns a warning in case of failure.

When the \code{summary} method is called on an object of class \code{gammadlm}, the order of auto-correlation of the residuals is estimated based on the Bayesian Information Criterion. If it is greater than 0, then the Heteroskedasticity and Autocorrelation Consistent (HAC, Newey & West, 1987) estimator of the covariance matrix of least squares estimates is applied to get robust standard errors.
The same holds for the \code{confint} method.
}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
%D. A. Dickey, and W. A. Fuller (1981). Likelihood ratio statistics for autoregressive time series with a unit root. \emph{Econometrica}, 49(4): 1057-1072.

C. W. J. Granger, and P. Newbold (1974). Spurious regressions in econometrics. \emph{Journal of Econometrics}, 2(2): 111-120.

A. Magrini (2022). A hill climbing algorithm for maximum likelihood estimation of the Gamma distributed-lag model with multiple explanatory variables. \emph{Austrian Journal of Statistics}, 51(2): 40-46.

W. K. Newey, K. D. West (1987). A simple, positive semi-definite, heteroskedasticity and autocorrelation consistent covariance matrix. \emph{Econometrica}, 55(3): 703-708
}
\seealso{\link{unirootTest}; \link{lagCoef}; \link{plot.gammadlm}.
%\link{residuals.gammadlm}.
}
\examples{
data(USstock)
mydata <- USstock[which(USstock$Date>="2020-04-01"),]
mydataLR <- tsDiff(time="Date", data=mydata, box.cox=0, ndiff=1)

## estimation with fixed values of delta and lambda parameters: 1-step OLS
dval <- list(DJA=0.85,IXIC=0.75,GSPC=0.55)
lval <- list(DJA=0.5, IXIC=0.35,GSPC=0.45)
mod <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR,
  control=list(delta.lim=dval, lambda.lim=lval))
summary(mod)  ## summary of estimation

## estimation through the hill climbing algorithm: NOT RUN
##
## * no constraints with 50 random restarts (by default)
#set.seed(100)
#m1 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR)
#summary(m1)
##
## * no contraints with 100 random restarts
#set.seed(100)
#m1a <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR,
#  control=list(nstart=100))
#summary(m1a)
##
## * constraints: peak>=1 and 3<=length<=10
#pklim <- list(DJA=c(1,Inf),IXIC=c(1,Inf),GSPC=c(1,Inf))
#lenlim <- list(DJA=c(3,10),IXIC=c(3,10),GSPC=c(3,10))
#set.seed(100)
#m2 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR,
#  control=list(peak.lim=pklim, length.lim=lenlim, nstart=100))
#summary(m2)
}