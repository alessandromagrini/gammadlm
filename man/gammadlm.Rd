\name{gammadlm}
\alias{gammadlm}
\title{Estimation of a Gamma distributed-lag model}
\description{Maximum likelihood estimation of a Gamma distributed-lag model with multiple explanatory variables using the hill climbing algorithm. Data may have a panel structure, even unbalanced.}
\usage{gammadlm(y.name, x.names, z.names=NULL, unit=NULL, time=NULL, data,
  offset=rep(0,length(x.names)), box.cox=1, ndiff=0, add.intercept=TRUE,
  control=list(nstart=NULL, delta.lim=NULL, lambda.lim=NULL, peak.lim=NULL,
  length.lim=NULL), quiet=FALSE)}
\arguments{
  \item{y.name}{Character including the name of the response variable, that must be a quantitative variable. If a vector with length greater than 1 is provided, only the first element is considered.}  
  \item{x.names}{Character vector of length 1 or greater including the names of the explanatory variables with lags, that must be quantitative variables. If the name of the response variable is indicated in \code{x.names}, then the autoregressive lag distribution will be estimated.}
  \item{z.names}{Character vector including the names of the explanatory variables without lags (optional). They may be either quantitative or qualitative variables. If \code{NULL} (the default), no explanatory variable without lags is included in the model.}
  \item{unit}{Character containing the name of the variable that identifies the units of observation. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), a single unit of observation is assumed.}
  \item{time}{Character containing the name of the time variable, which must be in numeric or date format. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), data are assumed to be temporally ordered.}
  \item{data}{Object of class \code{data.frame} containing the variables in the model. Variables in \code{y.name}, \code{x.names}, \code{z.names}, \code{unit} and \code{time} cannot contain missing values.}
  \item{offset}{Named vector with non-negative real values indicating the offset for variables in \code{x.names}. If \code{offset} is of length one and has no names, the same offset is used for all variables in \code{x.names}. If a variable in \code{x.names} has no name in \code{offset} or a negative value is provided, value 0 is assumed. If the name of the response variable is indicated in \code{x.names} and a value less than 1 is provided, value 1 is assumed. Default value of \code{offset} is 0.}
  \item{box.cox}{Named vector with non-negative real values indicating the parameters of the Box-Cox transformation (Box & Cox, 1964) for variables in \code{x.names}. If \code{box.cox} has no names and length greater than one, the same ordering as in \code{x.names} is assumed. If \code{box.cox} has no names and length equal to one, the same parameter is used for all variables in \code{x.names}. Value 0 of the parameter equates to the logarithmic transformation, while value 1 means no transformation. Default is 1 for all variables in \code{x.names}.}
  \item{ndiff}{Named vector with non-negative integer values indicating the number of differences for variables in \code{x.names}. If \code{ndiff} has no names and length greater than one, the same ordering as in \code{x.names} is assumed. If \code{ndiff} has no names and length equal to one, the same number of differences is used for all variables in \code{x.names}. Value 0 means no differencing. Default is 0 for all variables in \code{x.names}.}
  \item{add.intercept}{Logical value indicating whether the intercept should be included in the model. Default is \code{TRUE}. See 'Details'.}
  \item{control}{A list including control options for the hill climbing algorithm.
   \itemize{
    \item{\code{nstart}: }{positive integer value indicating the number of restarts. If equal to 1 or \code{NULL} (the default), shape parameters are initialized for one explanatory variable at a time based on OLS-grid search.}
    %\item{\code{grid.by}: }{positive value no greater than 0.1, indicating the increment in grid search. Default is 0.05.}
    \item{\code{delta.lim}: }{a named list with one component for each variable in \code{x.names}, that must be either a numerical vector of length 2 indicating the minimum and the maximum value of \eqn{\delta}, or a numerical value indicating the exact value of \eqn{\delta}. If there is no component in \code{delta.lim} for a certain variable in \code{x.names}, then the theoretical range [0,1) is assumed for \eqn{\delta}.
%If, for an explanatory variable, minimum and maximum values coincide, then \eqn{\delta} parameter is kept fixed at their values.
}
    \item{\code{lambda.lim}: }{the same as \code{delta.lim}, but it is about \eqn{\lambda} parameters.}
    \item{\code{peak.lim}: }{the same as \code{delta.lim}, but it is about the peak of the lag distributions.}
    \item{\code{length.lim}: }{the same as \code{delta.lim}, but it is about the 99.9th percentile of the lag distributions.}
    }}
  \item{quiet}{Logical value indicating whether prompt messages should be suppressed. Default is \code{FALSE}.}
  }
\value{An object of class \code{lm} and \code{gammadlm}, including all the components of an object of class \code{lm} plus the following components:
\itemize{
  \item{\code{par}: }{matrix including the shape parameters for each variable in \code{x.names} (by column);}
  \item{\code{offset}: }{vector including the offset of the lag distributions;}
  \item{\code{add.intercept}: }{value passed to argument \code{add.intercept};}
  \item{\code{local.max}: }{list including the best model for each restart.}
  \item{\code{variables}: }{list including the names of the variables provided to arguments \code{y.name}, \code{x.names}, \code{z.names}, \code{unit} and \code{time}.}
  \item{\code{control}: }{list including control options used in parameter estimation;}
  \item{\code{unit.id}: }{list including the row names of the observations for each unit. \code{NULL} if \code{unit} is \code{NULL};}
  \item{\code{box.cox}: }{named vector including the parameter  of the Box-Cox trasnformation for each time series;}
  \item{\code{ndiff}: }{named vector including the order fo differencing for each time series;}
\item{\code{lag.order}: }{lag orders of residuals and squared residuals estimated by minimizing the Bayesian Information Criterion (see 'Note'). If \code{unit} is not \code{NULL}, lag orders will be provided as a matrix with one row for each unit of observation;}
  \item{\code{data.orig}: }{data.frame including the data before Box-Cox transformation and differencing.}
  \item{\code{data.used}: }{data.frame including the data used for parameter estimation.}

  }
}
\details{
Second-order stationarity of all time series (expected value and autocorrelation function independent of time) is a basic assumption of the model, that is guaranteed if no time series contains unit roots.
Box-Cox transformation (argument \code{box.cox}) and/or differencing (argument \code{ndiff}) can be exploited in order to achieve second-order stationarity of each time series.

If the response variable is not differenced, the intercept captures the effect of time-invariant factors. Otherwise, time-invariant factors are canceled out and the intercept represents the coefficient of a linear trend (drift). When \code{unit} is not \code{NULL} (panel data), unit-specific intercepts are included if the response variable is not differenced, otherwise a unique intercept (drift) for all units is included. 

All S3 methods for class \code{lm} are also available for class \code{gammadlm}:

- \code{plot}: to display the estimated lag distributions;

- \code{summary}: to see the summary of parameter estimation;

- \code{residuals}: to get the residuals;

- \code{fitted.values}: to get the fitted values;

- \code{predict}: to obtain predictions.
%performs \eqn{h} step ahead forecast.

Furthermore, method \link{lagCoef} can be used to see the estimated dynamic multipliers.
}
\note{
%Before calling the function \link{gammadlm}, the user is strongly recommended to check the absence of unit roots in each time series through the function \link{unirootTest}, and to apply differencing through the function \link{preProcess} to the ones containing unit roots.
The first order difference of logarithmic values (\code{box.cox}=0 and \code{ndiff}=1) provides the log returns, which approximate the relative changes with respect to the previous time point.

If a variable contains negative values, the Box-Cox transformation will be not applied and a warning is returned.
If \code{box.cox}=0 (logarithmic transformation) and the variable contains positive and null values, argument \code{box.cox} will be set to 0.5 (squared root transformation) and a warning is returned.

If the number of differencing exceeds \code{n-5}, where \code{n} is the sample size, differencing will be not applied and a warning is returned.

The absence of unit roots in the residuals is automatically checked and a warning is returned in case of failure.

When the \code{summary} method is called on an object of class \code{gammadlm}, the Heteroskedasticity and Autocorrelation Consistent (HAC, Newey & West, 1987) estimator of the covariance matrix of least squares estimates is applied, based on the lag order of residuals and squared residuals (see \code{lag.order} in 'Values').
The same holds for the \code{confint} method.
}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
G. E. P. Box, and D. R. Cox (1964). An analysis of transformations. \emph{Journal of the Royal Statistical Society}, Series B (Methodological), 26(2): 211-252. DOI: 10.1111/j.2517-6161.1964.tb00553.x

A. Magrini (2022). A hill climbing algorithm for maximum likelihood estimation of the Gamma distributed-lag model with multiple explanatory variables. \emph{Austrian Journal of Statistics}, 51(2): 40-46. DOI: 10.17713/ajs.v51i2.1244

W. K. Newey, K. D. West (1987). A simple, positive semi-definite, heteroskedasticity and autocorrelation consistent covariance matrix. \emph{Econometrica}, 55(3): 703-708
}
\seealso{\link{unirootTest}; \link{tsEM}; \link{lagCoef}; \link{plot.gammadlm}.
}
\examples{
# load data on Bitcoin and US stock market price
data(btc)
mydata <- btc[which(btc$Date>="2020-04-01"),]

# model with pre-specified shape parameters (simple OLS)
#   on data in logarithmic differences (yearly relative changes)
dval <- list(DJA=0.85,IXIC=0.75,GSPC=0.55)
lval <- list(DJA=0.05, IXIC=0.35,GSPC=0.45)
m0 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), time="Date", data=mydata,
  box.cox=0, ndiff=1, control=list(delta.lim=dval, lambda.lim=lval))
summary(m0)  ## summary of estimation

# hill climbing algorithm without random restarts
#   (initialization based on independent OLS-grid searches)
m1 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydata,
  box.cox=0, ndiff=1)
summary(m1)

# add constraints: peak>=1 and 3<=length<=15
pklim <- list(DJA=c(1,Inf),IXIC=c(1,Inf),GSPC=c(1,Inf))
lenlim <- list(DJA=c(3,15),IXIC=c(3,15),GSPC=c(3,15))
m2 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydata,
  box.cox=0, ndiff=1, control=list(peak.lim=pklim, length.lim=lenlim))
summary(m2)

\dontrun{
# 50 random restarts without constraints
set.seed(100)
m3 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydata,
  box.cox=0, ndiff=1, control=list(nstart=50))
summary(m3)  ## better fit than m1

# 50 random restarts with constraints: NOT RUN
set.seed(100)
m4 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydata,
  box.cox=0, ndiff=1, control=list(peak.lim=pklim, length.lim=lenlim, nstart=50))
summary(m4)  ## better fit than m2}


###  Example with panel data  ###

\dontrun{
# load data on farm performance and subsidies in EU countries, 2004-2019
data(fadn)

# considered variables:
#  Y: productivity
#  X1-X4: subsidies by typology
#  Z1: utilized agricultural area (dimensional class),
#  Z2: total output (economic class)
y_name <- "TFP"
x_name <- c("Subs_prod","Subs_inv","Subs_rur","Subs_dec")
z_name <- c("Land","Output_total")

# model on data in differences: in this case, the intercept
#   represents the coefficient of a linear trend (drift)
m_fadn <- gammadlm(y.name=y_name, x.names=x_name, z.names=z_name,
  unit="Country", time="Year", data=fadn, box.cox=1, ndiff=1,
  control=list(peak.lim=c(1,Inf), length.lim=c(3,10)))
summary(m_fadn)}
}