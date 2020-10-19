\name{gammadlm}
\alias{gammadlm}
\title{Estimation of a Gamma distributed-lag model}
\description{Maximum likelihood estimation of a Gamma distributed-lag model with multiple explanatory variables using hill climbing algorithm.}
\usage{gammadlm(y.name, x.names, z.names=NULL, time.name=NULL, data,
  offset=rep(0,length(x.names)), control=list(nstart=50, grid.by=0.05,
  delta.lim=NULL, lambda.lim=NULL, peak.lim=NULL, length.lim=NULL), quiet=FALSE)}
\arguments{
  \item{y.name}{Character including the name of the response variable, that must be a quantitative variable. If a vector with length greater than 1 is provided, only the first element is considered.}  
  \item{x.names}{Character vector of length 1 or greater including the names of the explanatory variables with lags, that must be quantitative variables. If the name of the response variable is indicated in \code{x.names}, then the autoregressive lag distribution will be estimated.}
  \item{z.names}{Character vector including the names of the explanatory variables without lags (optional). They may be either quantitative or qualitative variables. If \code{NULL} (the default), no explanatory variable without lags is included in the model.}
  \item{time.name}{Character including the name of the time variable, which must be in numeric or date format. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), data are assumed to be temporally ordered.}
  \item{data}{Object of class \code{data.frame} containing the variables in the model. Variables in \code{y.name}, \code{x.names}, \code{z.names} and \code{time.name} cannot contain missing values.}
  %\item{par}{} %%%
  \item{offset}{Numerical vector indicating, in order, the offset for each variable in \code{x.names}. Default is 0. If the name of the response variable is indicated in \code{x.names} and its offset is less than 1, then it will be set automatically to 1.}
  \item{control}{A list including control options for the hill climbing algorithm.
   \itemize{
    \item{\code{nstart}: }{positive integer value indicating the number of restarts. Default is 50.}
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
  \item{\code{data}: }{data.frame including the data used for parameter estimation.}
  }
}
\details{
Let \eqn{Y} be the response variable and \eqn{X_1,\ldots,X_J} be \eqn{J} explanatory variables.
Also, let \eqn{y_t} and \eqn{x_{j,t}} be, respectively, the value of \eqn{Y} and of \eqn{X_j} (\eqn{j=1,\ldots,J}) at time \eqn{t}.
The Gamma distributed-lag model explaining \eqn{Y} from \eqn{X_1,\ldots,X_J} is defined as:
\deqn{y_t = \alpha+\sum_{j=1}^J\theta_j\sum_{k=0}^\infty w_{j,k}x_{j,t-k}+\varepsilon_t}
\deqn{w_{j,k} = \frac{(k+1-\eta_j)^{\frac{\delta_j}{1-\delta_j}}\lambda_j^{k-\eta_j}}{\sum_{l=0}^\infty(l+1-\eta_j)^{\frac{\delta_j}{1-\delta_j}}\lambda_j^{l-\eta_j}}}
where:
\itemize{
\item{\eqn{\alpha} is the intercept};
\item{\eqn{\theta_j} is the scale parameter of the lag distribution of \eqn{X_j}}; \item{\eqn{0\leq\delta_j<1} and \eqn{0\leq\lambda_j<1} are the shape parameters of the lag distribution of \eqn{X_j};}
\item{\eqn{\eta_j} is the offset of the lag distribution of \eqn{X_j}, typically set to 0;}
\item{\eqn{\varepsilon_t} is the random error at time \eqn{t}}.
}
All S3 methods for class \code{lm} are also available for class \code{gammadlm}.
Furthermore, method \link{lagCoef} can be used to obtain dynamic coefficients, method \code{plot} to obtain graphics for the estimated lag distributions, and method \code{residuals} to obtain graphical model diagnostics.
}
\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
A. Magrini. A hill climbing algorithm for maximum likelihood estimation of the Gamma distributed-lag model with multiple explanatory variables. \emph{Under review}.

J. M. Alston, M. A. Andersen, J. S. James, P. G. Pardey (2011). The economic returns to U.S. public agricultural research. \emph{American Journal of Agricultural Economics}, 93(5), 1257-1277.
}
\seealso{\link{lagCoef}; \link{plot.gammadlm}; \link{residuals.gammadlm}.}
\examples{
data(USstock)
mydata <- USstock[which(USstock$Date>="2020-04-01"),]
mydataLR <- tsDiff(time.name="Date", data=mydata, ndiff=1, log=TRUE)

## estimation with fixed values of delta and lambda parameters: 1-step OLS
dval <- list(DJA=0.85,IXIC=0.75,GSPC=0.55)
lval <- list(DJA=0.5, IXIC=0.35,GSPC=0.45)
mod <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR,
  control=list(delta.lim=dval, lambda.lim=lval))
summary(mod)  ## summary of estimation

## estimation through hill climbing algorithm: NOT RUN
##
## * no constraints
#set.seed(100)
#m1 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR)
#summary(m1)
##
## * same but with 100 random restarts
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