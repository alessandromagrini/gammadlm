\name{tsEM}
\alias{tsEM}
\title{Imputation of multivariate time series data}
\description{Imputation of a multivariate time series based on the Expectation-Maximization algorithm. Data may have a panel structure, even unbalanced.}
\usage{tsEM(var.names, unit=NULL, time=NULL, data, nlags=NULL, tol=1e-4, maxit=1000, quiet=FALSE)}
\arguments{
  \item{var.names}{Character vector including the name of the variables in the imputation model. See 'Details'.}
  \item{unit}{Character containing the name of the variable that identifies the units of observation. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), a single unit of observation is assumed.}
  \item{time}{Character containing the name of the time variable, which must be in numeric or date format. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), data are assumed to be temporally ordered.}
  \item{data}{Object of class \code{data.frame} containing the variables in \code{var.names}, \code{unit} and \code{time}.}
  \item{nlags}{Non-negative integer number indicating the lag order of the imputation model. It will be automatically bounded in order to keep the number of free parameters lower than half the sample size. If \code{NULL} (the default), it will be automatically determined based on empirical correlations among variables in \code{var.names}.}
  \item{tol}{Positive number indicating the tolerance. Default is 1e-4.}
  \item{maxit}{Positive integer indicating the number of iterations. Default is 1000.}
  \item{quiet}{Logical value indicating whether prompt messages should be suppressed. Default is \code{FALSE}.}
}
\value{The object provided to argument \code{data} where variables in  \code{var.names} have been imputed.
}
%\note{}
\details{The imputation model is a nonstationary vector autoregressive (VAR) model, where quantitative and categorical variables provided to argument \code{var.names} are considered, respectively, endogenous and exogenous.

Missing values of quantitative variables are replaced by their conditional expected value under the imputation model, which is estimated through the Expectation-Maximization algorithm (Dempster et al., 1977).

The optimal Box-Cox transformation (Box & Cox, 1964) for each quantitative variable is determined automatically, and the transformation is inverted before returning the imputed data.}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
Box, G.E.P., Cox, D.R. (1964). An analysis of transformations. \emph{Journal of the Royal Statistical Society}, Series B (Methodological), 26(2): 211-252. DOI: 10.1111/j.2517-6161.1964.tb00553.x

Dempster, A.P., Laird, N.M., Rubin, D.B. (1977). Maximum likelihood from incomplete data via the EM algorithm. \emph{Journal of the Royal Statistical Society}n Series B (Methodological), 39(1): 1-38. DOI: 10.1111/j.2517-6161.1977.tb01600.x
}
%\seealso{\link{gammadlm}.}
\examples{
data(btc)
mydata <- btc[which(btc$Date>="2020-04-01"),]

# add some missing values
mydata0 <- mydata
mydata0$BTC[sample(1:nrow(mydata),30)] <- NA
mydata0$DJA[sample(1:nrow(mydata),30)] <- NA 
mydata0$GSPC[sample(1:nrow(mydata),30)] <- NA 
mydata0$IXIC[sample(1:nrow(mydata),30)] <- NA 
summary(mydata0)
plot(ts(mydata0[,c("BTC","DJA","IXIC","GSPC")]), main="")

mydataI <- tsEM(var.names=c("BTC","DJA","IXIC","GSPC"), time="Date",
  data=mydata0)
summary(mydataI)
plot(ts(mydataI[,c("BTC","DJA","IXIC","GSPC")]), main="")
}