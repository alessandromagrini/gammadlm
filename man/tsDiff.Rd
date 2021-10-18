\name{tsDiff}
\alias{tsDiff}
\title{Time series differencing}
\description{Application of differencing to a multivariate time series, eventually structured as a panel.}
\usage{tsDiff(var.names=NULL, panel=NULL, time=NULL, data, box.cox=1, ndiff=0)}
\arguments{
  \item{var.names}{Character vector including the name of the variables to be differenced. If \code{NULL} (the default), all the quantitative variables in the dataset provided to argument \code{data} will be differenced, with the exception of the variable indicated in \code{time}.}
  \item{panel}{Character containing the name of the units of observation in case of a panel structure. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), a single unit of observation is assumed.}
  \item{time}{Character containing the name of the time variable, which must be in numeric or date format. If a vector with length greater than 1 is provided, only the first element is considered. If \code{NULL} (the default), data are assumed to be temporally ordered.}
  \item{data}{Object of class \code{data.frame} containing the variables in \code{var.names}, \code{panel} and \code{time}.}
  \item{box.cox}{Named vector indicating the parameter of the Box-Cox transformation for variables in \code{x.names}. If \code{box.cox} is of length one and has no names, the same parameter is used for all variables in \code{x.names}. If a variable in \code{x.names} has no name in \code{box.cox}, value 1 is assumed, meaning no transformation. Value 0 of the parameter equates to the logarithmic transformation. Default value of \code{box.cox} is 1, meaning no transformation for all variables in \code{x.names}.}
  \item{ndiff}{Named vector with non-negative integer values indicating the number of differences for variables in \code{x.names}. If \code{ndiff} is of length one and has no names, the same number of differences is used for all variables in \code{x.names}. If a variable in \code{x.names} has no name in \code{ndiff}, value 0 is assumed, meaning no differencing. Default value of \code{ndiff} is 0, meaning no differencing for all variables in \code{x.names}.}
}
\value{The object provided to argument \code{data} where variables in  \code{var.names} have been transformed and/or differenced, and the following two attributes are added:

- \code{inits}: if \code{panel}=\code{NULL}, named vector indicating the first observed value for each variable in \code{x.names}, otherwise a named matrix indicating the first observed value for each unit of observation (by row) and each variable in \code{x.names} (by column);

- \code{box.cox}: named vector indicating the parameter of the Box-Cox transformation for each variable in \code{x.names};

- \code{ndiff}: named vector indicating the number of differences for each variable in \code{x.names}.
}
\note{The first order difference of logarithmic values (\code{box.cox}=0 and \code{ndiff}=1) provides the log returns, which approximate the proportional changes with respect to the previous time point.

If the parameter of the Box-Cox transformation is set to 0 for a non-negative variable, the transformation will not be applied (parameter equal to 1) and a warning is returned.}
%\details{}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\link{adfTest}.}
\examples{
data(USstock)
mydata <- USstock[which(USstock$Date>="2020-04-01"),]

# first order differencing ('box.cox'=1 by default)
mydataD <- tsDiff(var.names=c("BTC","DJA","IXIC","GSPC"), time="Date",
  data=mydata, ndiff=1)
summary(mydataD)
plot(ts(mydataD[,c("BTC","DJA","IXIC","GSPC")]), main="")

# setting box.cox=0 and ndiff=1 produces the log returns
mydataLR <- tsDiff(var.names=c("BTC","DJA","IXIC","GSPC"), time="Date",
  data=mydata, box.cox=0, ndiff=1)
summary(mydataLR)
plot(ts(mydataLR[,c("BTC","DJA","IXIC","GSPC")]), main="")

# same result by omitting 'var.names'
mydataLR2 <- tsDiff(time="Date", data=mydata, box.cox=0, ndiff=1)
summary(mydataLR2)
plot(ts(mydataLR2[,c("BTC","DJA","IXIC","GSPC")]), main="")
}