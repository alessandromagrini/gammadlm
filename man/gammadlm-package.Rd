\name{gammadlm-package}
\alias{gammadlm-package}
\docType{package}
\title{The Gamma Distributed-Lag Model with Multiple Explanatory Variables}
\description{Maximum likelihood estimation and inference for the Gamma distributed-lag model with multiple explanatory variables. A panel structure can be taken into account.}
\details{
  \tabular{ll}{
  Package: \tab gammadlm\cr
  Type: \tab Package\cr
  Version: \tab 0.1.1\cr
  Date: \tab 2021-12-11\cr
  License: \tab GPL-2
  }
Let \eqn{Y} be the response variable and \eqn{X_1,\ldots,X_J} be \eqn{J} explanatory variables.
Also, let \eqn{y_t} and \eqn{x_{j,t}} be, respectively, the value of \eqn{Y} and of \eqn{X_j} (\eqn{j=1,\ldots,J}) observed at time \eqn{t}. %on unit \eqn{i}
Under the assumption that the time series of \eqn{Y} and of \eqn{X_1,\ldots,X_J} are all weakly stationary (i.e., expected value and autocorrelation function independent of time), the Gamma distributed-lag model explaining \eqn{Y} from \eqn{X_1,\ldots,X_J} is defined as:
\deqn{y_t = \alpha+\sum_{j=1}^J\sum_{k=0}^\infty\beta_{j,k}(\theta_j,\delta_j,\lambda_j,\eta_j)\hspace{0.07cm}x_{j,t-k}+\varepsilon_t}
\deqn{\beta_{j,k}(\theta_j,\delta_j,\lambda_j,\eta_j)=\theta_j \hspace{0.07cm} w_{j,k}(\delta_j,\lambda_j,\eta_j)}
\deqn{w_{j,k}(\delta_j,\lambda_j,\eta_j)=\frac{(k+1-\eta_j)^{\frac{\delta_j}{1-\delta_j}}\lambda_j^{k-\eta_j}}{\sum_{l=0}^\infty(l+1-\eta_j)^{\frac{\delta_j}{1-\delta_j}}\lambda_j^{l-\eta_j}}}
where:
\itemize{
\item{\eqn{\alpha} is the intercept};
\item{\eqn{\beta_{j,k}} is the dynamic coefficient for \eqn{X_j} at time lag \eqn{k}, equal to the scale parameter \eqn{\theta_j} times the weight \eqn{w_{j,k}}.
The set \eqn{\{w_{j,k}: k=0,1,\ldots,\infty \}} includes the weights for \eqn{X_j} and is defined by the shape parameters \eqn{0\leq\delta_j<1} and \eqn{0\leq\lambda_j<1} and the offset \eqn{\eta_j} (typically set to 0).
The set \eqn{\{ \beta_{j,k}: k=0,1,\ldots,\infty \}} is called \emph{lag distribution} of \eqn{X_j};}
\item{\eqn{\varepsilon_t} is the random error at time \eqn{t}}.
}
In case of a panel structure, one intercept is specified for each unit of observation, while the lag distributions are the same for all units:
\deqn{y_{i,t} = \alpha_i+\sum_{j=1}^J\sum_{k=0}^\infty\beta_{j,k}(\theta_j,\delta_j,\lambda_j,\eta_j)\hspace{0.07cm}x_{i,j,t-k}+\varepsilon_{i,t}}
where:
\itemize{
\item{\eqn{y_{i,t}} and \eqn{x_{i,j,t}} are, respectively, the value of \eqn{Y} and of \eqn{X_j} (\eqn{j=1,\ldots,J}) observed on unit \eqn{i} at time \eqn{t};}
\item{\eqn{\alpha_i} is the intercept for unit \eqn{i};}
\item{\eqn{\varepsilon_{i,t}} is the random error for unit \eqn{i} at time \eqn{t}.}
}
The main functions of the package are:
\itemize{
  \item{\link{adfTest} and \link{tsDiff}, }{to check weak stationarity and to apply differencing in order to achieve it;}
  \item{\link{gammadlm}, }{to estimate the model through the hill climbing algorithm;}
  \item{\link{lagCoef}, }{to see the estimated dynamic coefficients;}
  \item{\link{plot.gammadlm}, }{to display the estimated lag distributions.}
}
%method \code{residuals} to obtain graphical model diagnostics.
}
\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
A. Magrini (2022). A hill climbing algorithm for maximum likelihood estimation of the Gamma distributed-lag model with multiple explanatory variables. \emph{Austrian Journal of Statistics}, 51(2): 40-46.
}
