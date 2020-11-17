\name{gammadlm-package}
\alias{gammadlm-package}
\docType{package}
\title{The Gamma Distributed-Lag Model with Multiple Explanatory Variables}
\description{Maximum likelihood estimation and inference for the Gamma distributed-lag model with multiple explanatory variables.}
\details{
  \tabular{ll}{
  Package: \tab gammadlm\cr
  Type: \tab Package\cr
  Version: \tab 0.0\cr
  Date: \tab 2020-11-17\cr
  License: \tab GPL-2
  }
Let \eqn{Y} be the response variable and \eqn{X_1,\ldots,X_J} be \eqn{J} explanatory variables.
Also, let \eqn{y_t} and \eqn{x_{j,t}} be, respectively, the value of \eqn{Y} and of \eqn{X_j} (\eqn{j=1,\ldots,J}) at time \eqn{t}.
Under the assumption that the time series of \eqn{Y} and of \eqn{X_1,\ldots,X_J} are all stationary, the Gamma distributed-lag model explaining \eqn{Y} from \eqn{X_1,\ldots,X_J} is defined as:
\deqn{y_t = \alpha+\sum_{j=1}^J\theta_j\sum_{k=0}^\infty w_{j,k}x_{j,t-k}+\varepsilon_t}
\deqn{w_{j,k} = \frac{(k+1-\eta_j)^{\frac{\delta_j}{1-\delta_j}}\lambda_j^{k-\eta_j}}{\sum_{l=0}^\infty(l+1-\eta_j)^{\frac{\delta_j}{1-\delta_j}}\lambda_j^{l-\eta_j}}}
where:
\itemize{
\item{\eqn{\alpha} is the intercept};
\item{\eqn{\theta_j} is the scale parameter of the lag distribution of \eqn{X_j}}; \item{\eqn{0\leq\delta_j<1} and \eqn{0\leq\lambda_j<1} are the shape parameters of the lag distribution of \eqn{X_j};}
\item{\eqn{\eta_j} is the offset of the lag distribution of \eqn{X_j}, typically set to 0;}
\item{\eqn{\varepsilon_t} is the random error at time \eqn{t}}.
}
The main functions of the package are:
\itemize{
  \item{\link{gammadlm}, }{to estimate the model through the hill climbing algorithm;}
  \item{\link{lagCoef}, }{to obtain the estimates of dynamic coefficients.}
  }
Also, method \link{lagCoef} can be used to obtain dynamic coefficients, method \code{plot} to obtain graphics for the estimated lag distributions, and method \code{residuals} to obtain graphical model diagnostics.
}
\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
A. Magrini. A hill climbing algorithm for maximum likelihood estimation of the Gamma distributed-lag model with multiple explanatory variables. \emph{Under review}.

J. M. Alston, M. A. Andersen, J. S. James, P. G. Pardey (2011). The economic returns to U.S. public agricultural research. \emph{American Journal of Agricultural Economics}, 93(5), 1257-1277.
}
