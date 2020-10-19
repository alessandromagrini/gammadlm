# gammadlm
__The Gamma distributed-lag model with multiple explanatory variables__

`gammadlm` is a R package implementing functionalities for maximum likelihood estimation and inference for the Gamma distributed-lag model with multiple explanatory variables.
The reference paper is:

A. Magrini. A hill climbing algorithm for maximum likelihood estimation of the Gamma distributed-lag model with multiple explanatory variables. _Under review_.


R (The R Project for Statistical Computing) needs to be installed on your system in order
to use the `gammadlm` package. R can be downloaded from https://www.r-project.org/.

To install the `CobbDouglas` package, open the console of R and type:
```
install.packages("devtools")  ## do not run if package 'devtools' is already installed
library(devtools)
install_github("alessandromagrini/gammadlm")
```

For any request or feedback, please write to <alessandro.magrini@unifi.it> (Alessandro Magrini)

Below, you find some examples of use of the package.
_________________________________________________________________

Load Bitcoin and US stock market data
```
data(USstock)
```
Trasform the variables in log return to achieve stationarity
```
dataLR <- tsDiff(time.name="Date", data=USstock, ndiff=1, log=TRUE)
```
Estimation with fixed values of delta and lambda parameters: 1-step OLS
```
dval <- list(DJA=0.85,IXIC=0.75,GSPC=0.55)
lval <- list(DJA=0.5, IXIC=0.35,GSPC=0.45)
mod <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=dataLR, control=list(delta.lim=dval, lambda.lim=lval))
summary(mod)  ## summary of estimation
```
Estimation through hill climbing algorithm
```
# * no constraints
set.seed(100)
m1 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=dataLR)
summary(m1)

# * same but with 100 random restarts
set.seed(100)
m1a <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=dataLR, control=list(nstart=100))
summary(m1a)

# * constraints: peak>=1 and 3<=length<=10
pklim <- list(DJA=c(1,Inf),IXIC=c(1,Inf),GSPC=c(1,Inf))
lenlim <- list(DJA=c(3,10),IXIC=c(3,10),GSPC=c(3,10))
set.seed(100)
m2 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=dataLR, control=list(peak.lim=pklim, length.lim=lenlim, nstart=100))
summary(m2)
```
Graphical model diagnostics
```
residuals(mod, plot=T)
```
Estimated dynamic coefficients
```
lagCoef(mod)
lagCoef(mod, cumul=T)  ## cumulative coefficients
```
Graphic of the estimated lag distributions
```
plot(mod)
```
