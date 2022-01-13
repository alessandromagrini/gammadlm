# gammadlm
__Maximum likelihood estimation of the Gamma distributed-lag model with multiple explanatory variables__

`gammadlm` is an R package implementing the hill climbing algorithm for maximum likelihood estimation of the Gamma distributed-lag model with multiple explanatory variables. A panel structure can be taken into account. The reference paper is:

A. Magrini (2022). A hill climbing algorithm for maximum likelihood estimation of the Gamma distributed-lag model with multiple explanatory variables. _Austrian Journal of Statistics_, 51(2): 40-46. Online first: https://ajs.or.at/index.php/ajs/article/view/1244


R (The R Project for Statistical Computing) needs to be installed on your system in order
to use the `gammadlm` package. R can be downloaded from https://www.r-project.org/.

To install the `gammadlm` package, open the console of R and type:
```
install.packages("devtools")  ## do not run if package 'devtools' is already installed
library(devtools)
install_github("alessandromagrini/gammadlm")
```

For any request or feedback, please write to <alessandro.magrini@unifi.it> (Alessandro Magrini)

Below, you find some examples of use of the package.
_________________________________________________________________

Load the `gammadlm` package
```
library(gammadlm)
```
Load data on Bitcoin and US stock market price
```
data(USstock)
mydata <- USstock[which(USstock$Date>="2020-04-01"),]  ## select data from April 2020 onwards
```
Trasform the variables in log return to achieve weak stationarity
```
mydataLR <- preProcess(time="Date", data=mydata, box.cox=0, ndiff=1)
```
Estimation with fixed values of delta and lambda parameters (simple OLS):
```
dval <- list(DJA=0.85,IXIC=0.75,GSPC=0.55)
lval <- list(DJA=0.05, IXIC=0.35,GSPC=0.45)
m0 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR,
  control=list(delta.lim=dval, lambda.lim=lval))
summary(m0)  ## summary of estimation
```
Estimation with unknown delta and lambda parameters (hill climbing algorithm):
```
# no random restarts (default: initialization is based on independent OLS-grid searches)
m1 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR)
summary(m1)

# add constraints: peak>=1 and 3<=length<=15
pklim <- list(DJA=c(1,Inf),IXIC=c(1,Inf),GSPC=c(1,Inf))
lenlim <- list(DJA=c(3,15),IXIC=c(3,15),GSPC=c(3,15))
m2 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR,
  control=list(peak.lim=pklim, length.lim=lenlim))
summary(m2)

# 50 random restarts without constraints
set.seed(100)
m3 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR,
  control=list(nstart=50))
summary(m3)

# 50 random restarts with constraints
set.seed(100)
m4 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR,
  control=list(peak.lim=pklim, length.lim=lenlim, nstart=50))
summary(m4)  ## same as m2
```
Estimated dynamic coefficients:
```
lagCoef(m4)
lagCoef(m4, cumul=T)  ## cumulative coefficients
```
Graphic of the estimated lag distributions:
```
plot(m4)
```
