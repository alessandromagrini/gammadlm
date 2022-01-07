# gammadlm
__The Gamma distributed-lag model with multiple explanatory variables__

`gammadlm` is an R package implementing maximum likelihood estimation and inference for the Gamma distributed-lag model with multiple explanatory variables.
The reference paper is:

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
Estimation with fixed values of delta and lambda parameters (1 step of ordinary least squares):
```
dval <- list(DJA=0.85,IXIC=0.75,GSPC=0.55)
lval <- list(DJA=0.5, IXIC=0.35,GSPC=0.45)
mod <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR,
                control=list(delta.lim=dval, lambda.lim=lval))
summary(mod)  ## summary of estimation
```
Estimation with unknown delta and lambda parameters (hill climbing algorithm):
```
# * no constraints with 50 random restarts (default)
set.seed(100)
m1 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR)
summary(m1)

# * no constraints with 100 random restarts
set.seed(100)
m1a <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR,
                control=list(nstart=100))
summary(m1a)

# * constraints: peak>=1 and 3<=length<=10
pklim <- list(DJA=c(1,Inf),IXIC=c(1,Inf),GSPC=c(1,Inf))
lenlim <- list(DJA=c(3,10),IXIC=c(3,10),GSPC=c(3,10))
set.seed(100)
m2 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydataLR,
               control=list(peak.lim=pklim, length.lim=lenlim, nstart=100))
summary(m2)
```
Estimated dynamic coefficients:
```
lagCoef(mod)
lagCoef(mod, cumul=T)  ## cumulative coefficients
```
Graphic of the estimated lag distributions:
```
plot(mod)
```
