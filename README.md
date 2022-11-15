# gammadlm
__Maximum likelihood estimation of the Gamma distributed-lag model with multiple explanatory variables__

`gammadlm` is an R package implementing the hill climbing algorithm for maximum likelihood estimation of the Gamma distributed-lag model with multiple explanatory variables. Data may have a panel structure, even unbalanced. The reference paper is:

A. Magrini (2022). A hill climbing algorithm for maximum likelihood estimation of the Gamma distributed-lag model with multiple explanatory variables. _Austrian Journal of Statistics_, 51(2): 40-46. DOI: <a href="https://doi.org/10.17713/ajs.v51i2.1244">10.17713/ajs.v51i2.1244</a>


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
data(btc)
mydata <- btc[which(btc$Date>="2020-04-01"),]  ## select data from April 2020 onwards
```
Model with pre-specified shape parameters (simple OLS) on data in logarithmic differences (yearly relative changes):
```
dval <- list(DJA=0.85,IXIC=0.75,GSPC=0.55)
lval <- list(DJA=0.05, IXIC=0.35,GSPC=0.45)
m0 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), time="Date", data=mydata,
  box.cox=0, ndiff=1, control=list(delta.lim=dval, lambda.lim=lval))
summary(m0)  ## summary of estimation
```
Model with unknown delta and lambda parameters (hill climbing algorithm):
```
# no random restarts (default: initialization is based on independent OLS-grid searches)
m1 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydata, box.cox=0, ndiff=1)
summary(m1)

# add constraints: peak>=1 and 3<=length<=15
pklim <- list(DJA=c(1,Inf),IXIC=c(1,Inf),GSPC=c(1,Inf))
lenlim <- list(DJA=c(3,15),IXIC=c(3,15),GSPC=c(3,15))
m2 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydata,
  box.cox=0, ndiff=1, control=list(peak.lim=pklim, length.lim=lenlim))
summary(m2)

# 50 random restarts without constraints
set.seed(100)
m3 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydata,
  box.cox=0, ndiff=1, control=list(nstart=50))
summary(m3)  ## better fit than m1

# 50 random restarts with constraints
set.seed(100)
m4 <- gammadlm(y.name="BTC", x.names=c("DJA","IXIC","GSPC"), data=mydata,
  box.cox=0, ndiff=1, control=list(peak.lim=pklim, length.lim=lenlim, nstart=50))
summary(m4)  ## better fit than m2
```
Estimated dynamic multipliers:
```
lagCoef(m4)
lagCoef(m4, cumul=T)  ## cumulative multipliers
```
Graphic of the estimated lag distributions:
```
plot(m4)
```
Example with panel data:
```
# load data on farm performance and subsidies in EU countries (FADN database)
data(fadn)

# considered variables:
#  Y: productivity
#  X1-X4: subsidies by typology
#  Z1: utilized agricultural area (dimensional class),
#  Z2: total output (economic class)
y_name <- "TFP"
x_name <- c("Subs_prod","Subs_inv","Subs_rur","Subs_dec")
z_name <- c("Land","Output_total")

# model on data in differences: in this case, the intercept represents the coefficient of a linear trend (drift)
m_fadn <- gammadlm(y.name=y_name, x.names=x_name, z.names=z_name,
  unit="Country", time="Year", data=fadn, box.cox=1, ndiff=1,
  control=list(peak.lim=c(1,Inf), length.lim=c(3,10)))
summary(m_fadn)
