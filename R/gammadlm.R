### DA FARE
#
# - differencing automatico
# - diagnostiche grafiche
# - metodo predict (inversione differenze e box.cox)
# - quandt test
# - draw sample
# - gam_inits(): migliorare il campionamento
# - constraint segno
# - full covariance matrix


# apply box-cox transformation (auxiliary)
makeBoxCox <- function(x, par) {
  if(par==1) {
    x
    } else if(par==0) {
    #ind <- which(x<=0)
    #if(length(ind)>0) {
    #  x[ind] <- min(x[setdiff(1:length(x),ind)])/2
    #  }
    log(x)
    } else {
    (x^par-1)/par  
    }
  }

# invert box-cox transformation (auxiliary)
invertBoxCox <- function(z, par) {
  if(par==1) {
    z
    } else if(par==0) {
    #ind <- which(x<=0)
    #if(length(ind)>0) {
    #  x[ind] <- min(x[setdiff(1:length(x),ind)])/2
    #  }
    exp(z)
    } else {
    #z=(x^par-1)/par  
    (z*par+1)^(1/par)
    }
  }

# linear interpolation (auxiliary)
linInterp <- function(x) {
  if(sum(!is.na(x))>=2) {
    approx(x, xout=1:length(x))$y
    } else {
    x
    }
  }

# unit root test for one variabile (auxiliary)
oneTest <- function(x, unit=NULL, max.lag=NULL) {
  if(is.null(unit)) {
    x <- na.omit(linInterp(x))
    n <- length(x)
    #if(n<5) stop("At least 5 observations are required",call.=F)
    } else {
    gr <- levels(factor(unit))
    nvet <- c()
    for(w in gr) {
      ind <- which(unit==w)
      x[ind] <- linInterp(x[ind])
      nvet[w] <- length(na.omit(x[ind]))
      #if(nvet[w]<5) stop("At least 5 observations are required for each unit",call.=F)
      }
    isOK <- which(!is.na(x))
    x <- x[isOK]
    unit <- unit[isOK]
    n <- min(nvet)
    }
  max.lag <- max.lag[1]
  if(!is.null(max.lag)&&is.na(max.lag)) max.lag <- NULL
  if(!is.null(max.lag)&&!is.numeric(max.lag)) max.lag <- NULL
  if(is.null(max.lag)) {
    max.lag <- round(sqrt(n))
    } else {
    max.lag <- round(min(max(0,max.lag),sqrt(n)))
    }
  if(is.null(unit)) {
    res <- res1 <- adfFun(x=x, max.lag=max.lag)
    res2 <- kpssFun(x=x, max.lag=max.lag)
    for(i in 1:3) {
      res[[i]] <- c(adf=res1[[i]],kpss=res2[[i]])
      }
    } else {
    gr <- levels(factor(unit))
    res1 <- res2 <- vector("list",length=3)
    for(w in gr) {
      ind <- which(unit==w)
      iadf <- adfFun(x=x[ind], max.lag=max.lag)
      ikpss <- kpssFun(x=x[ind], max.lag=max.lag)
      for(j in 1:length(res1)) {
        res1[[j]] <- c(res1[[j]],iadf[[j]])
        res2[[j]] <- c(res2[[j]],ikpss[[j]])
        }
      }
    #
    pvalComb <- function(x) {
      x[which(x<=0)] <- 1e-8
      x[which(x>=1)] <- 1-1e-8
      ind <- which(!is.na(x))
      if(length(ind)>0) {
        m <- length(ind)
        logp <- qnorm(x[ind])
        rhat <- 1-var(logp)
        rstar <- max(rhat,-1/(m-1))
        auxz <- sum(logp)/sqrt(m*(1+(m-1)*(rstar+0.2*sqrt(2/(m+1))*(1-rstar))))
        #auxz <- sum(logp)/sqrt(m)
        c(x,'(combined)'=2*pnorm(-abs(auxz)))
        } else {
        c(x,'(combined)'=NaN)
        }
      }
    #
    res1 <- lapply(res1, function(z){names(z)<-gr; z})
    names(res1) <- names(iadf)
    res1$p.value <- pvalComb(res1$p.value)
    res2 <- lapply(res2, function(z){names(z)<-gr; z})
    names(res2) <- names(ikpss)
    res2$p.value <- pvalComb(res2$p.value)
    res <- res1
    for(i in 1:3) {
      res[[i]] <- cbind(adf=res1[[i]],kpss=res2[[i]])
      }
    }
  res
  }

# function for adf test (auxiliary)  
adfFun <- function(x, max.lag) {
  #
  doADF <- function(k) {
    y <- diff(x)
    n <- length(y)
    k <- k+1
    z <- embed(y,k)
    yt <- z[,1]
    xt1 <- x[k:n]
    tt <- k:n
    if(k>1) {
      yt1 <- z[,2:k,drop=F]
      res <- lm(yt~xt1+tt+yt1)
      } else {
      res <- lm(yt~xt1+tt)
      }
    suppressWarnings(res.sum <- summary.lm(res)$coefficients)
    if(nrow(res.sum)>=2) {
      STAT <- res.sum[2,1]/res.sum[2,2]
      table <- -1*cbind(c(4.38, 4.15, 4.04, 3.99, 3.98, 3.96),
                        c(3.95, 3.8, 3.73, 3.69, 3.68, 3.66),
                        c(3.6, 3.5, 3.45, 3.43, 3.42, 3.41),
                        c(3.24, 3.18, 3.15, 3.13, 3.13, 3.12),
                        c(1.14, 1.19, 1.22, 1.23, 1.24, 1.25),
                        c(0.8, 0.87, 0.9, 0.92, 0.93, 0.94),
                        c(0.5, 0.58, 0.62, 0.64, 0.65, 0.66),
                        c(0.15, 0.24, 0.28, 0.31, 0.32, 0.33))
      tablen <- dim(table)[2]
      tableT <- c(25, 50, 100, 250, 500, 1e+05)
      tablep <- c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)
      tableipl <- numeric(tablen)
      for(i in (1:tablen)) {
        tableipl[i] <- approx(tableT,table[,i],n,rule=2)$y
        }
      PVAL <- approx(tableipl,tablep,STAT,rule=2)$y
      } else {
      STAT <- PVAL <- NA
      }
    c(STAT,PVAL)
    }
  k <- 0
  if(length(x)>=5 & var(x)>0) {  ## <--
    #if(max.lag>0) k <- ar(x, order.max=max.lag)$order
    if(max.lag>0) k <- lagSelect(x, order.max=max.lag)
    res <- doADF(k)
    } else {
    res <- c(NaN,NaN)
    }
  list(statistic=res[1], lag.order=k, p.value=res[2])
  }

# function for kpss test (internal use only)
kpssFun <- function(x, max.lag) {
  #
  doKPSS <- function(lag) {
    n <- length(x)
    #if(trend==T) {
    t <- 1:n
    e <- residuals.lm(lm(x ~ t))
    table <- c(0.216, 0.176, 0.146, 0.119)
    #  } else {
    #  e <- residuals.lm(lm(x ~ 1))
    #  table <- c(0.739, 0.574, 0.463, 0.347)
    #  }
    tablep <- c(0.01, 0.025, 0.05, 0.1)
    s <- cumsum(e)
    eta <- sum(s^2)/(n^2)
    s2 <- sum(e^2)/n
    k <- 0
    for(i in 1:lag) {
      ik <- 0
      for(j in (i+1):n) {
        ik <- ik+e[j]*e[j-i]
        }
      k <- k+(1-i/(lag+1))*ik
      }
    STAT <- eta/(s2+2*k/n)
    PVAL <- approx(table,tablep,STAT,rule=2)$y
    c(statistic=STAT, p.value=PVAL)
    }
  #
  k <- 0
  if(length(x)>=5 & var(x)>0) {  ## <--
    #if(max.lag>0) k <- ar(x, order.max=max.lag)$order
    if(max.lag>0) k <- lagSelect(x, order.max=max.lag)
    res <- doKPSS(k)
    } else {
    res <- c(NaN,NaN)
    }
  list(statistic=unname(res[1]), lag.order=k, p.value=unname(res[2]))
  }

# recognize quantitative variable (auxiliary)
isQuant <- function(x) {
  is.numeric(x)&!setequal(unique(na.omit(x)),c(0,1))
  }

# perform unit root test
unirootTest <- function(var.names, unit=NULL, time=NULL, data, box.cox=1, ndiff=0, max.lag=NULL) {
  dataD <- preProcess(var.names=var.names, unit=unit, time=time, data=data, box.cox=box.cox, ndiff=ndiff)
  xcat <- c()
  for(i in 1:length(var.names)) {
    if(isQuant(data[,var.names[i]])==F) {
      xcat <- c(xcat,var.names[i])
      warning("Variable '",var.names[i],"' is not quantitative and has been ignored",call.=F)
      }
    }
  var.names <- setdiff(var.names,xcat)
  if(length(var.names)==0) stop("No quantitative variable provided to argument 'var.names'")
  #
  if(is.null(unit)) gr <- NULL else gr <- dataD[,unit]
  max.lag <- max.lag[1]
  if(!is.numeric(max.lag)) max.lag <- NULL else max.lag <- max(0,ceiling(max.lag))
  testList <- list()
  for(i in 1:length(var.names)) {
    iadf <- oneTest(x=dataD[,var.names[i]], unit=gr, max.lag=max.lag)
    iadf$box.cox <- unname(attr(dataD,"box.cox")[var.names[i]])
    iadf$ndiff <- unname(attr(dataD,"ndiff")[var.names[i]])
    testList[[i]] <- iadf
    }
  names(testList) <- var.names
  class(testList) <- "unirootTest"
  testList
  }

# print method for class 'unirootTest'
print.unirootTest <- function(x, ...) {
  cat("p-values","\n")
  cat("  ADF:  null hypothesis is 'unit root'","\n")
  cat("  KPSS: null hypothesis is 'no unit roots'","\n")
  tab <- sapply(x, function(z){
    pval <- round(z$p.value,4)
    if(is.matrix(pval)) pval["(combined)",] else pval
    })
  print(t(tab))
  }

# function for differencing (auxiliary)
diffFun <- function(var.names, data, ndiff) {
  newdat <- data
  n <- nrow(data)
  for(i in 1:length(var.names)) {
    if(ndiff[var.names[i]]>0) {
      idat <- linInterp(data[,var.names[i]])
      idat_lag <- c(rep(NA,ndiff[i]),idat)[1:length(idat)]
      newdat[,var.names[i]] <- idat-idat_lag
      }
    }
  attr(newdat,"ndiff") <- ndiff
  newdat
  }

# pre-processing (auxiliary)
preProcess <- function(var.names, unit=NULL, time=NULL, data, box.cox=1, ndiff=0) {
  #
  if(missing(data)) stop("Argument 'data' is missing")
  if(!identical(class(data),"data.frame")) stop("Argument 'data' must be a data.frame")
  if(missing(var.names)) stop("Argument 'var.names' is missing")
  if(!is.character(var.names)) {
    stop("Argument 'var.names' must be a character vector")
    } else {
    var.names[which(!is.na(var.names))]
    if(length(var.names)<1) stop("Argument 'var.names' must be a character vector of length 1 or greater")
    }
  auxchk <- setdiff(var.names,colnames(data))  
  if(length(auxchk)>0) stop("Variable '",auxchk[1],"' not found")
  xnum <- xcat <- c()
  for(i in 1:length(var.names)) {
    if(isQuant(data[,var.names[i]])) {
      xnum <- c(xnum,var.names[i])
      } else {
      if(sum(is.na(data[,var.names[i]]))>0) stop("Variable '",var.names[i],"' is categorical and contains missing values")
      xcat <- c(xcat,var.names[i])
      }
    }
  #
  unit <- unit[1]
  if(!is.null(unit)&&is.na(unit)) unit <- NULL
  if(!is.null(unit)) {
    if(length(setdiff(unit,colnames(data)))>0) stop("Variable '",unit,"' not found")
    if(length(intersect(unit,var.names))>0) stop("Variable '",unit,"' appears in both arguments 'var.names' and 'unit'")
    if(sum(is.na(data[,unit]))>0) stop("Variable '",unit,"' provided to argument 'unit' contains missing values")
    data[,unit] <- factor(data[,unit])
    }
  #
  time <- time[1]
  if(!is.null(time)&&is.na(time)) time <- NULL
  if(!is.null(time)) {
    if(length(setdiff(time,colnames(data)))>0) stop("Variable '",time,"' not found")
    if(length(intersect(time,var.names))>0) stop("Variable '",time,"' appears in both arguments 'var.names' and 'time'")
    if(length(intersect(time,unit))>0) stop("Variable '",time,"' appears in both arguments 'unit' and 'time'")
    if(!is.numeric(data[,time])&!identical(class(data[,time]),"Date")) stop("Variable '",time,"' must be numeric or of class 'Date'")
    if(sum(is.na(data[,time]))>0) stop("Variable '",time,"' provided to argument 'time' contains missing values")
    if(is.null(unit)) {
      if(sum(duplicated(data[,time]))>0) stop("Variable '",time,"' contains duplicated values")
      } else {
      if(sum(sapply(split(data,data[,unit]),function(x){sum(duplicated(x[,time]))}))>0) stop("Variable '",time,"' contains duplicated values")
      }
    }
  #
  if(!is.numeric(box.cox)) stop("Argument 'box.cox' must be a numeric value or vector")
  if(length(box.cox)==1) {
    box.cox <- rep(box.cox,length(xnum))
    names(box.cox) <- xnum
    } else if(is.null(names(box.cox))) {
    box.cox <- box.cox[1:length(xnum)]
    names(box.cox) <- xnum
    box.cox[which(is.na(box.cox))] <- 1
    } else {
    box.cox <- box.cox[xnum]
    names(box.cox) <- xnum
    box.cox[which(is.na(box.cox))] <- 1
    }
  if(sum(box.cox<0)>0) stop("Argument 'box.cox' must contain non-negative real values")
  #
  if(!is.numeric(ndiff)) stop("Argument 'ndiff' must be a numeric value or vector")
  if(length(ndiff)==1) {
    ndiff <- rep(ndiff,length(xnum))
    names(ndiff) <- xnum
    } else if(is.null(names(ndiff))) {
    ndiff <- ndiff[1:length(xnum)]
    names(ndiff) <- xnum
    ndiff[which(is.na(ndiff))] <- 0
    } else {
    ndiff <- ndiff[xnum]
    names(ndiff) <- xnum
    ndiff[which(is.na(ndiff))] <- 0
    }
  if(sum(ndiff<0|round(ndiff)!=ndiff)>0) stop("Argument 'ndiff' must contain non-negative integer values")
  # sort by time
  if(!is.null(time)) {
    if(is.null(unit)) {
      data <- data[order(data[,time]),]
      } else {
      data <- data[order(data[,unit],data[,time]),]
      }
    }
  # box-cox
  dataL <- data
  if(length(xnum)>0) {
    for(i in 1:length(xnum)) {
      ilam <- box.cox[[i]]
      if(ilam!=1 & sum(data[,xnum[i]]<0,na.rm=T)>0) {
        box.cox[[i]] <- ilam <- 1
        warning("No transformation applied to variable '",xnum[i],"'",call.=F)
        } else if(ilam==0 & sum(data[,xnum[i]]==0,na.rm=T)>0) {
        box.cox[[i]] <- ilam <- 0.5
        warning("Box-Cox parameter 0.5, instead of 0, applied to variable '",xnum[i],"'",call.=F)
        }
      dataL[,xnum[i]] <- makeBoxCox(data[,xnum[i]],ilam)
      }
    }
  # differencing
  if(length(xnum)>0) {
    if(is.null(unit)) {
      n <- nrow(data)
      for(i in 1:length(xnum)) {
        if(ndiff[xnum[i]]>n-5) {
          ndiff[xnum[i]] <- 0
          warning("Differencing not applied to variable '",xnum[i],"'",call.=F)
          }
        }
      dataD <- diffFun(var.names=xnum, data=dataL, ndiff=ndiff)
      attr(dataD,"box.cox") <- box.cox
      attr(dataD,"ndiff") <- ndiff
      if(max(ndiff)>0) {
        res <- dataD[setdiff(1:nrow(data),1:max(ndiff)),,drop=F]
        } else {
        res <- dataD
        }
      } else {
      data[,unit] <- factor(data[,unit])
      dataD <- dataL
      isNA <- c()
      gr <- levels(data[,unit])
      val0 <- matrix(nrow=length(gr),ncol=length(var.names))
      rownames(val0) <- gr
      colnames(val0) <- var.names
      n_gr <- c()
      for(w in 1:length(gr)) {
        ind <- which(data[,unit]==gr[w])
        n_gr[w] <- length(ind)
        if(max(ndiff)>0) isNA <- c(isNA, ind[1]:ind[max(ndiff)])
        dataD[ind,] <- diffFun(var.names=xnum, data=dataL[ind,], ndiff=ndiff)
        val0[w,xnum] <- as.numeric(data[ind[1],xnum])
        }
      n <- max(n_gr)
      for(i in 1:length(xnum)) {
        if(ndiff[xnum[i]]>n-5) {
          ndiff[xnum[i]] <- 0
          warning("Differencing not applied to variable '",xnum[i],"'",call.=F)
          }
        }
      attr(dataD,"box.cox") <- box.cox
      attr(dataD,"ndiff") <- ndiff
      res <- dataD[setdiff(1:nrow(data),isNA),,drop=F]
      }
    res
    } else {
    dataL
    }
  }

# unnormalized gamma weights (auxiliary)
gam_wei <- function(k, par, offset=0) {
  delta <- par[1]
  lambda <- par[2]
  wf <- function(x) {
    if(delta>0) {
      if(x>offset-1) (x-offset+1)^(delta/(1-delta))*lambda^(x-offset) else 0
      } else {
      if(x>=offset) lambda^(x-offset) else 0
      }
    #if(x>offset-1) (x-offset+1)^(delta/(1-delta))*lambda^(x-offset) else 0
    }
  sapply(k,wf)
  }

# gamma weights
gammaWeights <- function(k, par, offset=0, normalize=TRUE) {
  #
  if(missing(k)) stop("Argument 'k' is missing")
  if(!is.numeric(k)) stop("Argument 'k' must be a numerical vector")
  #
  if(missing(par)) stop("Argument 'par' is missing")
  if(length(par)>=2) par <- par[1:2] else stop("Argument 'par' must be of length 2")
  if(!is.numeric(par)) stop("Argument 'par' must be a numeric vector of length 2")
  par[1] <- min(max(0,par[1]),0.985)
  par[2] <- min(max(0,par[2]),0.985)
  #if(!is.numeric(par) || (par[1]<0 | par[1]>=1 | par[2]<0 | par[2]>=1)) stop("Components of argument 'par' must be values >=0 and <1")
  #
  offset <- offset[1]
  if(!is.numeric(offset)) offset <- 0
  #
  normalize <- normalize[1]
  if(is.na(normalize)||(!is.logical(normalize)|is.null(normalize))) normalize <- TRUE
  #
  auxwei <- gam_wei(k=k, par=par, offset=offset)
  if(normalize) {
    swei <- gam_const(par=par)
    } else {
    swei <- 1
    }
  auxwei/swei
  }

# single gamma kernel projection (auxiliary)
gamkern <- function(x, par, offset, normalize) {
  n <- length(x)
  wei <- gammaWeights(0:(n-1), par=par, offset=offset, normalize=normalize)
  res <- c()
  for(i in 1:n) {
    res[i] <- sum(wei[1:i]*x[i-(0:(i-1))])
    }
  res
  #x2 <- c(rep(mean(x,na.rm=T),length(x)),x)
  #n <- length(x2)
  #wei <- gammaWeights(0:(n-1), par=par, offset=offset, normalize=normalize)
  #res <- c()
  #for(i in 1:n) {
  #  res[i] <- sum(wei[1:i]*x2[i-(0:(i-1))])
  #  }
  #res[(length(x)+1):n]
  }

# gamma kernel projection
gammaKernel <- function(x, par, unit=NULL, offset=0, normalize=TRUE) {
  if(missing(x)) stop("Argument 'x' is missing")
  if(!is.numeric(x)) stop("Argument 'x' must be a numerical vector")
  #
  if(missing(par)) stop("Argument 'par' is missing")
  if(length(par)>=2) par <- par[1:2] else stop("Argument 'par' must be of length 2")
  if(!is.numeric(par)) stop("Argument 'par' must be a numeric vector of length 2")
  par[1] <- min(max(0,par[1]),0.985)
  par[2] <- min(max(0,par[2]),0.985)
  #if(!is.numeric(par) || (par[1]<0 | par[1]>=1 | par[2]<0 | par[2]>=1)) stop("Components of argument 'par' must be values >=0 and <1")
  #
  offset <- offset[1]
  if(!is.numeric(offset)) offset <- 0
  #
  normalize <- normalize[1]
  if(is.na(normalize)||(!is.logical(normalize)|is.null(normalize))) normalize <- TRUE
  #
  if(is.null(unit)) {
    gamkern(x=x, par=par, offset=offset, normalize=normalize)
    } else {
    gr <- levels(factor(unit))
    res <- rep(NA,length(x))
    for(w in gr) {
      ind <- which(unit==w)
      res[ind] <- gamkern(x=x[ind], par=par, offset=offset, normalize=normalize)
      }
    res
    }
  }

# gamma lag normalization constant (auxiliary)
gam_const <- function(par, tol=1e-12) {
  lval <- ifelse(par[2]>0, ceiling(par[1]/(par[1]-1)/log(par[2])), 0)
  swei <- sum(gam_wei(0:lval, par))
  fine <- 0
  while(fine==0) {
    lval <- lval+1
    wei <- gam_wei(lval, par)
    swei <- swei+wei
    if(wei<tol) fine <- 1
    }
  swei
  #if(par[2]==0) 1 else integrate(gam_wei, lower=-Inf, upper=Inf, par=par, ...)$value
  }

# find quantile of a gamma lag distribution
gammaQuantile <- function(prob, par, offset=0) {
  #
  if(missing(prob)) stop("Argument 'prob' is missing")
  if(!is.numeric(prob)) stop("Argument 'prob' must be a numeric vector")
  prob <- sapply(prob, function(x){min(max(1e-12,x),1-1e-12)})
  #if(!is.numeric(prob) || sum(prob<=0)>0 | sum(prob>=1)>0) stop("Argument 'prob' must be a vector of values >0 and <1")
  #
  if(missing(par)) stop("Argument 'par' is missing")
  if(length(par)>=2) par <- par[1:2] else stop("Argument 'par' must be of length 2")
  if(!is.numeric(par)) stop("Argument 'par' must be a numeric vector of length 2")
  par[1] <- min(max(0,par[1]),0.985)
  par[2] <- min(max(0,par[2]),0.985)
  #if(!is.numeric(par) || (par[1]<0 | par[1]>=1 | par[2]<0 | par[2]>=1)) stop("Components of argument 'par' must be values >=0 and <1")
  #
  offset <- offset[1]
  if(!is.numeric(offset)) offset <- 0
  #
  swei <- gam_const(par=par)
  qcalc <- function(p) {
    lval <- -1
    wei_old <- 0
    wei <- gam_wei(lval, par=par, offset=0)/swei
    while(wei<p) {
      lval <- lval+1
      w0 <- gam_wei(lval, par=par, offset=0)/swei
      wei_old <- wei
      wei <- wei+w0
      }
    max(0, approx(c(wei_old,wei),c(lval-1,lval),xout=p)$y+offset)
    }
  res <- sapply(prob, qcalc)
  names(res) <- paste(round(prob*100,3),"%",sep="")
  res
  }

# fit ols with gamma lag (auxiliary)
gam_olsFit <- function(y.name, x.names, z.names, unit, par, offset, data, normalize, add.intercept, multi.intercept, weights) {
  p <- length(x.names)
  if(normalize) normstr <- "" else normstr <- paste(",normalize=",normalize,sep="")
  if(is.null(unit)) {
    if(add.intercept) {
      form0 <- paste(y.name," ~ ",sep="")
      } else {
      form0 <- paste(y.name," ~ -1+",sep="")
      }
    gstr <- ""
    } else {
    if(add.intercept) {
      if(multi.intercept) {
        form0 <- paste(y.name," ~ -1+",unit,"+",sep="")
        } else {
        form0 <- paste(y.name," ~ ",sep="")
        }
      } else {
      form0 <- paste(y.name," ~ -1+",sep="")
      }
    gstr <- paste(",unit=",unit,sep="")
    }
  for(i in 1:p) {
    if(i>1) form0 <- paste(form0,"+",sep="")
    if(offset[i]==0) offstr <- "" else offstr <- paste(",offset=",offset[i],sep="")
    form0 <- paste(form0,"gammaKernel(",x.names[i],",par=c(",paste(par[,i],collapse=","),")",gstr,offstr,normstr,")",sep="")
    }
  if(length(z.names)>0) form0 <- paste(form0,"+",paste(z.names,collapse="+"),sep="")
  formOK <- formula(form0)
  mod <- lm(formOK,data=data,weights=weights)
  mod$call$formula <- formOK
  #mod$call$data <- deparse(substitute(data))
  #
  colnames(par) <- names(offset) <- x.names
  rownames(par) <- c("delta","lambda")
  mod$par <- par
  mod$offset <- offset
  mod$add.intercept <- add.intercept
  #mod$multi.intercept <- multi.intercept
  class(mod) <- c("gammadlm","lm")
  mod
  }

# compute peak (auxiliary)
peakCalc <- function(x) {
  if(x[1]==0|x[2]==0) 0 else x[1]/(x[1]-1)/log(x[2])-1
  }

# create grid (auxiliary)
gam_parGrid <- function(delta.lim, lambda.lim, peak.lim, length.lim, grid.by) {
  p <- length(delta.lim)
  stopped <- 0
  gridList <- vector("list",length=p)
  for(i in 1:p) {
    delval <- seq(max(delta.lim[[i]][1], 0),
                  min(delta.lim[[i]][2], 1-grid.by, 0.985), by=grid.by)
    lamval <- seq(max(lambda.lim[[i]][1], 0),
                  min(lambda.lim[[i]][2], 1-grid.by, 0.985), by=grid.by)
    gridMat <- as.matrix(expand.grid(delval,lamval))
    auxdelta <- which(gridMat[,2]==0)
    if(length(auxdelta)>0) gridMat <- rbind(c(0,0),gridMat[setdiff(1:nrow(gridMat),auxdelta),])
    ind2del <- c()
    for(j in 1:nrow(gridMat)) {
      if(peak.lim[[i]][1]> -Inf | peak.lim[[i]][2]<Inf) {
        xpeak <- peakCalc(gridMat[j,])
        if(xpeak<peak.lim[[i]][1] | xpeak>peak.lim[[i]][2]) ind2del <- c(ind2del,j)
        }
      if(length.lim[[i]][1]>0 | length.lim[[i]][2]<Inf) {
        xlen <- gammaQuantile(0.999, gridMat[j,])
        if(xlen<length.lim[[i]][1] | xlen>length.lim[[i]][2]) ind2del <- c(ind2del,j)
        }
      }
    indOK <- setdiff(1:nrow(gridMat),ind2del)
    if(length(indOK)>0) {
      gridList[[i]] <- gridMat[indOK,,drop=F]
      } else {
      stopped <- 1
      break()
      }
    }
  if(stopped==0) gridList
  }

# generate random inits (auxiliary)
gam_inits <- function(gridList, visitList, maxtry) {
  p <- length(gridList)  ## <----- migliorare efficienza
  ntry <- fine <- 0
  iniOK <- c()
  while(fine==0) {
    ntry <- ntry+1
    ini0 <- c()
    for(i in 1:p) {
      igr <- gridList[[i]]
      isam <- sample(1:nrow(igr),1)
      ini0 <- cbind(ini0,igr[isam,])
      }
    check0 <- visitFun(ini0, visitList)
    if(check0==0 | ntry>=maxtry) {
      if(check0==0) iniOK <- ini0
      fine <- 1
      }
    }
  iniOK
  }

# generate combinations of gamma and lambda (auxiliary)
gam_comb <- function(par, gridmat, grid.by) {
  combF <- function(s1,s2) {
    newpar <- c(par[1]+s1*grid.by,par[2]+s2*grid.by)
    aux <- which(gridmat[,1]==newpar[1] & gridmat[,2]==newpar[2])
    if(length(aux)>0) newpar
    }
  cbind(combF(0,0),
        combF(0,1),
        combF(0,-1),
        combF(1,0),
        combF(1,1),
        combF(1,-1),
        combF(-1,0),
        combF(-1,1),
        combF(-1,-1)
        )
  }

# check whether a parameter value was visited (auxiliary)
visitFun <- function(par, par.list) {
  if(length(par.list)>0) {
    sum(sapply(par.list, function(x){identical(x, par)}))
    } else {
    0
    }
  }

# function for hill climbing (auxiliary)
gam_hcFun <- function(y.name, x.names, z.names, unit, data, offset, inits, visitList, gridList, sign, grid.by, add.intercept, multi.intercept, weights) {
  p <- length(x.names)
  n <- nrow(data)
  logn <- log(n)
  modOK <- NULL
  rssOK <- Inf
  parOK <- inits
  fine <- 0
  testL <- list()
  while(fine==0) {
    ibest_rss <- Inf
    ibest_par <- inits
    for(i in 1:p) {
      icomb <- gam_comb(parOK[,i], gridmat=gridList[[i]], grid.by=grid.by)
      for(j in 1:ncol(icomb)) {
        ijpar <- parOK
        ijpar[,i] <- icomb[,j]
        check0 <- visitFun(ijpar, visitList)
        if(check0==0) {
          visitList <- c(visitList,list(ijpar))
          ijm <- gam_olsFit(y.name=y.name, x.names=x.names, z.names=z.names, unit=unit, par=ijpar, offset=offset, data=data, normalize=F, add.intercept=add.intercept, multi.intercept=multi.intercept, weights=weights)
          ijrss <- sum(ijm$residuals^2)
          testL <- c(testL,list(ijpar))
          #ijsign <- sign(ijm$coef[-1])   <------ gestire constraint segno
          #if(sum(abs(ijsign-sign)>=2)==0) {
          #  ijrss <- sum(ijm$residuals^2)
          #  } else {
          #  ijrss <- Inf
          #  }
          } else {
          ijrss <- Inf
          }
        if(ijrss<ibest_rss) {
          ibest_rss <- ijrss
          ibest_par <- ijpar
          ibest_mod <- ijm
          }
        }
      }
    if(ibest_rss<rssOK) {
      rssOK <- ibest_rss
      parOK <- ibest_par
      modOK <- ibest_mod
      } else {
      fine <- 1
      }
    }
  modFinal <- gam_olsFit(y.name=y.name, x.names=x.names, z.names=z.names, unit=unit, par=parOK, offset=offset, data=data, normalize=T, add.intercept=add.intercept, multi.intercept=multi.intercept, weights=weights)
  list(model=modFinal,par.tested=testL)
  }

# format control options (auxiliary)
optFormat <- function(optList, nomi, val) {
  auxopt <- optList[nomi]
  if(is.null(auxopt)) auxopt <- vector("list",length=length(nomi))
  val <- na.omit(val)[1:2]
  if(is.na(val[1])) val[1] <- -Inf
  if(is.na(val[2])) val[2] <- Inf
  for(i in 1:length(nomi)) {
    iopt <- na.omit(auxopt[[i]])[1:2]
    if(!is.numeric(iopt)) {
      iopt <- val
      } else {
      if(is.na(iopt[1])) iopt[1] <- val[1]
      if(is.na(iopt[2])) iopt[2] <- val[2]
      iopt[1] <- min(max(iopt[1],val[1]),val[2])
      iopt[2] <- min(max(iopt[2],val[1]),val[2])
      if(iopt[1]>iopt[2]) iopt[2] <- iopt[1]
      }
    auxopt[[i]] <- iopt
    }
  names(auxopt) <- nomi
  auxopt
  }

# MASTER FUNCTION
gammadlm <- function(y.name, x.names, z.names=NULL, unit=NULL, time=NULL, data,
  offset=rep(0,length(x.names)), box.cox=1, ndiff=0, add.intercept=TRUE, weights=NULL,
  control=list(nstart=NULL, delta.lim=NULL, lambda.lim=NULL, peak.lim=NULL, length.lim=NULL), quiet=FALSE) {
  #
  if(!identical(class(data),"data.frame")) stop("Argument 'data' must be a data.frame")
  #
  if(missing(y.name)) stop("Argument 'y.name' is missing")
  if(!is.character(y.name)) {
    stop("Argument 'y.name' must be a character vector of length 1")
    } else {
    y.name <- y.name[which(!is.na(y.name))]
    if(length(y.name)!=1) stop("Argument 'y.name' must be a character vector of length 1")
    }
  auxchk <- setdiff(y.name,colnames(data))  
  if(length(auxchk)>0) stop("Variable '",auxchk[1],"' not found")
  #
  if(missing(x.names)) stop("Argument 'x.names' is missing")
  if(!is.character(x.names)) {
    stop("Argument 'x.names' must be a character vector")
    } else {
    x.names <- x.names[which(!is.na(x.names))]
    if(length(x.names)<1) stop("Argument 'x.names' must be a character vector")
    for(i in 1:length(x.names)) {
      if(isQuant(data[,x.names[i]])==F) stop("Variable '",x.names[i],"' is not quantitative")
      }
    }
  auxchk <- setdiff(x.names,colnames(data))
  if(length(auxchk)>0) stop("Variable '",auxchk[1],"' not found")
  #if(y.name%in%x.names) stop("Variable '",y.name,"' appears in both arguments 'y.name' and 'x.names'")
  #
  if(!is.null(z.names)) {
    z.names <- z.names[which(!is.na(z.names))]
    if(!is.character(z.names)|length(z.names)==0) stop("Argument 'z.names' must be either a character vector or NULL")
    if(length(z.names)>0) {
      auxchk2 <- setdiff(z.names,colnames(data))
      if(length(auxchk2)>0) stop("Variable '",auxchk2[1],"' not found")
      if(y.name%in%z.names) stop("Variable '",y.name,"' appears in both arguments 'y.name' and 'z.names'")
      auxchk3 <- intersect(z.names,x.names)
      if(length(auxchk3)>0) stop("Variable '",auxchk3[1],"' appears in both arguments 'x.names' and 'z.names'")
      } else {
      z.names <- NULL
      }
    }
  #
  if(sum(is.na(data[,c(y.name,x.names,z.names)]))>0) stop("Data contain missing values: you can use tsEM() to impute them")
  #
  if(!is.null(unit)) {
    unit <- unit[which(!is.na(unit))]
    if(!is.character(unit)|length(unit)!=1) stop("Argument 'unit' must be either NULL or a character vector of length 1")
    if(length(unit)>0) {
      if(length(setdiff(unit,colnames(data)))>0) stop("Variable '",unit,"' not found")
      if(length(intersect(unit,y.name))>0) stop("Variable '",unit,"' appears in both arguments 'y.name' and 'unit'")
      if(length(intersect(unit,x.names))>0) stop("Variable '",unit,"' appears in both arguments 'x.names' and 'unit'")
      if(length(intersect(unit,z.names))>0) stop("Variable '",unit,"' appears in both arguments 'z.names' and 'unit'")
      if(sum(is.na(data[,unit]))>0) stop("Variable '",unit,"' provided to argument 'unit' contains missing values")
      data[,unit] <- factor(data[,unit])
      } else {
      unit <- NULL
      }
    }
  #
  if(!is.null(time)) {
    time <- time[which(!is.na(time))]
    if(!is.character(time)|length(time)!=1) stop("Argument 'time' must be either NULL or a character vector of length 1")
    if(length(time)>0) {
      if(length(setdiff(time,colnames(data)))>0) stop("Variable '",time,"' not found")
      if(length(intersect(time,unit))>0) stop("Variable '",time,"' appears in both arguments 'unit' and 'time'")
      if(length(intersect(time,y.name))>0) stop("Variable '",time,"' appears in both arguments 'y.name' and 'time'")
      if(length(intersect(time,x.names))>0) stop("Variable '",time,"' appears in both arguments 'x.names' and 'time'")
      if(length(intersect(time,z.names))>0) stop("Variable '",time,"' appears in both arguments 'z.names' and 'time'")
      if(!is.numeric(data[,time])&!identical(class(data[,time]),"Date")) stop("Variable '",time,"' must be numeric or of class 'Date'")
      if(sum(is.na(data[,time]))>0) stop("Variable '",time,"' provided to argument 'time' contains missing values")
      if(is.null(unit)) {
        if(sum(duplicated(data[,time]))>0) stop("Variable '",time,"' contains duplicated values")
        } else {
        if(sum(sapply(split(data,data[,unit]),function(x){sum(duplicated(x[,time]))}))>0) stop("Variable '",time,"' contains duplicated values")
        }
      } else {
      time <- NULL  
      }
    }
  if(!is.null(time)) {
    if(is.null(unit)) {
      data <- data[order(data[,time]),]
      } else {
      data <- data[order(data[,unit],data[,time]),]
      }
    }
  #
  #data <- data[complete.cases(data[,c(unit,time,y.name,x.names,z.names)]),]
  dataD <- preProcess(var.names=c(y.name,x.names,z.names), unit=unit, time=time, data=data, box.cox=box.cox, ndiff=ndiff)  ## <--
  if(!is.null(weights)&max(ndiff)>0) weights <- weights[-(1:max(ndiff))]
  #
  add.intercept <- add.intercept[1]
  if(is.na(add.intercept)||(!is.logical(add.intercept)|is.null(add.intercept))) add.intercept <- T
  multi.intercept <- ifelse(add.intercept & attr(dataD,"ndiff")[y.name]>0, F, T)
  #
  quiet <- quiet[1]
  if(is.na(quiet)||(!is.logical(quiet)|is.null(quiet))) quiet <- FALSE
  nstart <- control$nstart[1]
  if(is.null(nstart)) nstart <- 1
  if(!is.numeric(nstart)) nstart <- 1 else nstart <- max(1,ceiling(nstart))
  if(nstart==1) quiet <- T
  #
  #max.try <- control$max.try
  #if(is.null(max.try)) max.try <- 50
  #max.start <- control$max.start
  #if(is.null(max.start)) max.start <- max.try*10
  #grid.by <- control$grid.by
  #if(is.null(grid.by)) grid.by <- 0.05
  #if(length(grid.by)>1) grid.by <- grid.by[1]
  #if(!is.numeric(grid.by) || grid.by<=0 || grid.by>0.1) stop("Argument 'grid.by' must be a positive value no greater than 0.1")
  #
  grid.by <- 0.05
  delta.lim <- control$delta.lim
  lambda.lim <- control$lambda.lim
  peak.lim <- control$peak.lim
  length.lim <- control$length.lim
  sign <- NULL  ## <----- constraint segno
  max.try <- max.start <- nstart
  #
  if(!is.numeric(offset)) {
    offset <- rep(0,length(x.names))
    names(offset) <- x.names
    }
  if(length(offset)==1&is.null(names(offset))) {
    offset <- rep(max(0,offset),length(x.names))
    names(offset) <- x.names
    } else {
    offset <- offset[x.names]
    names(offset) <- x.names
    offset[which(is.na(offset)|offset<0)] <- 0
    }
  if(y.name%in%x.names && offset[y.name]<1) offset[y.name] <- 1
  #
  if(!is.list(delta.lim)&length(delta.lim)==2) delta.lim <- optFormat(list(), x.names, delta.lim)
  delta.lim <- optFormat(delta.lim, x.names, c(0,1))    
  if(!is.list(lambda.lim)&length(lambda.lim)==2) lambda.lim <- optFormat(list(), x.names, lambda.lim)
  lambda.lim <- optFormat(lambda.lim, x.names, c(0,1))
  if(sum(sapply(delta.lim,function(x){x[1]==x[2]})==F)==0 &
      sum(sapply(lambda.lim,function(x){x[1]==x[2]})==F)==0) {
    par <- rbind(sapply(delta.lim,function(x){x[1]}),sapply(lambda.lim,function(x){x[1]}))
    } else {
    par <- NULL
    }
  #
  if(!is.null(par)) {
    modOK <- gam_olsFit(y.name=y.name, x.names=x.names, par=par, z.names=z.names, unit=unit, offset=offset, data=dataD, normalize=T, add.intercept=add.intercept, multi.intercept=multi.intercept, weights=weights)
    } else {
    p <- length(x.names)
    if(!is.list(peak.lim)&length(peak.lim)==2) peak.lim <- optFormat(list(), x.names, peak.lim)
    peak.lim <- optFormat(peak.lim, x.names, c(0,Inf))
    if(!is.list(length.lim)&length(length.lim)==2) length.lim <- optFormat(list(), x.names, length.lim)
    length.lim <- optFormat(length.lim, x.names, c(0,Inf))
    #if(is.null(sign)) sign <- rep(0,p)
    #
    if(quiet==F) {
      cat("Looking for valid models ...")
      flush.console()
      }
    gridList <- gam_parGrid(delta.lim=delta.lim, lambda.lim=lambda.lim, peak.lim=peak.lim, length.lim=length.lim, grid.by=grid.by)
    pcombtot <- prod(sapply(gridList,nrow))
    if(quiet==F) cat("\r","Found ",signif(pcombtot)," valid models        ",sep="","\n")
    if(is.null(gridList)) {
      if(quiet==F) {
        cat("\n")
        cat("No solution with the current constraints","\n")
        }
      } else {
      rssOK <- Inf
      ntry_fit <- stopped <- 0
      modOK <- NULL
      searchList <- visitList <- list()
      if(nstart==1) {
        ini0 <- c()
        if(quiet==F) {
          cat('\r',"Generating starting values ...")
          flush.console()
          }
        for(i in 1:length(x.names)) {
          igrid <- gridList[[i]]
          irss <- Inf
          for(j in 1:nrow(igrid)) {
            suppressWarnings(
              ijm <- gam_olsFit(y.name=y.name, x.names=x.names[i],
                z.names=c(setdiff(x.names,x.names[i]),z.names),
                unit=unit, data=dataD, offset=offset[x.names[i]],
                par=cbind(igrid[j,]), add.intercept=add.intercept, multi.intercept=multi.intercept, normalize=F, weights=weights)
              )
            ijrss <- sum(ijm$residuals^2)
            if(ijrss<irss) {
              irss <- ijrss
              ipar <- ijm$par
              }
            }
          ini0 <- cbind(ini0, ipar)
          }
        if(quiet==F) {
          cat('\r',"Running hill climbing ...     ")
          flush.console()
          }
        gs0 <- gam_hcFun(y.name=y.name, x.names=x.names, z.names=z.names, unit=unit, data=dataD,
          offset=offset, inits=ini0, visitList=visitList, gridList=gridList,
          sign=sign, grid.by=grid.by, add.intercept=add.intercept, multi.intercept=multi.intercept, weights=weights)
        modOK <- gs0$model
        if(quiet==F) {
          cat('\r',"Explored ",length(gs0$par.tested)," valid models       ",sep="")
          }
        }
      #
      else
      #
      for(i in 1:max.start) {
        ini0 <- gam_inits(gridList=gridList, visitList=visitList, maxtry=max.start)
        if(!is.null(ini0)) {
          gs0 <- gam_hcFun(y.name=y.name, x.names=x.names, z.names=z.names, unit=unit, data=dataD,
            offset=offset, inits=ini0, visitList=visitList, gridList=gridList,
            sign=sign, grid.by=grid.by, add.intercept=add.intercept, multi.intercept=multi.intercept, weights=weights)
          } else {
          gs0 <- NULL
          stopped <- 2
          break()
          }
        visitList <- c(visitList,gs0$par.tested)
        if(quiet==F) {
          cat('\r',"Restart ",i,"/",nstart,": explored ",length(visitList)," valid models",sep="")
          flush.console()
          }
        m0 <- gs0$model
        if(is.null(m0)) {
          rss0 <- Inf
          } else {
          rss0 <- sum(m0$residuals^2)
          searchList <- c(searchList,list(m0))
          }
        if(rss0<rssOK) {
          modOK <- m0
          rssOK <- rss0
          ntry_fit <- 0
          } else {
          ntry_fit <- ntry_fit+1
          }
        if(ntry_fit>=max.try) {
          stopped <- 1
          break()
          }
        }
      if(quiet==F) {
        cat("\n")
        if(stopped==2) cat("No more valid models found. End",sep="","\n")
        #if(stopped==1) {
        #  cat("No improvement found in the last ",max.try, " restarts. End",sep="","\n")
        #  } else if(stopped==2) {
        #    cat("No more valid starting values found. End",sep="","\n")
        #  } else {
        #  cat("Maximum number of restarts reached. End",sep="","\n")
        #  }
        }
      #colnames(visitList[[1]]) <- x.names
      #rownames(visitList[[1]]) <- c("delta","lambda")
      #modOK$inits <- visitList
      rssVal <- sapply(searchList,function(x){sum(x$residuals^2)})
      searchOK <- searchList[order(rssVal)]
      modOK$local.max <- searchOK
      }
    }
  modOK$variables <- list(y.name=y.name,x.names=x.names,z.names=z.names,unit=unit,time=time)
  if(is.null(unit)) {
    idg <- NULL
    } else {
    idg <- lapply(split(dataD, dataD[,unit]), rownames)  
    }
  modOK$control <- list(nstart=nstart,delta.lim=delta.lim,lambda.lim=lambda.lim,peak.lim=peak.lim,length.lim=length.lim)
  modOK$unit.id <- idg
  modOK$box.cox <- attr(dataD,"box.cox")
  modOK$ndiff <- attr(dataD,"ndiff")
  modOK$lag.order <- lagOrderCalc(modOK$residuals,modOK$unit.id)
  modOK$data.orig <- data[,c(unit,time,y.name,x.names,z.names)]
  modOK$data.used <- dataD[,c(unit,time,y.name,x.names,z.names)]
  if(is.null(unit)) {
    pval <- oneTest(modOK$residuals)$p.value
    } else {
    pval <- oneTest(modOK$residuals, unit=dataD[,unit])$p.value["(combined)",]
    }
  if(pval["adf"]>0.05) warning("ADF test on residuals is not significant: regression could be spurious", call.=F)
  if(pval["kpss"]<0.05) warning("KPSS test on residuals is significant: regression could be spurious", call.=F)
  modOK
  }

# get the number of intercepts (auxiliary)
nInterc <- function(object) {
  if(object$add.intercept) {
    if(is.null(object$variables$unit)) {
      1
      } else {
      #ifelse(multi.intercept, length(object$unit.id), 1)
      ifelse(object$ndiff[object$variables$y.name]==0, length(object$unit.id), 1)
      }
    } else {
    0
    }
  }

# function for p-value notation (auxiliary)
pStars <- function(p) {
  res <- rep("",length(p))
  res[which(p<0.1)] <- "."
  res[which(p<0.05)] <- "*"
  res[which(p<0.01)] <- "**"
  res[which(p<0.001)] <- "***"
  res
  }

# summary method for class 'gammadlm'
summary.gammadlm <- function(object, ...) {
  summ <- summary.lm(object, ...)
  ttab <- summ$coefficients
  SS <- hacCalc(Xmat=model.matrix(object), resid=object$residuals, uS=summary.lm(object)$cov.unscaled, unitID=object$unit.id, lagsel=object$lag.order)
  bse <- sqrt(diag(SS))
  ttab[,2] <- bse
  ttab[,3] <- ttab[,1]/ttab[,2]
  ttab[,4] <- round(2*pt(-abs(ttab[,3]),object$df.residual),6)
  #
  ttab <- data.frame(ttab,pStars(ttab[,4]))
  colnames(ttab) <- c("Estimate","S.E.","t value","Pr(>|t|)","")
  n_alpha <- nInterc(object)
  if(n_alpha>0) alphaTab <- ttab[1:n_alpha,,drop=F] else alphaTab <- NULL
  nomi <- object$variables$x.names
  quan <- matrix(nrow=ncol(object$par),ncol=5)
  for(i in 1:ncol(object$par)) {
    quan[i,] <- c(peakCalc(object$par[,i])+object$offset[i],
                  gammaQuantile(c(.5,.95,.99,.999),object$par[,i],object$offset[i]))
    }
  rownames(quan) <- colnames(object$par)
  colnames(quan) <- c("peak","50%","95%","99%","99.9%")
  par <- cbind(t(object$par),offset=object$offset,round(quan,1))
  xTab <- ttab[(n_alpha+1):(n_alpha+length(nomi)),]
  rownames(xTab) <- nomi
  colnames(xTab)[1:2] <- c("theta","S.E.(theta)")
  if(!is.null(object$variables$z.names)) {
    nomiZ <- rownames(ttab)[sapply(paste0("^",object$variables$z.names),grep,rownames(ttab))]
    zTab <- ttab[nomiZ,,drop=F]
    } else {
    zTab <- NULL
    }
  if(!is.null(object$unit.id)) {
    n_unit <- length(object$unit.id)
    } else {
    n_unit <- 1
    }
  pred <- object$fitted.values
  obs <- object$data.used[,object$variables$y.name]
  err <- c(rmse=sqrt(mean((obs-pred)^2)),
           mae=mean(abs(obs-pred)),
           mape=mean(abs((obs-pred)/obs)))
  summ <- list(formula=object$call$formula,
               variables=object$variables,
               n_unit=n_unit, n=nobs(object),
               intercept=alphaTab, par=par, x=xTab, z=zTab,
               sigma=summ$sigma, df=object$df,
               fstat=c(summ$fstatistic,p=1-pf(summ$fstatistic[1],summ$fstatistic[2],summ$fstatistic[3])),
               rsq=summ$r.squared, error=err)
  class(summ) <- "summary.gammadlm"
  summ
  }

# print method for class 'summary.gammadlm'
print.summary.gammadlm <- function(x, ...) {
  cat("Gamma distributed-lag model estimated through hill-climbing","\n")
  cat("  Number of units: ",x$n_unit,sep="","\n")
  cat("  Average number of time points per unit: ",x$n/x$n_unit,"\n",sep="")
  cat("  Response variable: ",x$variables$y.name,"\n",sep="")
  cat("\n")
  cat("Gamma lag distributions","\n")
  print(x$par)
  print(x$x)
  cat("\n")
  if(!is.null(x$z)) {
    cat("Explanatory variables without lags","\n")
    print(x$z)
    cat("\n")
    }
  if(!is.null(x$intercept)) {
    cat("Intercepts","\n")
    print(x$intercept)
    cat("\n")
    }
  cat("In-sample prediction error","\n")
  print(x$error)
  cat("\n")
  cat("Residual std. error: ",round(x$sigma,6)," on ",x$df," DFs (R-squared: ",round(x$rsq,6),")","\n",sep="")
  cat("F statistic: ",round(x$fstat[1],6)," on ",x$fstat[2]," and ",x$fstat[3]," DFs (p-value: ",round(x$fstat[4],4),")","\n",sep="")  
  }

# vcov method for class 'gammadlm'
vcov.gammadlm <- function(object, ...) {
  hacCalc(Xmat=model.matrix(object), resid=object$residuals, uS=summary.lm(object)$cov.unscaled, unitID=object$unit.id, lagsel=object$lag.order)
  }

# asterisk notation (auxiliary)
starFun <- function(x) {
  res <- rep("",length(x))
  res[which(x<0.1)] <- "."
  res[which(x<0.05)] <- "*"
  res[which(x<0.01)] <- "**"
  res[which(x<0.001)] <- "***"
  res
  }

# selection of lag order (auxiliary)
lagSelect <- function(y, x=NULL, order.max=NULL) {
  n <- length(y)
  if(is.null(order.max)) order.max <- sqrt(n)
  if(is.null(x)) {
    Xlag <- rep(1,n)
    for(i in 1:order.max) Xlag <- cbind(Xlag, c(rep(NA,i),y[1:(n-i)]))
    } else {
    Xlag <- x
    for(i in 1:order.max) Xlag <- cbind(Xlag, c(rep(NA,i),x[1:(n-i)]))
    }
  Xlag[1:order.max,] <- NA
  bic <- c()
  for(i in 1:ncol(Xlag)) {
    bic[i] <- BIC(lm(y~Xlag[,1:i,drop=F]))
    }
  which.min(bic)-1
  }

# lag order of residuals (auxiliary)
lagOrderCalc <- function(resid, unitID) {
  if(is.null(unitID)) {
    #max.lag <- ar(resid)$order
    max.lag <- lagSelect(resid)
    max.lag2 <- lagSelect(resid^2)
    c('residuals'=max.lag,'residuals^2'=max.lag2)
    } else {
    max.lag <- max.lag2 <- c()
    for(i in 1:length(unitID)) {
      #max.lag[i] <- ar(resid[unitID[[i]]])$order
      max.lag[i] <- lagSelect(resid[unitID[[i]]])
      max.lag2[i] <- lagSelect(resid[unitID[[i]]]^2)
      }
    names(max.lag) <- names(max.lag2) <- names(unitID)
    cbind('residuals'=max.lag,'residuals^2'=max.lag2)
    }
  }

# hac covariance matrix (auxiliary)
hacCalc <- function(Xmat, resid, uS, unitID, lagsel) {
  xdel <- setdiff(colnames(Xmat),colnames(uS))
  if(length(xdel)>0) Xmat <- Xmat[,colnames(uS)]
  if(is.null(unitID)) {
    max.lag <- lagsel[1]
    max.lag2 <- lagsel[2]
    } else {
    max.lag <- lagsel[,1]
    max.lag2 <- lagsel[,2]
    }
  if(sum(max.lag)+sum(max.lag2)==0) {
    s2 <- sum(resid^2)/(length(resid)-ncol(Xmat))
    W <- diag(rep(s2,length(resid)))
    } else {
    W <- diag(resid^2)
    rownames(W) <- colnames(W) <- names(resid)
    if(is.null(unitID)) {
      n <- length(resid)
      if(max.lag>0) {
        for(i in 1:(n-max.lag)) {
          for(j in (i+1):(i+max.lag)) {
            W[i,j] <- W[j,i] <- resid[i]*resid[j]*(1-abs(i-j)/(max.lag+1))
            }
          }
        }
      } else {
      W <- diag(resid^2)
      rownames(W) <- colnames(W) <- names(resid)
      for(w in 1:length(unitID)) {
        ind <- unitID[[w]]
        wn <- length(ind)
        wres <- resid[ind]
        wei <- matrix(0,nrow=length(ind),ncol=length(ind))
        diag(wei) <- wres^2
        if(max.lag[w]>0) {
          for(i in 1:(wn-max.lag[w])) {
            for(j in (i+1):(i+max.lag[w])) {
              wei[i,j] <- wei[j,i] <- wres[i]*wres[j]*(1-abs(i-j)/(max.lag[w]+1))
              }
            }
          }
        W[ind,ind] <- wei
        }
      }
    }
  uS%*%t(Xmat)%*%W%*%Xmat%*%uS
  }

# extract lag coefficients
lagCoef <- function(x, conf=0.95, cumulative=FALSE, max.lag=NULL, max.quantile=0.999) {
  if(missing(x)) stop("Argument 'x' is missing")
  if(!identical(class(x),c("gammadlm","lm"))) stop("Argument 'x' must be an object of class 'gammadlm'")
  #
  cumulative <- cumulative[1]
  if(is.na(cumulative)||(!is.logical(cumulative)|is.null(cumulative))) cumulative <- FALSE
  #
  if(!is.numeric(conf)) {
    conf <- NULL
    } else {
    conf <- conf[1]
    if(conf<=0|conf>=1) conf <- NULL
    }
  #
  max.lag <- max.lag[1]
  if(!is.numeric(max.lag)) max.lag <- NULL else max.lag <- max(0,ceiling(max.lag))
  #
  max.quantile <- max.quantile[1]
  if(!is.numeric(max.quantile)) max.quantile <- 0.999 else max.quantile <- min(max(1e-12,max.quantile),1-1e-12)
  #
  gpar <- x$par
  offs <- x$offset
  p <- ncol(gpar)
  n_alpha <- nInterc(x)
  Smat <- vcov(x)
  res <- vector("list",length=p)
  names(res) <- x$variables$x.names
  for(i in 1:p) {
    if(is.null(max.lag)) {
      imaxlag <- gammaQuantile(prob=max.quantile, par=gpar[,i], offset=offs[i])
      } else {
      imaxlag <- max.lag
      }
    ires <- data.frame(matrix(nrow=imaxlag+1,ncol=2))
    lagwei <- gammaWeights(0:imaxlag, par=gpar[,i], offset=offs[i], normalize=T)
    ires[,1] <- x$coef[i+n_alpha]*lagwei
    ires[,2] <- sqrt(Smat[i+n_alpha,i+n_alpha])*lagwei
    rownames(ires) <- 0:imaxlag
    colnames(ires) <- c("Estimate","Std. Error")
    if(cumulative==F) {
      res[[i]] <- ires
      } else {
      iresC <- ires
      iresC[,1] <- cumsum(ires[,1])
      #
      for(j in 1:nrow(ires)) {
        ijse <- ires[1:j,2]
        iresC[j,2] <- sqrt(sum(ijse%*%t(ijse)))
        }
      res[[i]] <- iresC
      }
    }
  if(!is.null(conf)) {
    zval <- qnorm((1+conf)/2)
    lapply(res, function(x) {
      x[,paste0(round(100*conf,1),"%_lower")] <- x[,1]-zval*x[,2]
      x[,paste0(round(100*conf,1),"%_upper")] <- x[,1]+zval*x[,2]
      x
      })
    } else {
    res
    }
  }

# plot method for class 'gammadlm'
plot.gammadlm <- function(x, x.names=NULL, conf=0.95, max.lag=NULL, max.quantile=0.999,
  xlim=NULL, ylim=NULL, add.legend=TRUE, cex.legend=1, digits=4, grid.length=100, main=NULL, ylab=NULL, xlab=NULL, ...) {
  #
  if(!is.numeric(conf)) {
    conf <- NULL
    } else {
    conf <- conf[1]
    if(conf<=0|conf>=1) conf <- NULL
    }
  #
  max.lag <- max.lag[1]
  if(!is.numeric(max.lag)) max.lag <- NULL else max.lag <- max(0,ceiling(max.lag))
  #
  max.quantile <- max.quantile[1]
  if(!is.numeric(max.quantile)) max.quantile <- 0.999 else max.quantile <- min(max(1e-12,max.quantile),1-1e-12)
  #
  add.legend <- add.legend[1]
  if(is.na(add.legend)||(!is.logical(add.legend)|is.null(add.legend))) add.legend <- TRUE
  if(add.legend) {
    grid.length <- grid.length[1]
    if(!is.numeric(grid.length)) grid.length <- 100 else grid.length <- max(100,grid.length)
    digits <- digits[1]
    if(!is.numeric(digits)) digits <- 4 else digits <- max(digits,1)
    cex.legend <- cex.legend[1]
    if(!is.numeric(cex.legend)) cex.legend <- 1 else cex.legend <- max(cex.legend,1e-12)
    }
  #
  makePlot <- function(i,main) {
    gpar <- x$par[,i]
    offs <- x$offset[i]
    auxlen <- ceiling(gammaQuantile(prob=max.quantile, par=gpar, offset=offs))
    if(is.null(max.lag)) {
      laglen <- auxlen
      } else {
      laglen <- max.lag
      }
    if(is.null(xlim)) xlim <- c(-1,laglen)
    lseq <- seq(xlim[1],xlim[2],length=grid.length)
    lseq <- sort(c(offs,offs-1e-12,lseq))
    lagwei <- gammaWeights(lseq, par=gpar, offset=offs, normalize=T)
    if(sum(lagwei==Inf)>0) {
      lagwei[which(lagwei==Inf)] <- 0
      auxsx <- which(lseq==max(lseq[which(lseq<=offs)]))
      auxdx <- which(lseq==min(lseq[which(lseq>=offs)]))
      lagwei <- c(lagwei[1:auxsx],1,lagwei[auxdx:(length(lagwei)-1)])
      }
    bcoef <- x$coef[i+n_alpha]*lagwei
    bcoef_se <- sqrt(Smat[i+n_alpha,i+n_alpha])*lagwei
    auxinf <- which(abs(bcoef_se)==Inf|is.na(abs(bcoef_se)))
    if(length(auxinf)>0) {
      auxbse <- bcoef_se
      auxbse[auxinf] <- NA
      #auxpred <- approx(auxbse,n=length(auxbse))$y
      auxpred <- spline(1:length(auxbse),auxbse,xout=1:length(auxbse))$y
      bcoef_se[auxinf] <- auxpred[auxinf]
      }
    #
    if(is.null(conf)) {
      lagco <- cbind(bcoef)
      } else {
      lagco <- cbind(bcoef,bcoef-tquan*bcoef_se,bcoef+tquan*bcoef_se)
      }
    if(is.null(ylim)) ylim <- range(lagco)
    plot(lseq, lagco[,1], type="n", ylim=ylim, ylab=ylab, xlab=xlab, main=main, ...)
    if(!is.null(conf)) {
      polygon(c(lseq,rev(lseq)), c(lagco[,2],rev(lagco[,3])),
              border=NA, col=adjustcolor("grey75",alpha.f=0.5))
      }
    grid()
    #abline(h=0,lty=3,col=1)
    lines(lseq,lagco[,1])
    if(add.legend) {
      bcum <- x$coef[i+n_alpha]
      if(is.null(conf)) {
        legtxt <- paste("Multiplier up to ",auxlen," lags: ",round(signif(bcum),digits))
        } else {
        bcum_se <- sqrt(Smat[i+n_alpha,i+n_alpha])
        bcumco <- c(bcum,bcum-tquan*bcum_se,bcum+tquan*bcum_se)
        bcumcoOK <- signif(bcumco)
        legtxt <- paste("Multiplier up to ",auxlen," lags: ",round(bcumcoOK[1],digits),"\n",
          "   ",100*conf,"% CI: (",round(bcumcoOK[2],digits),", ",round(bcumcoOK[3],digits),")",sep="")
        }
      legend("topright",legend=legtxt,cex=cex.legend,bty="n")
      }
    box()
    }
  xnam <- x$variables$x.names
  if(!is.null(x.names)) {
    chk0 <- setdiff(x.names,xnam)
    #if(length(chk0)>0) warning("Variables: ",paste(chk0,collapse=", ")," not found",call.=F)
    xOK <- setdiff(x.names,chk0)
    if(length(xOK)==0) xOK <- xnam
    } else {
    xOK <- xnam
    }
  if(is.null(ylab)) ylab <- "Multiplier"
  if(is.null(xlab)) xlab <- "Time lag"
  Smat <- vcov(x)
  n_alpha <- nInterc(x)
  if(is.null(main)) main <- xOK
  #if(!is.null(conf)) tquan <- qt((1+conf)/2,x$df.residual)
  if(!is.null(conf)) tquan <- qnorm((1+conf)/2)
  oldpar <- par(no.readonly=T)
  on.exit(par(oldpar))
  par(mfrow=n2mfrow(length(xOK)))
  for(i in 1:length(xOK)) {
    ind <- which(xnam==xOK[i])
    if(length(main)>=i) imain <- main[i] else imain <- ""
    makePlot(ind,imain)
    }
  }

# fitted method for class 'gammadlm'
fitted.gammadlm <- function(object, ...) {
  if(is.null(object$variables$unit)) {
    tab <- data.frame(object$data.used[,object$variables$time], object$fitted.values)
    colnames(tab) <- c(object$variables$time, object$variables$y.name)
    } else {
    tab <- data.frame(object$data.used[,c(object$variables$unit,object$variables$time)], object$fitted.values)
    colnames(tab) <- c(object$variables$unit, object$variables$time, object$variables$y.name)
    }
  tab
  }

# residuals method for class 'gammadlm'
residuals.gammadlm <- function(object, ...) {
  if(is.null(object$variables$unit)) {
    tab <- data.frame(object$data.used[,object$variables$time], object$residuals)
    colnames(tab) <- c(object$variables$time, object$variables$y.name)
    } else {
    tab <- data.frame(object$data.used[,c(object$variables$unit,object$variables$time)], object$residuals)
    colnames(tab) <- c(object$variables$unit, object$variables$time, object$variables$y.name)
    }
  tab
  }

# generate lags (auxiliary)
LAG <- function(x, p, unit=NULL, ...) {
  if(p>0) {
    #
    lfun <- function(v) {
      #v0 <- rep(0,p)
      v0 <- rep(mean(v[1:p],na.rm=T),p)  ## <--
      n <- length(v)
      res <- matrix(nrow=n,ncol=p)
      for(i in 1:p) {
        #res[,i] <- c(rep(NA, i), v[1:(n-i)])
        res[,i] <- c(v0[(p-i+1):p], v[1:(n-i)])
        }
      colnames(res) <- 1:p
      res
      }
    #
    if(is.null(unit)) {
      lfun(x)
      } else {
      res <- matrix(nrow=length(x),ncol=p)
      gr <- levels(factor(unit))
      for(i in 1:length(gr)) {
        ind <- which(unit==gr[i])
        res[ind,] <- lfun(x[ind])
        }
      colnames(res) <- 1:p
      res
      }
    } else {
    x  
    }
  }

# find cross-correlation order (auxiliary)
crossCorOrder <- function(var.names, unit=NULL, data, maxlag) {
  corMat <- array(dim=c(length(var.names),length(var.names),maxlag+1))
  dimnames(corMat) <- list(var.names,var.names,0:maxlag)
  if(is.null(unit)) {
    res <- data[,var.names,drop=F]
    uvar <- NULL
    } else {
    datlist <- split(data[,var.names,drop=F],data[,unit])
    res <- do.call(rbind,lapply(datlist,function(x){apply(x,2,function(z){z-mean(z,na.rm=T)})}))
    uvar <- data[,unit]
    }
  for(i in 1:length(var.names)) {
    for(j in i:length(var.names)) {
      y <- res[,i]
      x <- res[,j]
      lx <- cbind(x,LAG(x,maxlag,uvar))
      corMat[i,j,] <- cor(y,lx,use="pairwise.complete.obs")[1,]
      }
    }
  corMat
  }

# function for EM imputation
tsEM <- function(var.names, unit=NULL, time=NULL, data, nlags=NULL, tol=1e-4, maxit=1000, quiet=FALSE) {
  #
  if(missing(data)) stop("Argument 'data' is missing")
  if(!identical(class(data),"data.frame")) stop("Argument 'data' must be a data.frame")
  if(missing(var.names)) stop("Argument 'var.names' is missing")
  if(!is.character(var.names)) {
    stop("Argument 'var.names' must be a character vector")
    } else {
    var.names[which(!is.na(var.names))]
    if(length(var.names)<1) stop("Argument 'var.names' must be a character vector of length 1 or greater")
    }
  auxchk <- setdiff(var.names,colnames(data))  
  if(length(auxchk)>0) stop("Variable '",auxchk[1],"' not found")
  xnum <- xcat <- c()
  for(i in 1:length(var.names)) {
    if(isQuant(data[,var.names[i]])) {
      xnum <- c(xnum,var.names[i])
      } else {
      if(sum(is.na(data[,var.names[i]]))>0) stop("Variable '",var.names[i],"' is categorical and contains missing values")
      xcat <- c(xcat,var.names[i])
      }
    }
  #
  unit <- unit[1]
  if(!is.null(unit)&&is.na(unit)) unit <- NULL
  if(!is.null(unit)) {
    if(length(setdiff(unit,colnames(data)))>0) stop("Variable '",unit,"' not found")
    if(length(intersect(unit,var.names))>0) stop("Variable '",unit,"' appears in both arguments 'var.names' and 'unit'")
    if(sum(is.na(data[,unit]))>0) stop("Variable '",unit,"' provided to argument 'unit' contains missing values")
    data[,unit] <- factor(data[,unit])
    }
  #
  time <- time[1]
  if(!is.null(time)&&is.na(time)) time <- NULL
  if(!is.null(time)) {
    if(length(setdiff(time,colnames(data)))>0) stop("Variable '",time,"' not found")
    if(length(intersect(time,var.names))>0) stop("Variable '",time,"' appears in both arguments 'var.names' and 'time'")
    if(length(intersect(time,unit))>0) stop("Variable '",time,"' appears in both arguments 'unit' and 'time'")
    if(!is.numeric(data[,time])&!identical(class(data[,time]),"Date")) stop("Variable '",time,"' must be numeric or of class 'Date'")
    if(sum(is.na(data[,time]))>0) stop("Variable '",time,"' provided to argument 'time' contains missing values")
    if(is.null(unit)) {
      if(sum(duplicated(data[,time]))>0) stop("Variable '",time,"' contains duplicated values")
      } else {
      if(sum(sapply(split(data,data[,unit]),function(x){sum(duplicated(x[,time]))}))>0) stop("Variable '",time,"' contains duplicated values")
      }
    }
  #
  quiet <- quiet[1]
  if(is.na(quiet)||(!is.logical(quiet)|is.null(quiet))) quiet <- FALSE
  if(!is.numeric(tol)|is.null(tol)) tol <- 1e-4
  if(tol<=0) tol <- 1e-4
  if(!is.numeric(maxit)|is.null(maxit)) maxit <- 1000
  if(maxit<=0) maxit <- 1000 else maxit <- ceiling(maxit)
  #
  if(!is.null(time)) {
    if(is.null(unit)) {
      data <- data[order(data[,time]),]
      } else {
      data <- data[order(data[,unit],data[,time]),]
      }
    }
  if(!is.null(unit)) {
    data[,unit] <- factor(data[,unit])
    if(nlevels(data[,unit])<=1) unit <- NULL
    }
  if(is.null(unit)) {
    n <- nrow(data)
    } else {
    n <- min(sapply(split(data, data[,unit]), nrow))
    }
  dataI <- data
  lambda <- c()
  for(i in 1:length(xnum)) {
    if(sum(data[,xnum[i]]<=0,na.rm=T)==0) {
      lambda[i] <- 0
      dataI[,xnum[i]] <- log(data[,xnum[i]]) 
      } else if(sum(data[,xnum[i]]<0,na.rm=T)==0) {
      lambda[i] <- 0.5
      dataI[,xnum[i]] <- sqrt(data[,xnum[i]]) 
      } else {
      lambda[i] <- 1  
      }
    }
  names(lambda) <- xnum
  isNA <- list()
  for(i in 1:length(xnum)) {
    isNA[[i]] <- which(is.na(dataI[,xnum[i]]))
    dataI[isNA[[i]],xnum[i]] <- mean(dataI[,xnum[i]],na.rm=T) 
    }
  names(isNA) <- xnum
  #
  nlags <- nlags[1]
  if(!is.null(nlags)&&is.na(nlags)) nlags <- NULL
  if(!is.null(nlags)&&!is.numeric(nlags)) nlags <- NULL
  if(is.null(nlags)) {
    cmat <- crossCorOrder(xnum, unit=unit, data=data, maxlag=sqrt(n))
    cut <- qnorm(0.975)/sqrt(n)
    k0 <- quantile(which(abs(cmat)>cut,arr.ind=T)[,3],prob=0.9)-1
    if(!is.numeric(k0)) k0 <- sqrt(n)
    } else {
    k0 <- max(0,nlags)
    }
  nlags <- round(min(k0, 0.5*(nrow(data)-1-length(xcat))/length(xnum)))
  if(quiet==F) cat("Selected ",nlags," lag orders","\n",sep="")
  #
  ll <- -Inf
  if(is.null(unit)) {
    fstr <- ""
    ustr <- ""
    } else {
    fstr <- paste0(unit,"+")
    ustr <- paste0(",",unit)
    }
  if(length(xcat)>0) {
    zstr <- paste0("+",paste(xcat,collapse="+"))
    } else {
    zstr <- ""  
    }
  fine <- ind <- 0
  if(quiet==F) cat("EM iteration 0. Log likelihood: -")
  flush.console()
  while(fine==0) {
    ind <- ind+1
    data_new <- dataI
    ll0 <- 0
    for(i in 1:length(xnum)) {
      if(nlags==0) {
        if(length(xnum)==1) {
          x0str <- 1
          } else {
          x0str <- paste0("LAG(",xnum[-i],",",nlags,ustr,")")
          }
        iform <- formula(paste0(xnum[i],"~",fstr,paste(x0str,collapse="+"),zstr))
        } else {
        iform <- formula(paste0(xnum[i],"~",fstr,paste(paste0("LAG(",xnum,",",nlags,ustr,")"),collapse="+"),zstr))
        }
      imod <- lm(iform,data=dataI)
      imod_rev <- lm(iform,data=dataI[nrow(dataI):1,,drop=F])
      #
      #ll0 <- ll0+logLik(imod)[1]
      ll0 <- ll0+(logLik(imod)[1]+logLik(imod_rev)[1])/2
      #
      aux <- isNA[[xnum[i]]]
      if(length(aux)>0) {
        suppressWarnings(
          #data_new[aux,xnum[i]] <- predict(imod)[aux]
          data_new[aux,xnum[i]] <- (predict(imod)[aux]+rev(predict(imod_rev))[aux])/2
          )
        }
      }
    if(ll0>=ll) {
      if((ll0-ll)<tol | ind>=maxit) fine <- 1
      if(quiet==F) {
        cat("\r","EM iteration ",ind,". Log likelihood: ",ll0,sep="")
        flush.console()
        }
      ll <- ll0
      dataI <- data_new
      } else {
      #warning("Likelihood has decreased",call.=F)
      ind <- ind-1
      fine <- 1  
      }
    }
  if(quiet==F) {
    cat("\n")
    if(ind<maxit) {
      cat("Converged after ",ind," iterations",sep="")
      } else {
      cat("Convergence not achieved. Try to increase 'maxit' or decrease 'tol'")
      }
    cat("\n")
    }
  for(i in 1:length(xnum)) {
    if(lambda[xnum[i]]==0) {
      dataI[,xnum[i]] <- exp(dataI[,xnum[i]])
      } else if(lambda[xnum[i]]==0.5) {
      dataI[,xnum[i]] <- dataI[,xnum[i]]^2 
      }
    }
  dataI
  }
