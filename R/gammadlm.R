### DA FARE
#
# - residui e fitted values nel caso panel
#
# - gam_inits(): migliorare efficienza
# - gammadlm(): random restarts overdispersi
# - adfTest(): gestione missing
# - predict(): tenere conto dell'autocorrelazione degli errori
# - constraint segno
# - draw sample
# - full covariance matrix
#


# apply box-cox transformation (auxiliary)
makeBoxCox <- function(x, par) {
  if(par==1) {
    x
    } else if(par==0) {
    log(x)  
    } else {
    (x^par-1)/par  
    }
  }

# unconstrained kernel (auxiliary)
unconsKernel <- function(x, nlag, imputation=F) {
  if(nlag>0) {
    n <- length(x)
    res <- x
    for(i in 1:nlag) {
      ilx <- c(rep(NA,i),x[1:(n-i)])
      res <- cbind(res,ilx)
      }
    colnames(res) <- 0:nlag
    if(imputation==T) res[which(is.na(res))] <- mean(x,na.rm=T)
    } else {
    res <- matrix(x,ncol=1)
    colnames(res) <- 0
    }
  res
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
  if(length(par)!=2) stop("Argument 'par' must be of length 2")
  if(!is.numeric(par) || (par[1]<0 | par[1]>=1 | par[2]<0 | par[2]>=1)) stop("Both components of argument 'par' must be values >=0 and <1")
  #
  if(length(offset)>1) offset <- offset[1]
  if(!is.numeric(offset)) stop("Argument 'offset' must be a numerical value")
  #
  if(length(normalize)>1) normalize <- normalize[1]
  if(!is.logical(normalize)) stop("Argument 'normalize' must be a logical value")
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
  if(length(par)!=2) stop("Argument 'par' must be of length 2")
  if(!is.numeric(par) || (par[1]<0 | par[1]>=1 | par[2]<0 | par[2]>=1)) stop("Both components of argument 'par' must be values >=0 and <1")
  #
  if(length(offset)>1) offset <- offset[1]
  if(!is.numeric(offset)) stop("Argument 'offset' must be a numerical value")
  #
  if(length(normalize)>1) normalize <- normalize[1]
  if(!is.logical(normalize)) stop("Argument 'normalize' must be a logical value")
  #
  if(is.null(unit)) {
    gamkern(x=x, par=par, offset=offset, normalize=normalize)
    } else {
    gr <- unique(na.omit(unit))
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
  if(!is.numeric(prob) || sum(prob<=0)>0 | sum(prob>=1)>0) stop("Argument 'prob' must be a vector of values >0 and <1")
  #
  if(missing(par)) stop("Argument 'par' is missing")
  if(length(par)!=2) stop("Argument 'par' must be of length 2")
  if(!is.numeric(par) || (par[1]<0 | par[1]>=1 | par[2]<0 | par[2]>=1)) stop("Both components of argument 'par' must be values >=0 and <1")
  #
  if(length(offset)>1) offset <- offset[1]
  if(!is.numeric(offset)) stop("Argument 'offset' must be a numerical value")
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

# perform ADF test
adfTest <- function(var.names=NULL, unit=NULL, time=NULL, data, box.cox=1, ndiff=0, max.lag=NULL) {
  if(!identical(class(data),"data.frame")) stop("Argument 'data' must be a data.frame")
  #
  is.dummy <- function(x) {
    dom <- sort(unique(na.omit(x)))
    length(dom)==2&dom[1]==0&dom[2]==1
    }
  if(is.null(var.names)) {
    var.names <- setdiff(colnames(data),time)
    if(length(var.names)==0) stop("No quantitative variable found")
    x2del <- c()
    for(i in 1:length(var.names)) {
      if(!is.numeric(data[,var.names[i]])|is.dummy(data[,var.names[i]])) x2del <- c(x2del,var.names[i])
      }
    if(length(x2del)>0) var.names <- setdiff(var.names,x2del)
    if(length(var.names)==0) stop("No quantitative variable found")
    } else {
    auxchk <- setdiff(var.names,colnames(data))  
    if(length(auxchk)>0) stop("Unknown variable '",auxchk[1],"' in argument 'var.names'")
    #var.names <- setdiff(intersect(var.names,colnames(data)),time)
    }
  #
  if(!is.null(time)) {
    if(length(time)>1) time <- time[1]
    if((time%in%colnames(data))==F) stop("Unknown variable '",time,"' in argument 'time'")
    if(time%in%var.names) stop("Variable '",time,"' appears in both arguments 'var.names' and 'time'")
    if(!is.numeric(data[,time])&!identical(class(data[,time]),"Date")) stop("Variable '",time,"' must be numeric or of class 'Date'")
    }
  #
  testList <- list()
  dataD <- tsDiff(var.names=var.names, unit=unit, time=time, data=data, box.cox=box.cox, ndiff=ndiff)
  if(is.null(unit)) gr <- NULL else gr <- dataD[,unit]
  for(i in 1:length(var.names)) {
    iadf <- adfOneTest(x=dataD[,var.names[i]], unit=gr, max.lag=max.lag)
    iadf$box.cox <- unname(attr(dataD,"box.cox")[var.names[i]])
    iadf$ndiff <- unname(attr(dataD,"ndiff")[var.names[i]])
    testList[[i]] <- iadf
    }
  names(testList) <- var.names
  class(testList) <- "adfTest"
  testList
  }

# print method for class 'adfTest'
print.adfTest <- function(x, ...) {
  #cat("Box-Cox parameters:","\n")  
  #print(sapply(x,function(z){z$box.cox}))
  #cat("Number of differences:","\n")
  #print(sapply(x,function(z){z$ndiff}))
  cat("p-values (null hypothesis is unit root):","\n")
  print(sapply(x, function(z){
    pval <- round(z$p.value,4)
    if(length(pval)==1) pval else pval["(combined)"]
    }))
  }

# adf test for one variabile (auxiliary)
adfOneTest <- function(x, unit=NULL, max.lag=NULL) {
  if(missing(x)) stop("Argument 'x' is missing")
  if(!is.numeric(x)) stop("Argument 'x' must be a numerical vector")
  nmiss <- sum(is.na(x))
  if(nmiss>0) {
    x <- na.omit(x)  ## <----- gestione buchi
    warning(nmiss," missing values have been deleted")
    #
    #spline(x,xout=1:length(x),method="natural")$y
    }
  if(is.null(unit)) {
    n <- length(x)
    } else {
    n <- min(sapply(split(x,unit),length))
    }
  if(n<5) stop("At least 5 observations are required")
  if(is.null(max.lag)) {
    #max.lag <- min(n-3,trunc((n-1)^(1/3)))
    max.lag <- round(sqrt(n))
    } else {
    if(length(max.lag)>1) max.lag <- max.lag[1]
    if(!is.numeric(max.lag) || max.lag!=round(max.lag) || max.lag<0) stop("Argument 'max.lag' must be a non-negative integer value")
    if(max.lag>n-3) {
      max.lag <- n-3
      #warning("Argument 'max.lag' was set to ",n-3)
      }
    }
  if(is.null(unit)) {
    res <- adfFun(x=x, max.lag=max.lag)
    } else {
    gr <- unique(na.omit(unit))
    res <- vector("list",length=3)
    for(w in gr) {
      ind <- which(unit==w)
      iadf <- adfFun(x=x[ind], max.lag=max.lag)
      for(j in 1:length(iadf)) {
        res[[j]] <- c(res[[j]],iadf[[j]])
        }
      }
    res <- lapply(res, function(z){names(z)<-gr; z})
    names(res) <- names(iadf)
    m <- length(res$p.value)
    logp <- qnorm(res$p.value)
    rhat <- 1-var(logp)
    rstar <- max(rhat,-1/(m-1))
    auxz <- sum(logp)/sqrt(m*(1+(m-1)*(rstar+0.2*sqrt(2/(m+1))*(1-rstar))))
    #auxz <- sum(logp)/sqrt(m)
    res$p.value <- c(res$p.value,'(combined)'=2*pnorm(-abs(auxz)))
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
    res.sum <- summary(res)$coefficients
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
  if(max.lag>0) {
    k <- ar(x,order.max=max.lag)$order
    } else {
    k <- 0
    }
  res <- doADF(k)
  list(statistic=res[1],lag.selected=k,p.value=res[2])
  }

# function for kpss test (internal use only)
kpssFun <- function (x,trend,lshort) {
  n <- length(x)
  if(trend==T) {
    t <- 1:n
    e <- residuals.lm(lm(x ~ t))
    table <- c(0.216, 0.176, 0.146, 0.119)
    } else {
    e <- residuals.lm(lm(x ~ 1))
    table <- c(0.739, 0.574, 0.463, 0.347)
    }
  tablep <- c(0.01, 0.025, 0.05, 0.1)
  s <- cumsum(e)
  eta <- sum(s^2)/(n^2)
  s2 <- sum(e^2)/n
  if(lshort==T) {
    l <- trunc(4*(n/100)^0.25)
    } else {
    l <- trunc(12*(n/100)^0.25)
    }
  k <- 0
  for(i in 1:l) {
    ik <- 0
    for(j in (i+1):n) {
      ik <- ik+e[j]*e[j-i]
      }
    k <- k+(1-i/(l+1))*ik
    }
  STAT <- eta/(s2+2*k/n)
  PVAL <- approx(table,tablep,STAT,rule=2)$y
  list(statistic=STAT,lag.order=l,p.value=PVAL)
  }

# apply differencing
tsDiff <- function(var.names=NULL, unit=NULL, time=NULL, data, box.cox=1, ndiff=0) {
  if(!identical(class(data),"data.frame")) stop("Argument 'data' must be a data.frame")
  #
  is.dummy <- function(x) {
    dom <- sort(unique(na.omit(x)))
    length(dom)==2&dom[1]==0&dom[2]==1
    }
  if(is.null(var.names)) {
    var.names <- setdiff(colnames(data),time)
    if(length(var.names)==0) stop("No quantitative variable found")
    x2del <- c()
    for(i in 1:length(var.names)) {
      if(!is.numeric(data[,var.names[i]])|is.dummy(data[,var.names[i]])) x2del <- c(x2del,var.names[i])
      }
    if(length(x2del)>0) var.names <- setdiff(var.names,x2del)
    if(length(var.names)==0) stop("No quantitative variable found")
    } else {
    auxchk <- setdiff(var.names,colnames(data))  
    if(length(auxchk)>0) stop("Unknown variable '",auxchk[1],"' in argument 'var.names'")
    #var.names <- setdiff(intersect(var.names,colnames(data)),time)
    }
  #
  if(!is.null(time)) {
    if(length(time)>1) time <- time[1]
    if((time%in%colnames(data))==F) stop("Unknown variable '",time,"' in argument 'time'")
    if(time%in%var.names) stop("Variable '",time,"' appears in both arguments 'var.names' and 'time'")
    if(!is.numeric(data[,time])&!identical(class(data[,time]),"Date")) stop("Variable '",time,"' must be numeric or of class 'Date'")
    }
  #
  if(!is.vector(box.cox) || !is.numeric(box.cox)) stop("Argument 'box.cox' must be a numeric value or vector")
  if(length(box.cox)==1 & is.null(names(box.cox))) {
    box.cox <- rep(box.cox,length(var.names))
    names(box.cox) <- var.names
    } else {
    box.coxOK <- rep(1,length(var.names))
    names(box.coxOK) <- var.names
    box.coxOK[names(box.cox)] <- box.cox
    box.cox <- box.coxOK[var.names]
    }
  #
  if(!is.vector(ndiff) || !is.numeric(ndiff) || (sum(ndiff<0)>0 | sum(ndiff!=round(ndiff))>0)) stop("Argument 'ndiff' must be a non-negative integer value or vector")
  if(length(ndiff)==1 & is.null(names(ndiff))) {
    ndiff <- rep(ndiff,length(var.names))
    names(ndiff) <- var.names
    } else {
    ndiffOK <- rep(0,length(var.names))
    names(ndiffOK) <- var.names
    ndiffOK[names(ndiff)] <- ndiff
    ndiff <- ndiffOK[var.names]
    } 
  #
  dataL <- data
  for(i in 1:length(var.names)) {
    ilam <- box.cox[var.names[i]]
    if(ilam==0 & sum(data[,var.names[i]]<=0)>0) {
      box.cox[var.names[i]] <- ilam <- 1
      warning("Logarithmic transformation not applied to variable '",var.names[i],"'",call.=F)
      }
    dataL[,var.names[i]] <- makeBoxCox(data[,var.names[i]],ilam)
    }
  #
  if(is.null(unit)) {
    dataD <- diffFun(var.names=var.names, time=time, data=dataL, ndiff=ndiff)
    attr(dataD,"inits") <- unlist(data[1,var.names])
    attr(dataD,"box.cox") <- box.cox
    attr(dataD,"ndiff") <- ndiff
    if(max(ndiff)>0) {
      dataD[setdiff(1:nrow(data),1:max(ndiff)),]
      } else {
      dataD
      }
    } else {
    if(length(unit)>1) unit <- unit[1]
    if((unit%in%colnames(data))==F) stop("Unknown variable '",unit,"' in argument 'unit'")
    if(unit%in%time) stop("Variable '",unit,"' appears in both arguments 'time' and 'unit'")
    if(unit%in%var.names) stop("Variable '",unit,"' appears in both arguments 'var.names' and 'unit'")
    dataD <- dataL
    isNA <- c()
    gr <- unique(na.omit(data[,unit]))
    val0 <- matrix(nrow=length(gr),ncol=length(var.names))
    rownames(val0) <- gr
    colnames(val0) <- var.names
    for(w in 1:length(gr)) {
      ind <- which(data[,unit]==gr[w])
      if(max(ndiff)>0) isNA <- c(isNA, ind[1]:ind[max(ndiff)])
      dataD[ind,] <- diffFun(var.names=var.names, time=time, data=dataL[ind,], ndiff=ndiff)
      val0[w,] <- as.numeric(data[ind[1],var.names])
      }
    attr(dataD,"inits") <- val0
    attr(dataD,"box.cox") <- box.cox
    attr(dataD,"ndiff") <- ndiff
    dataD[setdiff(1:nrow(data),isNA),]
    }
  }

# function for differencing (auxiliary)
diffFun <- function(var.names, time, data, ndiff) {
  newdat <- data
  if(!is.null(time)) newdat <- newdat[order(newdat[,time]),]
  n <- nrow(data)
  for(i in 1:length(var.names)) {
    idat <- data[,var.names[i]]
    if(ndiff[i]>0) {
      newdat[,var.names[i]] <- idat-c(rep(NA,ndiff[i]),idat[1:(n-ndiff[i])])
      } else {
      newdat[,var.names[i]] <- idat
      }
    }
  newdat
  }

# fit ols con gamma lag (auxiliary)
gam_olsFit <- function(y.name, x.names, z.names, unit, par, offset, data, normalize) {
  p <- length(x.names)
  if(normalize) normstr <- "" else normstr <- paste(",normalize=",normalize,sep="")
  if(is.null(unit)) {
    form0 <- paste(y.name," ~ ",sep="")
    gstr <- ""
    } else {
    form0 <- paste(y.name," ~ -1+",unit,"+",sep="")
    gstr <- paste(",unit=",unit,sep="")
    }
  for(i in 1:p) {
    if(i>1) form0 <- paste(form0,"+",sep="")
    if(offset[i]==0) offstr <- "" else offstr <- paste(",offset=",offset[i],sep="")
    form0 <- paste(form0,"gammaKernel(",x.names[i],",par=c(",paste(par[,i],collapse=","),")",gstr,offstr,normstr,")",sep="")
    }
  if(!is.null(z.names)) form0 <- paste(form0,"+",paste(z.names,collapse="+"),sep="")
  formOK <- formula(form0)
  mod <- lm(formOK,data=data)
  mod$call$formula <- formOK
  #mod$call$data <- deparse(substitute(data))
  #
  colnames(par) <- names(offset) <- x.names
  rownames(par) <- c("delta","lambda")
  mod$par <- par
  mod$offset <- offset
  #
  #rownames(par) <- c("delta","lambda")
  #parlist <- list()
  #for(i in 1:ncol(par)) parlist[[i]] <- c(par[,i],offset=offset[i])
  #names(parlist) <- x.names
  #mod$par <- parlist
  #
  mod$variables <- list(y.name=y.name,x.names=x.names,z.names=z.names)
  if(is.null(unit)) {
    idg <- NULL
    } else {
    idg <- lapply(split(data, data[,unit]), rownames)  
    }
  mod$unit.id <- idg
  mod$data <- data
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
gam_hcFun <- function(y.name, x.names, z.names, unit, data, offset, inits, visitList, gridList, sign, grid.by) {
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
          ijm <- gam_olsFit(y.name=y.name, x.names=x.names, z.names=z.names, unit=unit, par=ijpar, offset=offset, data=data, normalize=F)
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
  #parOK <- modOK$par
  modFinal <- gam_olsFit(y.name=y.name, x.names=x.names, z.names=z.names, unit=unit, par=parOK, offset=offset, data=data, normalize=T)
  list(model=modFinal,par.tested=testL)
  }

# format control options (auxiliary)
optFormat <- function(optList, nomi, val) {
  auxopt <- optList[nomi]
  if(is.null(auxopt)) auxopt <- vector("list",length=length(nomi))
  for(i in 1:length(nomi)) {
    iopt <- auxopt[[i]]
    if(!is.numeric(iopt)) iopt <- NULL 
    if(length(iopt)==1) {
      iopt <- rep(iopt,2)
      } else if(length(iopt)>2) {
      iopt <- iopt[1:2]
      } else if(length(iopt)<1) {
      iopt <- val
      }
    auxopt[[i]] <- iopt
    }
  names(auxopt) <- nomi
  auxopt
  }

# MASTER FUNCTION
gammadlm <- function(y.name, x.names, z.names=NULL, unit=NULL, time=NULL, data, offset=rep(0,length(x.names)),
  control=list(nstart=50, grid.by=0.05, delta.lim=NULL, lambda.lim=NULL, peak.lim=NULL, length.lim=NULL),
  quiet=FALSE) {
  #
  if(!identical(class(data),"data.frame")) stop("Argument 'data' must be a data.frame")
  #
  if(missing(y.name)) stop("Argument 'y.name' is missing")
  if(length(y.name)>1) y.name <- y.name[1]
  if((y.name%in%colnames(data))==F) stop("Unknown variable '",y.name,"' in argument 'y.name'")
  if(sum(is.na(data[,y.name]))>0) stop("Variable '",y.name,"' contains missing values")
  #
  if(missing(x.names)) stop("Argument 'x.names' is missing")
  auxchk1 <- setdiff(x.names,colnames(data))
  if(length(auxchk1)>0) stop("Unknown variable '",auxchk1[1],"' in argument 'x.names'")
  for(i in 1:length(x.names)) {
    if(sum(is.na(data[,x.names[i]]))>0) stop("Variable '",x.names[i],"' contains missing values")
    }
  if(!is.null(z.names)) {
    auxchk2 <- setdiff(z.names,colnames(data))
    if(length(auxchk2)>0) stop("Unknown variable '",auxchk2[1],"' in argument 'z.names'")
    for(i in 1:length(z.names)) {
      if(sum(is.na(data[,z.names[i]]))>0) stop("Variable '",z.names[i],"' contains missing values")
      }
    if(y.name%in%z.names) stop("Variable '",y.name,"' appears in both arguments 'y.name' and 'z.names'")
    auxchk3 <- intersect(z.names,x.names)
    if(length(auxchk3)>0) stop("Variable '",auxchk3[1],"' appears in both arguments 'x.names' and 'z.names'")
    }
  #
  if(!is.null(time)) {
    if(length(time)>1) time <- time[1]
    if((time%in%colnames(data))==F) stop("Unknown variable '",time,"' in argument 'time'")
    if(time%in%y.name) stop("Variable '",time,"' appears in both arguments 'y.name' and 'time'")
    if(time%in%x.names) stop("Variable '",time,"' appears in both arguments 'x.names' and 'time'")
    if(!is.null(z.names) && time%in%z.names) stop("Variable '",time,"' appears in both arguments 'z.names' and 'time'")
    if(!is.numeric(data[,time])&!identical(class(data[,time]),"Date")) stop("Variable '",time,"' must be numeric or of class 'Date'")
    }
  #
  if(is.null(unit)) {
    if(!is.null(time)) data <- data[order(data[,time]),]
    } else {
    if(length(unit)>1) unit <- unit[1]
    if((unit%in%colnames(data))==F) stop("Unknown variable '",unit,"' in argument 'unit'")
    if(unit%in%y.name) stop("Variable '",unit,"' appears in both arguments 'y.name' and 'unit'")
    if(unit%in%x.names) stop("Variable '",unit,"' appears in both arguments 'x.names' and 'unit'")
    if(unit%in%z.names) stop("Variable '",unit,"' appears in both arguments 'z.names' and 'unit'")
    if(unit%in%time) stop("Variable '",unit,"' appears in both arguments 'time' and 'unit'")
    if(!is.null(time)) {
      gr <- unique(na.omit(unit))
      for(w in gr) {
        ind <- which(data[,unit]==w)
        idat <- data[ind,]
        data[ind,] <- idat[order(idat[,time]),]
        }
      }
    }
  #
  nstart <- control$nstart
  if(is.null(nstart)) nstart <- 50
  if(length(nstart)>1) nstart <- nstart[1]
  if(!is.numeric(nstart) || nstart!=round(nstart) || nstart<=0) stop("Argument 'nstart' must be a positive integer value")
  #max.try <- control$max.try
  #if(is.null(max.try)) max.try <- 50
  #max.start <- control$max.start
  #if(is.null(max.start)) max.start <- max.try*10
  grid.by <- control$grid.by
  if(is.null(grid.by)) grid.by <- 0.05
  if(length(grid.by)>1) grid.by <- grid.by[1]
  if(!is.numeric(grid.by) || grid.by<=0 || grid.by>0.1) stop("Argument 'grid.by' must be a positive value no greater than 0.1")
  delta.lim <- control$delta.lim
  lambda.lim <- control$lambda.lim
  peak.lim <- control$peak.lim
  length.lim <- control$length.lim
  sign <- NULL  ## <----- constraint segno
  if(!is.null(nstart)) max.try <- max.start <- nstart
  #
  if(!is.numeric(offset)) stop("Argument 'offset' must be a numerical vector")
  if(length(offset)==1) offset <- rep(offset,length(x.names))
  if(length(offset)>length(x.names)) offset <- offset[1:length(x.names)]
  if(length(offset)<length(x.names)) offset <- c(offset,rep(0,length(x.names)-length(offset)))    
  if(y.name%in%x.names && offset[which(x.names==y.name)]<1) offset[which(x.names==y.name)] <- 1
  #
  delta.lim <- optFormat(delta.lim, x.names, c(0,1))
  lambda.lim <- optFormat(lambda.lim, x.names, c(0,1))
  if(sum(sapply(delta.lim,function(x){x[1]==x[2]})==F)==0 &
      sum(sapply(lambda.lim,function(x){x[1]==x[2]})==F)==0) {
    par <- rbind(sapply(delta.lim,function(x){x[1]}),sapply(lambda.lim,function(x){x[1]}))
    } else {
    par <- NULL
    }
  if(!is.null(par)) {
    modOK <- gam_olsFit(y.name=y.name, x.names=x.names, par=par, z.names=z.names, unit=unit, offset=offset, data=data, normalize=T)
    } else {
    p <- length(x.names)
    peak.lim <- optFormat(peak.lim, x.names, c(0,Inf))
    length.lim <- optFormat(length.lim, x.names, c(0,Inf))
    #if(is.null(sign)) sign <- rep(0,p)
    #
    if(quiet==F & nstart>1) {
      cat("Scanning valid models ...")
      flush.console()
      }
    gridList <- gam_parGrid(delta.lim=delta.lim, lambda.lim=lambda.lim, peak.lim=peak.lim, length.lim=length.lim, grid.by=grid.by)
    pcombtot <- prod(sapply(gridList,nrow))
    if(quiet==F & nstart>1) cat("\r","Found ",signif(pcombtot)," valid models     ",sep="","\n")
    if(is.null(gridList)) {
      if(quiet==F & nstart>1) {
        cat("\n")
        cat("No solution with the current constraints","\n")
        }
      } else {
      rssOK <- Inf
      ntry_fit <- stopped <- 0
      modOK <- NULL
      searchList <- visitList <- list()
      #
      if(nstart==1) {
        ini0 <- matrix(0,nrow=2,ncol=length(x.names))
        gs0 <- gam_hcFun(y.name=y.name, x.names=x.names, z.names=z.names, unit=unit, data=data,
                         offset=offset, inits=ini0, visitList=visitList, gridList=gridList,
                         sign=sign, grid.by=grid.by)
        modOK <- gs0$model
        #modOK$local.max <- NULL
        }
      else
      #
      for(i in 1:max.start) {
        #if(quiet==F) {
        #  cat('\r',"Restart ",i,"/",nstart,": explored ",length(visitList)," valid models",sep="")
        #  flush.console()
        #  }
        ini0 <- gam_inits(gridList=gridList, visitList=visitList, maxtry=max.start)
        if(!is.null(ini0)) {
          gs0 <- gam_hcFun(y.name=y.name, x.names=x.names, z.names=z.names, unit=unit, data=data,
                           offset=offset, inits=ini0, visitList=visitList, gridList=gridList,
                           sign=sign, grid.by=grid.by)
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
      if(quiet==F & nstart>1) {
        cat("\n")
        if(stopped==2) cat("No more valid models found. End",sep="","\n")
        #if(stopped==1) {
        #  cat("No improvement found in the last ",max.try, " restarts. End",sep="","\n")
        #  } else if(stopped==2) {
        #    cat("No more valid initial values found. End",sep="","\n")
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
  #
  #if(!is.null(modOK)) {  ## <----- fitted values aggiustati per l'autocorrelazione
  #  resid <- modOK$residuals
  #  arRes <- ar(resid)
  #  if(arRes$order>0) {
  #    epsFit <- unconsKernel(resid, nlag=arRes$order, imputation=T)%*%c(0,arRes$ar)
  #    modOK$adj.fitted.values <- modOK$fitted.values-c(epsFit)
  #    } else {
  #    modOK$adj.fitted.values <- modOK$fitted.values
  #    }
  #  }
  #
  #
  if(is.null(unit)) {
    pval <- adfOneTest(modOK$residuals)$p.value
    } else {
    pval <- adfOneTest(modOK$residuals, unit=data[,unit])$p.value["(combined)"]
    }
  if(pval>0.05) warning("ADF test on residuals is not significant: regression could be spurious", call.=F)
  modOK
  }

# summary method for class 'gammadlm'
summary.gammadlm <- function(object, ...) {
  summ <- summary.lm(object, ...)
  ttab <- summ$coefficients
  S <- hacCalc(Xmat=model.matrix(object), resid=object$residuals, uS=summary.lm(object)$cov.unscaled, unitID=object$unit.id)
  bse <- sqrt(diag(S))
  ttab[,2] <- bse
  ttab[,3] <- ttab[,1]/ttab[,2]
  ttab[,4] <- round(2*pt(-abs(ttab[,3]),object$df.residual),6)
  summ$coefficients <- ttab
  summ
  }

# vcov method for class 'gammadlm'
vcov.gammadlm <- function(object, ...) {
  hacCalc(Xmat=model.matrix(object), resid=object$residuals, uS=summary.lm(object)$cov.unscaled, unitID=object$unit.id)
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

# AR test for residuals (auxiliary)
arTest <- function(resid, max.order=NULL) {
  n <- length(resid)
  if(is.null(max.order)) max.order <- trunc((length(resid)-1)^(1/3))
  ind <- (max.order+1):n
  y <- resid[ind]
  p.current <- max.order
  bicOK <- Inf
  fine <- 0
  while(fine==0) {
    if(p.current>0) {
      X <- unconsKernel(resid, nlag=p.current, imputation=F)[ind,-1,drop=F]
      mod <- lm(y~-1+X)
      } else {
      mod <- lm(y~-1)
      }
    #ichisq <- (n-i)*summary(im)$r.squared
    #ipval <- 1-pchisq(ichisq,i)
    bic <- extractAIC(mod,k=log(n))[2]
    if(bic<bicOK) {
      pOK <- p.current
      bicOK <- bic
      p.current <- p.current-1
      if(p.current<0) fine <- 1
      } else {
      fine <- 1
      }
    }
  if(pOK>0) {
    mOK <- lm(resid~unconsKernel(resid, nlag=pOK, imputation=F)[,-1,drop=F])
    bhat <- mOK$coef
    names(bhat) <- 0:(length(bhat)-1)
    } else {
    mOK <- lm(resid~-1)
    bhat <- c()
    }
  list(order=pOK,ar=bhat,var=summary(mOK)$sigma^2)
  }

# white test (auxiliary)
whitest <- function(Xmat, resid, max.degree) {
  n <- length(resid)
  p <- ncol(Xmat)-1
  res2 <- resid^2
  if(max.degree>1) {
    auxmat <- Xmat[,-1]
    for(i in 2:max.degree) {
      if((1+p*i)/n<=2/3) Xmat <- cbind(Xmat,auxmat^i)
      }
    }
  m0 <- lm(res2~Xmat)
  fstat <- summary(m0)$fstatistic
  pval <- pf(fstat[1],fstat[2],fstat[3])
  c(fstat,p.value=1-unname(pval))
  }

# hac covariance matrix (auxiliary)
hacCalc <- function(Xmat, resid, uS=NULL, unitID=NULL) {
  #
  #max.degree <- 3  ## <----- max degree for white test
  #wtest <- whitest(Xmat=Xmat, resid=resid, max.degree=max.degree)
  #homo <- wtest["p.value"]>0.05
  #
  homo <- ar(resid^2)$order==0
  #
  if(is.null(uS)) uS <- solve(t(Xmat)%*%Xmat)
  if(is.null(unitID)) {
    max.lag <- ar(resid)$order
    } else {
    max.lag <- c()
    for(i in 1:length(unitID)) {
      max.lag[i] <- ar(resid[unitID[[i]]])$order
      }
    }
  if(sum(max.lag)==0 & homo==T) {
    uS*sum(resid^2)/(length(resid)-ncol(Xmat))
    } else {
    #
    #Wcalc <- function(k) {
    #  W <- matrix(0,nrow=p,ncol=p)
    #  for(i in (k+1):n) {
    #    W <- W+resid[i]*resid[i-k]*(Xmat[i,]%*%t(Xmat[i-k,]))
    #    }
    #  W
    #  }
    #HH <- Wcalc(0)
    #if(max.lag>0) {
    #  for(k in 1:max.lag) {
    #    #iHH <- matrix(0,nrow=p,ncol=p)
    #    #for(i in (k+1):n) {
    #    #  iHH <- iHH+res[i]*res[i-k]*(Xmat[i,]%*%t(Xmat[i-k,])+Xmat[i-k,]%*%t(Xmat[i,]))
    #    #  }
    #    HH <- HH+(1-k/(max.lag+1))*2*Wcalc(k)
    #    }
    #  }
    #n/(n-p)*S%*%HH%*%S
    #
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
    uS%*%t(Xmat)%*%W%*%Xmat%*%uS
    }
  }

# extract lag coefficients
lagCoef <- function(x, cumulative=FALSE, max.lag=NULL, max.quantile=0.999) {
  if(missing(x)) stop("Argument 'x' is missing")
  if(!identical(class(x),c("gammadlm","lm"))) stop("Argument 'x' must be an object of class 'gammadlm'")
  #
  if(length(cumulative)>1) cumulative <- cumulative[1]
  if(!is.logical(cumulative)) stop("Argument 'cumulative' must be a logical value")
  #
  if(!is.null(max.lag)) {
    if(length(max.lag)>1) max.lag <- max.lag[1]
    if(!is.numeric(max.lag) || max.lag!=round(max.lag) || max.lag<0) stop("Argument 'max.lag' must be a non-negative integer value")
    }
  #
  if(length(max.quantile)>1) max.quantile <- max.quantile[1]
  if(!is.numeric(max.quantile) || (max.quantile<=0 | max.quantile>=1)) stop("Argument 'max.quantile' must be a value >0 and <1")
  #
  gpar <- x$par
  offs <- x$offset
  p <- ncol(gpar)
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
    ires[,1] <- x$coef[i+1]*lagwei
    ires[,2] <- sqrt(Smat[i+1,i+1])*lagwei
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
  res
  }

# plot method for class 'gammadlm'
plot.gammadlm <- function(x, x.names=NULL, conf=0.95, max.lag=NULL, max.quantile=0.999,
  xlim=NULL, ylim=NULL, add.legend=TRUE, cex.legend=1, digits=4, grid.length=100, main=NULL, ylab=NULL, xlab=NULL, ...) {
  #
  if(length(conf)>1) conf <- conf[1]
  if(!is.numeric(conf) || (conf<=0 | conf>=1)) stop("Argument 'conf' must be a value >0 and <1")
  #
  if(!is.null(max.lag)) {
    if(length(max.lag)>1) max.lag <- max.lag[1]
    if(!is.numeric(max.lag) || max.lag!=round(max.lag) || max.lag<0) stop("Argument 'max.lag' must be a non-negative integer value")
    }
  #
  if(length(max.quantile)>1) max.quantile <- max.quantile[1]
  if(!is.numeric(max.quantile) || (max.quantile<=0 | max.quantile>=1)) stop("Argument 'max.quantile' must be a value >0 and <1")
  #
  if(length(add.legend)>1) add.legend <- add.legend[1]
  if(!is.logical(add.legend)) stop("Argument 'add.legend' must be a logical value")
  #
  if(length(grid.length)>1) grid.length <- grid.length[1]
  if(!is.numeric(grid.length) || (grid.length!=round(grid.length) | grid.length<100)) stop("Argument 'grid.length' must be a value >=100")
  #
  makePlot <- function(i,main) {
    gpar <- x$par[,i]
    offs <- x$offset[i]
    if(is.null(max.lag)) {
      laglen <- gammaQuantile(prob=max.quantile, par=gpar, offset=offs)
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
    bcoef <- x$coef[i+1]*lagwei
    bcoef_se <- sqrt(Smat[i+1,i+1])*lagwei
    auxinf <- which(abs(bcoef_se)==Inf|is.na(abs(bcoef_se)))
    if(length(auxinf)>0) {
      auxbse <- bcoef_se
      auxbse[auxinf] <- NA
      #auxpred <- approx(auxbse,n=length(auxbse))$y
      auxpred <- spline(1:length(auxbse),auxbse,xout=1:length(auxbse))$y
      bcoef_se[auxinf] <- auxpred[auxinf]
      }
    #
    lagco <- cbind(bcoef,bcoef-tquan*bcoef_se,bcoef+tquan*bcoef_se)
    if(is.null(ylim)) ylim <- range(lagco)
    plot(lseq, lagco[,1], type="n", ylim=ylim, ylab=ylab, xlab=xlab, main=main, ...)
    polygon(c(lseq,rev(lseq)), c(lagco[,2],rev(lagco[,3])),
            border=NA, col=adjustcolor("grey75",alpha.f=0.5))
    grid()
    #abline(h=0,lty=3,col=1)
    lines(lseq,lagco[,1])
    if(add.legend) {
      auxw <- gammaWeights(0:(xlim[2]-1), par=gpar, offset=offs, normalize=T)
      bcum <- x$coef[i+1]*sum(auxw)
      auxse <- sqrt(Smat[i+1,i+1])*auxw
      bcum_se <- sqrt(sum(auxse%*%t(auxse)))
      bcumco <- cbind(bcum,bcum-tquan*bcum_se,bcum+tquan*bcum_se)
      if(is.null(cex.legend)) cex.legend <- 1
      bcumcoOK <- signif(bcumco)
      #siglags <- which(abs(bcoef/bcoef_se)>tquan)-1
      legtxt <- paste(#"Significant lags: ",min(siglags)," to ",max(siglags),"\n",
                      #", peak at ",ifelse(gpar[2]>0,round(gpar[1]/(gpar[1]-1)/log(gpar[2])-1,1),0)+offs,"\n",
                      "Cumulative coefficient: ",round(bcumcoOK[1],digits),"\n",
                      "   ",100*conf,"% CI: (",round(bcumcoOK[2],digits),", ",round(bcumcoOK[3],digits),")",sep="")
      legend("topright",legend=legtxt,cex=cex.legend,bty="n")
      }
    box()
    }
  xnam <- x$variables$x.names
  if(!is.null(x.names)) {
    chk0 <- setdiff(x.names,xnam)
    if(length(chk0)>0) warning(paste("Unknown variables: ",paste(chk0,collapse=", "),sep=""))
    xOK <- setdiff(x.names,chk0)
    if(length(xOK)==0) xOK <- xnam
    } else {
    xOK <- xnam
    }
  if(is.null(ylab)) ylab <- "Coefficient"
  if(is.null(xlab)) xlab <- "Time lag"
  Smat <- vcov(x)
  if(is.null(main)) main <- xOK
  #tquan <- qt((1+conf)/2,x$df.residual)
  tquan <- qnorm((1+conf)/2)
  mfrow0 <- par()$mfrow
  par(mfrow=n2mfrow(length(xOK)))
  for(i in 1:length(xOK)) {
    ind <- which(xnam==xOK[i])
    if(length(main)>=i) imain <- main[i] else imain <- ""
    makePlot(ind,imain)
    }
  par(mfrow=mfrow0)
  }
