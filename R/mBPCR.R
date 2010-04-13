.packageName <- "mBPCR"

logAdd <- function(x)
 {
    if (NCOL(x) == 1) {
        y <- max(x)
    if (y == -Inf) return(y)
    else return(y+log(sum(exp(x-y))))
    } else {
    y <- array(dim=NCOL(x))
    for (i in 1:NCOL(x)) {
        y[i] <- max(x[, i])
    }
    Y <- t(matrix(y, nrow=NCOL(x), ncol=NROW(x)))
    return(y+log(colSums(exp(x-Y))))
    }
 }

indexLA0 <- function(r, c, n)
 {
    if (length(r) == 1) {
	  if (length(c) == 1) {
            index <- c+(r-1)*(n-r/2)
	  } else {
            c1 <- c[1]:c[2]
            c1 <- c1[c1 >= r]
            index <- c1+(r-1)*(n-r/2)
        }
        return(index)
    }else{
        r1 <- r[1]:r[2]
        r1 <- r1[r1 <= c]
        index <- c+(r1-1)*(n-r1/2)
        return(index)
    }
 }

estGlobParam <- function(y, nu=NULL, rhoSquare=NULL, sigmaSquare=NULL, typeEstRho=1)
 {
    message('Estimation of global  parameters')
    n <- length(y)
    y <- c(y, y[1]) 
    m <- sum(y)
    s <- sum(y^2)
    l <- sum((y[1:(n)]-y[2:(n+1)])^2) 
    if (length(nu) == 0) nu <- m/n
    if (length(sigmaSquare) == 0) sigmaSquare <- l/(2*(n))
    if (length(rhoSquare) == 0) {
        if (typeEstRho == 1) {
            ##computation of rho^2 hat 1 estimator
            rhoSquare <- abs(sum((y[1:(n)]-m/n)*(y[2:(n+1)]-m/n)))/(n)
        } else { 
            if (typeEstRho == 0) {
                ##computation of rho^2 hat estimator
                rhoSquare <- s/n-(m/n)^2
            } else {
                stop('wrong value for the parameter typeEstRho')
            }
        }
    }    
    return(list(nu=nu, rhoSquare=rhoSquare, sigmaSquare=sigmaSquare))
 }

computeLA0Vect <- function(y, nu, rhoSquare, sigmaSquare)
 {
    n<-length(y)    
    message('Computation of log(A^0)')
    lA0 <- array(dim=sum(1:(n+1)))
    for (i in 1:(n)) {
        dist <- (seq(1, n+1)-i)[(i+1):(n+1)]
        ysum <- cumsum(y[i:n]-nu)
        y2sum <- cumsum((y[i:n]-nu)^2)
    lA0[indexLA0(i, c(i+1, n+1), n+1)] <- -0.5*dist*log(2*pi*sigmaSquare)-0.5*log(1+dist*rhoSquare/sigmaSquare)+0.5/sigmaSquare*(ysum^2/(dist+sigmaSquare/rhoSquare)-y2sum)
    lA0[indexLA0(i, i, n+1)] <- -Inf   
    }
    lA0[length(lA0)] <- -Inf
    return(lA0)
 } 

computeA10 <- function(i, j, y, nu, rhoSquare, sigmaSquare)
 {
    if (length(j) > 1) {
	  if (length(j) == 2) {
            ysum <- cumsum(y[(i+1):j[2]]-nu)
            dist <- (j[1]-i):(j[2]-i)
            a <- (rhoSquare*(ysum[(j[1]-i):(j[2]-i)]+dist*nu)+sigmaSquare*nu)/(dist*rhoSquare+sigmaSquare)
        } else {
            stop('wrong value for parameter j')
            return(NULL)
        }
    } else { 
        a <- (rhoSquare*sum(y[(i+1):j])+sigmaSquare*nu)/((j-i)*rhoSquare+sigmaSquare)
    }
    return(a)
 }

computeRecursions <- function(lA0, n, kMax=50)
 {  
    message('Computation of left and right recursions')
    lL <- matrix(-Inf, nrow=kMax+1, ncol=n+1)
    lR <- matrix(-Inf, nrow=kMax+1, ncol=n+1)
    lL[1, 1] <- 0
    lR[1, n+1] <- 0
    for ( k in 1:kMax) {                 
        for (j in 1:(n+1)) {                
        if ((j-1) >= k) {
                lL[k+1, j] <- logAdd(lL[k, k:(j-1)]+lA0[indexLA0(c(k, (j-1)), j, n+1)])
            } 
        if((j+1) <= (n+2-k)) {
            lR[k+1, j] <- logAdd(lA0[indexLA0(j, c((j+1), (n+2-k)), n+1)]+lR[k, (j+1):(n+2-k)])
        }
        }
    }
    return(list(lL=lL, lR=lR))  
 } 

computeRegrCurve <- function(y, typeRegr="BRC", n, kMax=50, lL, lR, lA0, nu, rhoSquare, sigmaSquare, option)
 {  
    if (typeRegr == "BRC") {
        message('Computation of Bayesian regression Curve')
        kml <- option
        f=lA0-lL[kml+1, n+1]
        for(i in 1:(n)){
            a <- c(0, computeA10(i-1, c(i, n), y, nu, rhoSquare, sigmaSquare))
            if (kml > 1) { 
                f[indexLA0(i, c(i, n+1), n+1)] <- a*exp(f[indexLA0(i, c(i,n+1), n+1)]+as.vector(log(exp(lL[1:kml, i])%*%exp(lR[kml-seq(1, kml)+1, (i):(n+1)]))))            
            }else{
                f[indexLA0(i, c(i, n+1), n+1)] <- a*exp(f[indexLA0(i, c(i,n+1), n+1)]+as.vector(log(exp(lL[kml, i])*exp(lR[1, (i):(n+1)]))))
            } 
        }
        regrEst <- array(dim=n+1)
        regrEst[1] <- 0  
        for (h in 1:n) {   
            regrEst[h+1] <- regrEst[h]+sum(f[indexLA0(h, c(h+1, n+1), n+1)])-sum(f[indexLA0(c(1, h), h, n+1)])
        }
        remove(f)
    } else {
        message('Computation of Bayesian regression Curve Ak')
        lC <- option
        ff <- matrix(nrow=kMax, ncol=n+1)
        for (k in 1:kMax) {
        f <- lA0-lL[k+1, n+1]
        for (i in 1:(n)) {
            a <- c(0, computeA10(i-1, c(i, n), y, nu, rhoSquare, sigmaSquare))
        if (k > 1) { 
            f[indexLA0(i, c(i, n+1), n+1)] <- a*exp(f[indexLA0(i, c(i, n+1), n+1)]+as.vector(log(exp(lL[1:k, i])%*%exp(lR[k-seq(1, k)+1, (i):(n+1)]))))         
        }else{
            f[indexLA0(i, c(i, n+1), n+1)] <- a*exp(f[indexLA0(i,c(i, n+1), n+1)]+as.vector(log(exp(lL[k, i])*exp(lR[1, (i):(n+1)]))))
        } 
        }
        S1 <- array(dim=n+1)
        S1[1] <- 0 
        for (h in 1:n) {
            S1[h+1] <- S1[h]+sum(f[indexLA0(h, c(h+1, n+1), n+1)])-sum(f[indexLA0(c(1, h), h, n+1)])
        }
        ff[k, ] <- S1
        remove(f)
        }
        regrEst <- t(exp(lC))%*%ff
        remove(ff)
    }
    return(regrEst[2:(n+1)])
 }  

computePCReg <- function(y, lA0, lL, lR, nu, rhoSquare, sigmaSquare, kMax=50, regr=NULL)
 {
    n=length(y)
    ##set a uniform distribution for k
    probK <- 1/kMax
    message('Determination of PC Regression');
    ##compute the logarithm of the posterior probability of k (lC[k]=log(p(k|y)))
    lC <- lL[seq(1, kMax)+1, n+1]-log(choose(n-1, seq(1, kMax)-1))+log(probK)
    lE <- logAdd(lC) 
    ##lE=logarithm of the evidence
    lC <- lC-lE
    ##estimator of k which minimizes the square error (k_2)
    ek <- sum(seq(1, kMax)*exp(lC))
    w <- which(lC > -Inf)
    err <- w^2-2*w*ek
    kml <- w[which.min(err)]
    ##estimator of the boundaries T_BinErrAk
    if (kml > 1) {  
        dd <- matrix(nrow=kMax-1, ncol=n-1)
        for (kk in 2:kMax) {
            for (i in 2:n) {
                dd[kk-1, i-1] <- logAdd(lL[2:kk, i]+lR[kk:2, i]-lL[kk+1, n+1])
            }
        }
        d1 <- array(dim=n-1)
        for (i in 1:(n-1)) {
            d1[i] <- exp(logAdd(dd[, i]+lC[2:kMax]))
        }
        s1 <- sort(d1)
        boundaries <- array(dim=kml-1)
        i <- 1
        while(i <= kml-1){
            if (length(which(d1 == s1[n-i])) == 1) {
                boundaries[i] <- which(d1 == s1[n-i])
                i <- i+1
            }else{
                ll <- length(which(d1 == s1[n-i]))
                boundaries[i:(i+ll-1)] <- which(d1 == s1[n-i])
                i <- i+ll
            } 
        }
        boundaries <- sort(boundaries)
        boundaries <- c(0, boundaries, n)
    }else{
        boundaries <- c(0, n)
        d1 <- array(0, dim=n-1)
    }
    ## Bayesian regression curve
    if (length(regr) != 0) {
        if (regr == "BRC")  {
            estRegr <- computeRegrCurve(y, regr, n, kMax, lL, lR, lA0, nu, rhoSquare, sigmaSquare, kml)
        }else{
            estRegr <- computeRegrCurve(y, regr, n, kMax, lL, lR, lA0, nu, rhoSquare, sigmaSquare, lC)
        } 
    }else{
        estRegr=NULL
    }
    ##estimation of the segment levels 
    mT <- array(dim=n)
    for (k in 1:kml) {
        m1 <- computeA10(boundaries[k], boundaries[k+1], y, nu, rhoSquare, sigmaSquare) 
        mT[(boundaries[k]+1):boundaries[k+1]] <- m1
    }     
    return(list(kml=kml, boundaries=boundaries, postProbT=d1, estPC=mT, estRegr=estRegr))  
 } 

computeMBPCR <- function(y, kMax=50, nu=NULL, rhoSquare=NULL, sigmaSquare=NULL, typeEstRho=1, regr=NULL)
 { 
    n <- length(y)
    kMax <- min(kMax, n)
    if(length(nu) == 0 || length(rhoSquare) == 0 || length(sigmaSquare) == 0) {
        results <- estGlobParam(y, nu, rhoSquare, sigmaSquare, typeEstRho)
        nu <- results$nu
        rhoSquare <- results$rhoSquare
        sigmaSquare <- results$sigmaSquare
    }
    lA0 <- computeLA0Vect(y, nu, rhoSquare, sigmaSquare)  
    results <- computeRecursions(lA0, n, kMax)
    lL <- results$lL
    lR <- results$lR
    remove(results)
    if (length(regr)!=0 && (regr != "BRC" & regr != "BRCAk")) {
        stop('wrong value for parameter regr')
        return(NULL)
    }
    results <- computePCReg(y, lA0, lL, lR, nu, rhoSquare, sigmaSquare, kMax, regr)
    return(list(estK=results$kml, estBoundaries=results$boundaries[-1], estPC=results$estPC, regrCurve=results$estRegr, nu=nu, rhoSquare=rhoSquare, sigmaSquare=sigmaSquare, postProbT=results$postProbT))
 }  

writeEstProfile <- function(path='', sampleName='', snpName, chr, position, logratio, chrToBeWritten, estPC, estBoundaries=NULL, postProbT=NULL, regrCurve=NULL, regr=NULL)
 {
    data11 <- data.frame('SNPname', 'chromosome', 'position', 'rawLog2ratio', 'mBPCRestimate')
    if (length(path) == 0) {
        nRowEst <- 0
        for (j in 1:length(chrToBeWritten)) {
            nRowEst <- nRowEst + length(which(chr == chrToBeWritten[j]))
        }
        mBPCRest <- matrix(nrow = nRowEst, ncol = 5)
        colnames(mBPCRest) <- as.matrix(data11)
        rowCountEst <- 0
    } else {
        PathResults1 <- paste(path, sampleName, '_mBPCRestimate.txt', sep='')
        write.table(data11, PathResults1, row.names=FALSE, col.names=FALSE, sep='\t')
    }
    if (length(estBoundaries) != 0) {
        if (length(postProbT) != 0) {
            data21 <- data.frame('SNPname(start)', 'SNPname(end)', 'chromosome', 'position(start)', 'position(end)', 'nProbes', 'mBPCRestimate', 'breakpointPostProb')
        } else {
            data21 <- data.frame('SNPname(start)', 'SNPname(end)', 'chromosome', 'position(start)', 'position(end)', 'nProbes', 'mBPCRestimate')
        }
        if (length(path) == 0) {
            nColEst <- length(data21)
            nRowEst <- length(unlist(estBoundaries))
            estBounds <- matrix(nrow = nRowEst, ncol = nColEst)
            colnames(estBounds) <- as.matrix(data21)
            rowCountBounds <- 0
        } else {
            PathResults2 <- paste(path, sampleName, '_mBPCRbreakpoints.txt', sep='')
            write.table(data21, PathResults2, row.names=FALSE, col.names=FALSE, sep='\t')
        }
    } else {
        if (length(postProbT) != 0) {
            stop('estBoundaries=NULL while posteriorProbT!=NULL')
            return(NULL)
        }
    }
    if(length(regrCurve) != 0 && length(which(!is.na(regrCurve))) != 0) {
        if (length(regr) == 0 || (regr != "BRC" & regr != "BRCAk") ) {
            stop('wrong value for parameter regr')
            return(NULL)
        } else {
            if (regr == "BRC") {
                data31 <- data.frame('SNPname', 'chromosome', 'position', 'rawLog2ratio', 'mBRC_estimate')
                if (length(path) > 0) {
                    PathResults3 <- paste(path, sampleName, '_mBRCestimate.txt', sep='')
                    write.table(data31, PathResults3, row.names=FALSE, col.names=FALSE, sep='\t')
                }
            } else {
                data31 <- data.frame('SNPname', 'chromosome', 'position', 'rawLog2ratio', 'BRCAk_estimate')
                if (length(path) > 0) {           
                    PathResults3 <- paste(path, sampleName, '_BRCAkestimate.txt', sep='')
                    write.table(data31, PathResults3, row.names=FALSE, col.names=FALSE, sep='\t')
                }
            }
            if (length(path) == 0) {
                mBRCest <- matrix(nrow = dim(mBPCRest)[1],ncol = dim(mBPCRest)[2])
                colnames(mBRCest) <- as.matrix(data31)
            }
        }
    } else {
        regr <- NULL
    }
    for (j in 1:length(chrToBeWritten)) {
        data12 <- data.frame(snpName[chr == chrToBeWritten[j]], chr[chr == chrToBeWritten[j]], position[chr == chrToBeWritten[j]], logratio[chr == chrToBeWritten[j]], estPC[chr == chrToBeWritten[j]])    
        if (length(path) == 0) {
            mBPCRest[(rowCountEst + 1):(rowCountEst + dim(data12)[1]),] <- as.matrix(data12)
            if (length(regr) == 0) rowCountEst <- rowCountEst + dim(data12)[1]
        } else {
            write.table(data12, PathResults1, row.names=FALSE, col.names=FALSE, sep='\t', append=TRUE)    
        }
        if (length(estBoundaries) != 0) {
            startBoundaries <- estBoundaries[[j]]+1
            if (length(which(!is.na(estPC[chr == chrToBeWritten[j]]))) != length(which(chr == chrToBeWritten[j]))) {
                if (estBoundaries[[j]][length(estBoundaries[[j]])] == length(which(chr == chrToBeWritten[j]))) {
                    startBoundaries <- c(1, length(which(is.na(estPC[chr == chrToBeWritten[j]])))+1, startBoundaries[-length(estBoundaries[[j]])])
                    estBoundaries[[j]] <- c(length(which(is.na(estPC[chr == chrToBeWritten[j]]))), estBoundaries[[j]])
                    if (length(postProbT) != 0) postProbT[[j]] <- c(NA, postProbT[[j]])
                } else {
                    startBoundaries <- c(1, startBoundaries[-length(estBoundaries[[j]])], length(which(is.na(estPC[chr == chrToBeWritten[j]])))+1)
                    estBoundaries[[j]] <- c(estBoundaries[[j]], length(which(chr == chrToBeWritten[j])))
                    if (length(postProbT) != 0) postProbT[[j]] <- c(postProbT[[j]], NA)
                } 
            } else {
                startBoundaries <- c(1, startBoundaries[-length(estBoundaries[[j]])])
            } 
            if (length(postProbT) != 0) {
                data22 <- data.frame(snpName[chr == chrToBeWritten[j]][startBoundaries], snpName[chr == chrToBeWritten[j]][estBoundaries[[j]]], array(chrToBeWritten[j], dim=length(estBoundaries[[j]])), position[chr == chrToBeWritten[j]][startBoundaries], position[chr == chrToBeWritten[j]][estBoundaries[[j]]], estBoundaries[[j]]-startBoundaries+1, estPC[chr == chrToBeWritten[j]][estBoundaries[[j]]], postProbT[[j]])
            } else {
                data22 <- data.frame(snpName[chr == chrToBeWritten[j]][startBoundaries], snpName[chr == chrToBeWritten[j]][estBoundaries[[j]]], array(chrToBeWritten[j], dim=length(estBoundaries[[j]])), position[chr == chrToBeWritten[j]][startBoundaries], position[chr == chrToBeWritten[j]][estBoundaries[[j]]], estBoundaries[[j]]-startBoundaries+1, estPC[chr == chrToBeWritten[j]][estBoundaries[[j]]])
            }
            if (length(path) == 0) {
                estBounds[(rowCountBounds + 1):(rowCountBounds + dim(data22)[1]),] <- as.matrix(data22)
                rowCountBounds <- rowCountBounds + dim(data22)[1]
            } else {
                write.table(data22, PathResults2, row.names=FALSE, col.names=FALSE, sep='\t', append=TRUE)
            }
        }
        if(length(regrCurve) != 0 && length(which(!is.na(regrCurve))) != 0) {
            data32 <- data.frame(snpName[chr == chrToBeWritten[j]], chr[chr == chrToBeWritten[j]], position[chr == chrToBeWritten[j]], logratio[chr == chrToBeWritten[j]], regrCurve[chr == chrToBeWritten[j]])
            if (length(path) == 0) {
                 mBRCest[(rowCountEst + 1):(rowCountEst + dim(data32)[1]),] <- as.matrix(data32)
                 rowCountEst <- rowCountEst + dim(data32)[1]
            } else {
                 write.table(data32, PathResults3, row.names=FALSE, col.names=FALSE, sep='\t', append=TRUE)
            }
        }
    }
    if (length(path) == 0) {
        if (length(regr) != 0 & length(estBoundaries) != 0) {
            estBounds <- data.frame(estBounds)
            colnames(estBounds)=as.matrix(data21)
            return(list(mBPCRestimate=data.frame(mBPCRest), mBPCRbreakpoints=estBounds , regrCurveEstimate=data.frame(mBRCest)))
        } else {
            if (length(regr) != 0 | length(estBoundaries) != 0) {
                if (length(regr) != 0) {
                    return(list(mBPCRestimate=data.frame(mBPCRest), regrCurveEstimate=data.frame(mBRCest)))
                } else {
                    estBounds <- data.frame(estBounds)
                    colnames(estBounds)=as.matrix(data21)
                    return(list(mBPCRestimate=data.frame(mBPCRest), mBPCRbreakpoints=estBounds))
                }
            } else {
                return(list(mBPCRestimate=data.frame(mBPCRest)))    
            }
        }
    }
 }

importCNData <- function(path, NRowSkip, ifLogRatio=1)
 {
    results <- read.table(file=path, as.is=TRUE, row.names=NULL, flush=TRUE, skip=NRowSkip)
    if ( ifLogRatio == 1) return(list(snpName=results[[1]], chr=results[[2]], position=results[[3]], logratio=results[[4]]))
    if ( ifLogRatio == 0) return(list(snpName=results[[1]], chr=results[[2]], position=results[[3]], logratio=log2(results[[4]])-1))
 }

estProfileWithMBPCR <- function(snpName, chr, position, logratio, chrToBeAnalyzed, maxProbeNumber, rhoSquare=NULL, kMax=50, nu=NULL, sigmaSquare=NULL, typeEstRho=1, regr=NULL)
 {
    if (length(nu) == 0 || length(rhoSquare) == 0 || length(sigmaSquare) == 0) {
        results <- estGlobParam(logratio, nu, rhoSquare, sigmaSquare, typeEstRho)
	  nu <- results$nu
        rhoSquare <- results$rhoSquare
        sigmaSquare <- results$sigmaSquare
    }
    indexNo <- NULL
    estPC <- array(dim=length(snpName))
    estBoundaries <- list(dim=length(chrToBeAnalyzed))
    postProbT <- list(dim=length(chrToBeAnalyzed))
    regrCurve <- array(dim=length(snpName))
    for (j in chrToBeAnalyzed) {
        y <- logratio[which(chr == j)]
        n <- length(y)
        if (n <= maxProbeNumber) {
            message(paste('Estimation of the profile of chromosome ', j, sep=''))
	      results <- computeMBPCR(y, kMax, nu, rhoSquare, sigmaSquare, typeEstRho, regr)
	      estPC[chr == j] <- results$estPC
            if (length(regr) != 0) regrCurve[chr == j] <- results$regrCurve
            estBoundaries[[which(chrToBeAnalyzed == j)]] <- results$estBoundaries
            postProbT[[which(chrToBeAnalyzed == j)]] <- c(results$postProbT[results$estBoundaries[-results$estK]],1)
	      remove(results)
        }else{
	  a <- mean(unlist(centromere(j)))		
          a <- which(position[chr == j] > a)[1]-1
          bounds1 <- NULL
          postProb1 <- NULL
          if (a > maxProbeNumber & length(which(chr == j))-a > maxProbeNumber){
              warning(paste('Warning: the profile of chromosome ',j,' has not been estimated because of its size', sep=''))
              indexNo <- c(indexNo, which(chrToBeAnalyzed == j))
          } else {
              message(paste('Estimation of the profile of chromosome ', j, sep=''))
          }
          if (a <= maxProbeNumber) { 
	        y <- y[1:a]
	        results <- computeMBPCR(y, kMax, nu, rhoSquare, sigmaSquare, typeEstRho, regr)
	        estPC[chr == j][1:a] <- results$estPC
		  if (length(regr) != 0) regrCurve[chr == j][1:a] <- results$regrCurve
              bounds1 <- c(bounds1, results$estBoundaries)
              postProb1 <- c(postProb1, c(results$postProbT[results$estBoundaries[-results$estK]],1))
              remove(results)
          }
          if (length(which(chr == j))-a <= maxProbeNumber) {
	        y <- logratio[which(chr == j)][(a+1):n]
	        results <- computeMBPCR(y, kMax, nu, rhoSquare, sigmaSquare, typeEstRho, regr)
	        estPC[chr == j][(a+1):n] <- results$estPC
              if (length(regr) != 0) regrCurve[chr == j][(a+1):n] <- results$regrCurve
              bounds1 <- c(bounds1, a+results$estBoundaries)
              postProb1 <- c(postProb1, c(results$postProbT[results$estBoundaries[-results$estK]],1))
              remove(results)
          }
          if (a <= maxProbeNumber | length(which(chr == j))-a <= maxProbeNumber){
              estBoundaries[[which(chrToBeAnalyzed == j)]] <- bounds1
              postProbT[[which(chrToBeAnalyzed == j)]] <- postProb1   
          } 
          if (a <= maxProbeNumber & length(which(chr == j))-a > maxProbeNumber){
              warning(paste('Warning: the profile of arm q of chromosome ',j,' has not been estimated because of its size', sep=''))
          } 
          if (a > maxProbeNumber & length(which(chr == j))-a <= maxProbeNumber){
              warning(paste('Warning: the profile of arm p of chromosome ',j,' has not been estimated because of its size', sep=''))
          }   
	}
    } 
    if (length(regr) != 0) {
        return(list(estPC=estPC, estBoundaries=estBoundaries, postProbT=postProbT, regrCurve=regrCurve))  
    } else {
        return(list(estPC=estPC, estBoundaries=estBoundaries, postProbT=postProbT))  
    }
 }
 
plotEstProfile <- function(sampleName='', chr, position, logratio, chrToBePlotted, estPC, maxProbeNumber, legendPosition='bottomleft', regrCurve=NULL, regr=NULL)
 {
    for (j in 1:length(chrToBePlotted)) {
        if(j>1) x11()
        if (length(estPC) != 0) { 
		plot(position[chr == chrToBePlotted[j]], logratio[chr == chrToBePlotted[j]], pch='.', cex=2, col='grey', xlab=paste('chromosome', chrToBePlotted[j], sep=' '), ylab='log2ratio')
        	title(paste(sampleName, sep=''))
            if (length(regr) == 0) {
                legend(x=legendPosition, legend=c('mBPCR'), lty=c(1), col=c(4))
            } else {
                if (regr == "BRC") legend(x=legendPosition, legend=c('mBPCR', 'BRC with K_2'), lty=c(1, 1), col=c(4, 2))
                if (regr == "BRCAk") legend(x=legendPosition, legend=c('mBPCR', 'BRCAk'), lty=c(1, 1), col=c(4, 2))
                if (regr != "BRC" & regr != "BRCAk"){
                    stop('wrong value for parameter regr')
                    return(NULL)
                }
            }
            n <- length(which(chr == chrToBePlotted[j]))
            if(n <= maxProbeNumber){
                points(position[chr == chrToBePlotted[j]], estPC[chr == chrToBePlotted[j]], type='l', col=4)
                if (length(regr) != 0 && (regr == "BRC" | regr == "BRCAk")) points(position[chr == chrToBePlotted[j]], regrCurve[chr == chrToBePlotted[j]], type='l', col='red')
            } else {
		a <- mean(unlist(centromere(chrToBePlotted[j]))) 		
                a <- which(position[chr == chrToBePlotted[j]] > a)[1]-1
                points(position[chr == chrToBePlotted[j]][1:a], estPC[chr == chrToBePlotted[j]][1:a], type='l', col=4)
                points(position[chr == chrToBePlotted[j]][(a+1):n], estPC[chr == chrToBePlotted[j]][(a+1):n], type='l', col=4)
                if (length(regr) != 0 && (regr == "BRC" | regr == "BRCAk")) points(position[chr == chrToBePlotted[j]][1:a], regrCurve[chr == chrToBePlotted[j]][1:a], type='l', col='red')
                if (length(regr) != 0 && (regr == "BRC" | regr == "BRCAk")) points(position[chr == chrToBePlotted[j]][(a+1):n], regrCurve[chr == chrToBePlotted[j]][(a+1):n], type='l', col='red')
            }
        } else {
		plot(position[chr == chrToBePlotted[j]], logratio[chr == chrToBePlotted[j]], pch='.', cex=2, col='grey', xlab=paste('chromosome', chrToBePlotted[j], sep=' '), ylab='log2ratio')
        	title(paste(sampleName, sep=''))
            if (regr == "BRC") legend(x=legendPosition, legend=c('BRC with K_2'), lty=c(1), col=c(4))
            if (regr == "BRCAk") legend(x=legendPosition, legend=c('BRCAk'), lty=c(1), col=c(4))
            if (regr != "BRC" & regr != "BRCAk"){
                stop('wrong value for parameter regr')
                return(NULL)
            }
            n <- length(which(chr == chrToBePlotted[j]))
            if(n <= maxProbeNumber){
                points(position[chr == chrToBePlotted[j]], regrCurve[chr == chrToBePlotted[j]], type='l', col=4)
            } else {
                a <- mean(unlist(centromere(chrToBePlotted[j])))	
                a <- which(position[chr == chrToBePlotted[j]] > a)[1]-1
                points(position[chr == chrToBePlotted[j]][1:a], regrCurve[chr == chrToBePlotted[j]][1:a], type='l', col=4)
                points(position[chr == chrToBePlotted[j]][(a+1):n], regrCurve[chr == chrToBePlotted[j]][(a+1):n], type='l', col=4)
            }
        }
    }
 }
estProfileWithMBPCRforOligoSnpSet <- function(sampleData, sampleToBeAnalyzed, chrToBeAnalyzed, maxProbeNumber, ifLogRatio=1, rhoSquare=NULL, kMax=50, nu=NULL, sigmaSquare=NULL, typeEstRho=1, regr=NULL)
 {
    if (class(sampleData) !="oligoSnpSet" || (length(assayData(sampleData)$copyNumber) == 0 | length(featureData(sampleData)$chromosome) == 0 | length(featureData(sampleData)$position) == 0)) {
        stop("wrong object-data in input")
    } else {
        log2ratios <- assayData(sampleData)$copyNumber
        snpName <- featureNames(featureData(sampleData))
        chr <- featureData(sampleData)$chromosome
        position <- featureData(sampleData)$position
        chr[chr == "X"] <- 23
        chr[chr == "Y"] <- 24
        chr <- as.numeric(chr)
        o <- order(chr)
        chr  <- chr[o]
        position  <- position[o]
        snpName  <- snpName[o]
        for (i in unique(chr)) {
            o <- order(position[chr == i])
            position[chr == i] <- position[chr == i][o]
            snpName[chr == i] <- snpName[chr == i][o]
        }
        log2ratios <- log2ratios[snpName,]
        if (ifLogRatio == 0) log2ratios <- log2(log2ratios) - 1
        chr[chr == 23] <- "X"
        chr[chr == 24] <- "Y"
        o <- order(snpName)
        cnMatrix <- matrix(nrow=length(chr), ncol=dim(log2ratios)[2]) 
        if (length(regr) > 0)  regrMatrix <- matrix(nrow=length(chr), ncol=dim(log2ratios)[2])  
        for (i in sampleToBeAnalyzed) {
            r <- estProfileWithMBPCR(snpName, chr, position, log2ratios[,i], chrToBeAnalyzed, maxProbeNumber, rhoSquare, kMax, nu, sigmaSquare, typeEstRho, regr)
            cnMatrix[,i] <- r$estPC[o]
            if (length(regr) > 0) regrMatrix[,i] <- r$regrCurve[o]
        } 
        resultPC <- new("oligoSnpSet", call=assayData(sampleData)$call, callProbability =assayData(sampleData)$callProbability, cnConfidence = matrix(NA, length(chr), dim(log2ratios)[2]), copyNumber=cnMatrix, annotation = annotation(sampleData) ,phenoData = phenoData(sampleData), featureData = featureData(sampleData))
	for(i in varLabels(featureData(resultPC))){
	    featureData(resultPC)[[i]] <- featureData(sampleData)[[i]]
        }
        if (length(regr) == 0) {
            return(list(estPC=resultPC))
        } else {
            resultRegr <- new("oligoSnpSet", call=assayData(sampleData)$call, callProbability =assayData(sampleData)$callProbability, cnConfidence = matrix(NA, length(chr), dim(log2ratios)[2]),  copyNumber=regrMatrix,annotation = annotation(sampleData) ,phenoData = phenoData(sampleData), featureData = featureData(sampleData))
	    for(i in varLabels(featureData(resultRegr))){
	        featureData(resultRegr)[[i]] <- featureData(sampleData)[[i]]
            }
            return(list(estPC=resultPC, regrCurve=resultRegr))
        }
    }
}
