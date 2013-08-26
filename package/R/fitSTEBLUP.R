#' fitSTEBLUP
#' 
#' @description Wrapper function for the estimation of the spatio-temporal-FH model as in Marhuenda (2012). Uses the
#' function \code{FitSpatioTemporalFH} provided by Sample Deliverable 22: Software on Small Area Estimation
#' 
#' @param formula formula for the fixed effects part
#' @param dat data with dependent and independent variable as they are used in formula. The data should 
#' contain a variable \code{Domain} and \code{Time} indicating the domain and time period for each observation respectively
#' @param beta, sigma, rho start values for the parameters
#' @param sigmaSamplingError True variances for the sampling errors
#' @param nDomains the number of domains. Will be determined by default, see \code{dat}
#' @param nTime the number of time periods. Will be determined by default, see \code{dat}
#' @param w0 proximity matrix - not needed in the call to \code{FitSpatioTemporalFH}
#' @param w at the moment row-standardized proximity matrix
#' @param tol numerical tolerance for algorithm
#' @param method method used by the function \code{\link{optim}} - not needed in the call to \code{FitSpatioTemporalFH}
#' @param maxIter maximal number of iterations
#'
#'@examples require(spatioTemporalData); require(SAE)
#' set.seed(3)
#' setup <- simSetupMarhuenda(nDomains=30, nTime=5, sarCorr=0.5, arCorr=0.5, n = 200)
#' output <- simRunMarhuenda(setup) 
#' dat <- slot(output[[1]], "data")[[1]]
#' result <- fitSTEBLUP(y~x, dat, c(0,1), c(1,1), c(0.5,0.5))
#' summary(result)
fitSTEBLUP <- function(formula, dat, beta, sigma, rho, 
                                sigmaSamplingError = unlist(lapply(1:nDomains, seSigmaClosure(nDomains, nTime), t = 1:nTime)),
                                nDomains = getNDomains(dat),
                                nTime = getNTime(dat),
                                w0 = w0Matrix(nDomains), 
                                w = wMatrix(nDomains),
                                tol = 1e-3, method = "Nelder-Mead", maxIter = 500) {
  
  modelSpecs <- prepareData(formula, dat, nDomains, nTime, beta, sigma, rho, sigmaSamplingError, w0, w, tol, method, maxIter)
  
  modelFit <- try(FitSpatioTemporalFH("B", modelSpecs$x, modelSpecs$y, modelSpecs$nDomains, modelSpecs$nTime,
                                  modelSpecs$sigmaSamplingError,
                                  theta0 = matrix(c(modelSpecs$sigma[1], modelSpecs$rho[1], modelSpecs$sigma[2], modelSpecs$rho[2])),
                                  W = modelSpecs$w, MAXITER = modelSpecs$maxIter, PRECISION = modelSpecs$tol,
                                  confidence = 0.95),
                  silent = TRUE)
  
  if (class(modelFit) == "try-error") {
      modelFit <- 
        list(estimates = numeric(length(nDomains*nTime)),
             beta = list("coef" = numeric(length(modelSpecs$beta))),
             theta <- data.frame("estimate" = numeric(length(c(modelSpecs$sigma, modelSpecs$rho)))))
    }
  
  output <- list(estimates = data.frame(yHat = modelFit$estimates),
                 beta = as.matrix(modelFit$beta["coef"]),
                 sigma = as.numeric(modelFit$theta[c(1, 3), "estimate"]),
                 rho = as.numeric(modelFit$theta[c(2, 4), "estimate"]))
  class(output) <- c("fitSTEBLUP", "fitSTREBLUP", "list")
  output
}

FitSpatioTemporalFH <- function(model,X,y,nD,nT,sigma2dt,theta0,
                                W,MAXITER,PRECISION,confidence)
{
  
  diagonalizematrix <- function(A,ntimes)
  {
    nrowA <- nrow(A)
    ncolA <- ncol(A)
    
    Adiag <- matrix(0,nrow=nrowA*ntimes, ncol=ncolA*ntimes)
    
    firsti <- 1
    firstj <- 1
    for (n in 1:ntimes)
    {
      lasti <- firsti+nrowA-1
      lastj <- firstj+ncolA-1
      Adiag[firsti:lasti,firstj:lastj]<-A
      firsti <- lasti+1
      firstj <- lastj+1
    } 
    return (Adiag)
  }
  
  #cat("FitSpatioTemporalFH",model,": Start:",date(),"\n")
  
  result <- list(model=model, convergence=TRUE, iterations=0,
                 validtheta=FALSE, theta=0, beta=0, 
                 goodnessoffit=0, estimates=0)
  
  if (model!="A" && model!="B")
  {
    cat("Error REML_Model: Model must be A or B:",model)
    result$convergence <- FALSE
    return (result)
  }      
  
  M <- nD*nT
  nparam <- nrow(theta0)
  # Initialization
  invA   <- matrix(0,nrow=M,ncol=M)
  invAZ1 <- matrix(0,nrow=M,ncol=nD)
  tZ1PZ1 <- matrix(0,nrow=nD, ncol=nD) 
  tZ1P   <- matrix(0,nrow=nD,ncol=M)  
  S <- trPV   <- matrix(0,nrow=nparam,ncol=1)
  F <- trPVPV <- matrix(0,nrow=nparam,ncol=nparam)
  
  Va     <- list()
  ### Calculate Z1
  vector1T <- matrix(1,nrow=nT, ncol=1)
  Z1       <- matrix(0,nrow=nD*nT, ncol=nD)
  first <- 1
  for (d in 1:nD)
  {
    last <- first+nT-1
    Z1[first:last,d]<-vector1T
    first <- last+1
  } 
  
  tZ1 <- t(Z1)                   
  ty  <- t(y)
  tX  <- t(X)
  tWW <- crossprod(W)
  p1derivrho1 <- - W - t(W) 
  Id      <- diag(1,nrow=nD,ncol=nD)
  Tmen1   <- nT-1
  
  if (model=="A")
    Va[[3]] <- diag(1,nrow=M,ncol=M)                          
  else
  {
    PV <- list()                                                   
    Omega2drho2 <- derivOmega2drho2 <- matrix(0,nrow=nT,ncol=nT)     
    seqTmen1_1 <- sequence(Tmen1:1)                                
    Ve <- diag(sigma2dt)                                           
  }
  
  thetakmas1 <- thetak <- theta0
  k <- 0
  diff <- PRECISION+1
  while (diff>PRECISION & k<MAXITER)
  {
    k <- k+1
    thetak    <- thetakmas1
    sigma21_k <- thetak[1]
    rho1_k    <- thetak[2]
    sigma22_k <- thetak[3]
    
    #cat("It:", k, "Thetak:", thetak, "H:",date(),"\n")
    
    Omega1rho1_k <- solve(crossprod(Id-rho1_k*W))        
    Vu1 <- sigma21_k*Omega1rho1_k
    
    if (model=="A")
    {
      invAvec <- 1/(sigma22_k+sigma2dt)                      
      invA <- diag(invAvec)                                  
      first <- 1
      for (i in 1:nD)
      {
        last <- first+Tmen1
        firstlast <- first:last
        invAZ1[firstlast,i]<-invAvec[firstlast]
        first <- first + nT         
      }
    }
    else
    {
      rho2_k        <- thetak[4]                                      
      Unomenrho22_k <- 1-(rho2_k^2)                             
      Omega2drho2_k <- matrix(0,nrow=nT,ncol=nT)                       
      Omega2drho2_k[lower.tri(Omega2drho2_k)] <- rho2_k^seqTmen1_1
      Omega2drho2_k        <- Omega2drho2_k+t(Omega2drho2_k)    
      diag(Omega2drho2_k)  <- 1                                 
      Omega2drho2_k        <- (1/Unomenrho22_k)*Omega2drho2_k   
      sigma22Omega2drho2_k <- sigma22_k*Omega2drho2_k           
      
      #inverse in blocks
      first <- 1
      for (i in 1:nD)
      {
        last <- first+Tmen1
        firstlast <- first:last
        Ved <- Ve[first:last,first:last]
        Ad <- sigma22Omega2drho2_k + Ved
        invAd <- solve(Ad)              
        invA[first:last,first:last]<-invAd
        first <- first + nT         
      }
      invAZ1 = invA%*%Z1              
    }
    invVu1 <- solve(Vu1)
    
    invV <- invA - invAZ1%*%solve(invVu1+tZ1%*%invAZ1)%*%t(invAZ1)
    tXinvV     <- tX %*% invV
    inv_tXinVX <- solve(tXinvV %*% X) 
    P          <- invV - t(tXinvV) %*% inv_tXinVX %*% tXinvV
    
    # calculate S and F
    derivrho1_k <- p1derivrho1 + 2*rho1_k*tWW
    # calculate Va
    sigmaOmegaderivrho1Omega <- (-sigma21_k)*(Omega1rho1_k %*% 
                                                derivrho1_k %*% Omega1rho1_k)
    
    Va[[1]] <- Z1 %*% Omega1rho1_k %*% tZ1 
    Va[[2]] <- Z1 %*% sigmaOmegaderivrho1Omega %*% tZ1
    
    if (model=="A")
    {
      tZ1P    <- tZ1%*%P               
      tZ1PZ1  <- tZ1P%*%Z1           
      auxV1   <- tZ1PZ1%*%Omega1rho1_k       
      trPV[1] <- sum(diag(auxV1))   
      auxV2   <- tZ1PZ1%*%sigmaOmegaderivrho1Omega  
      trPV[2] <- sum(diag(auxV2))
      trPV[3] <- sum(diag(P))    
      Py      <- P %*% y                                    
      tyP     <- t(Py)                         # P is symmetric
      trPVPV[1,1] <- sum(auxV1*t(auxV1))
      trPVPV[2,2] <- sum(auxV2*t(auxV2))
      trPVPV[3,3] <- sum(P*t(P))        
      trPVPV[1,2] <- sum(auxV1*t(auxV2))
      tZ1PPZ1     <- tZ1P%*%t(tZ1P)         
      trPVPV[1,3] <- sum(tZ1PPZ1*t(Omega1rho1_k))
      trPVPV[2,3] <- sum(tZ1PPZ1*t(sigmaOmegaderivrho1Omega))
    }
    else
    {
      Va[[3]] <- diagonalizematrix(Omega2drho2_k,ntimes=nD)
      
      derivOmega2drho2_k <- matrix(0,nrow=nT,ncol=nT)
      derivOmega2drho2_k[lower.tri(derivOmega2drho2_k)]<-seqTmen1_1*
        rho2_k^(seqTmen1_1-1) 
      derivOmega2drho2_k <- derivOmega2drho2_k +
        t(derivOmega2drho2_k)
      derivOmega2drho2_k <- (1/Unomenrho22_k)*derivOmega2drho2_k +
        (2*rho2_k/Unomenrho22_k)*Omega2drho2_k
      sigma22derivOmega2drho2_k <- sigma22_k * derivOmega2drho2_k
      
      Va[[4]] <- diagonalizematrix(sigma22derivOmega2drho2_k,
                                   ntimes=nD)
      
      for (i in 1:nparam)
      {
        PV[[i]] <- P %*% Va[[i]]
        trPV[i] <- sum(diag(PV[[i]])) 
      }
      
      for (j in 1:nparam)
      {
        tPVj <- t(PV[[j]])
        for (i in 1:j)
          trPVPV[i,j] <- sum(PV[[i]]*tPVj)
      }
      Py  <- P %*% y  
      tyP <- t(Py)    
    }
    
    for (a in 1:nparam)
    {
      S[a] <- (-0.5)*trPV[a] + 0.5*(tyP %*% Va[[a]] %*% Py) 
      for (b in a:nparam)
        F[a,b] <- 0.5*trPVPV[a,b]
    }     
    
    for (a in 2:nparam)     # symmetric
    {
      for (b in 1:(a-1))
        F[a,b] <- F[b,a]
    }
    
    Finv <- ginv(F)                     
    thetakmas1 <- thetak + Finv %*% S 
    
    # Test values!=0 to avoid division errors
    if (any(thetak==0))
    {
      ##cat("Warning: Iteration ",k,". Thetak:",thetak )
      for (i in 1:nparam)
      {
        if (thetak[i]==0)
          thetak[i]<- 0.0001
      }
      ##cat("  New values:",thetak )
    }
    diff <- max( abs(thetak - thetakmas1)/thetak )
    
  } #while (diff>PRECISION & k<MAXITER)
  
  result$iterations <- k
  
  # validate output
  if (k>=MAXITER && diff>=PRECISION)
  {
    result$convergence <- FALSE
    return (result)
  }
  niter     <- k
  sigma21_k <- thetakmas1[1]
  rho1_k    <- thetakmas1[2]
  sigma22_k <- thetakmas1[3]
  
  ##cat("ThetakMAS1:", thetakmas1,"\n")
  
  if (sigma21_k<0 || rho1_k < (-1) || rho1_k>1 || sigma22_k<0 || 
        (model=="B" && (thetakmas1[4]<(-1) || thetakmas1[4]>1)) )
  {
    result$theta <- thetakmas1
    return(result)
  }
  
  # calculate estimates beta, u, mu
  
  Omega1rho1_k <- solve(crossprod(Id-rho1_k*W)) 
  Vu1 <- sigma21_k*Omega1rho1_k
  
  if (model=="A")
  {
    invAvec <- 1/(sigma22_k+sigma2dt)
    invA <- diag(invAvec)            
    first <- 1
    for (i in 1:nD)
    {
      last <- first+Tmen1
      firstlast <- first:last
      invAZ1[firstlast,i]<-invAvec[firstlast]
      first <- first + nT         
    }
  }
  else
  {
    rho2_k        <- thetakmas1[4] 
    Unomenrho22_k <- 1-(rho2_k^2)  
    Omega2drho2_k <- matrix(0,nrow=nT,ncol=nT)
    Omega2drho2_k[lower.tri(Omega2drho2_k)] <- rho2_k^seqTmen1_1
    Omega2drho2_k        <- Omega2drho2_k+t(Omega2drho2_k)      
    diag(Omega2drho2_k)  <- 1                              
    Omega2drho2_k        <- (1/Unomenrho22_k)*Omega2drho2_k
    sigma22Omega2drho2_k <- sigma22_k*Omega2drho2_k        
    
    first <- 1
    for (i in 1:nD)
    {
      last <- first+Tmen1
      firstlast <- first:last
      Ad <- sigma22Omega2drho2_k + Ve[first:last,first:last]
      invAd <- solve(Ad)                                    
      invA[first:last,first:last]<-invAd                    
      first <- first + nT         
    }
    invAZ1 <- invA%*%Z1
  }
  invVu1  <- solve(Vu1)
  invV    <- invA - invAZ1%*%solve(invVu1+tZ1%*%invAZ1)%*%t(invAZ1)
  
  tXinvV  <- tX %*% invV
  Q       <- solve(tXinvV %*% X)                  
  betaest <- Q %*% (tXinvV %*% y)             
  ymenXbetaest  <- ( y - X %*% betaest )
  invVymenXBest <- invV %*% ymenXbetaest
  
  parte1  <- Vu1 %*% tZ1
  if (model=="A")
    parte2 <- diag(sigma22_k,nrow=M)
  else
    parte2 <- diagonalizematrix(sigma22Omega2drho2_k,ntimes=nD)
  
  u1est   <- parte1 %*% invVymenXBest
  u2dtest <- parte2 %*% invVymenXBest
  u1dtest <- matrix(data=rep(u1est, each=nT),nrow=M,ncol=1) 
  mudtest <- X%*%betaest + u1dtest + u2dtest                   
  
  V       <- solve(invV)
  loglike <- (-0.5) * ( M*log(2*pi) +
                          determinant(V,logarithm=TRUE)$modulus +
                          t(ymenXbetaest)%*%invV%*%ymenXbetaest )
  AIC     <- (-2)*loglike + 2*(length(betaest)+nparam)
  BIC     <- (-2)*loglike + log(M)*(length(betaest)+nparam)
  
  # calculate confidence intervals and pvalues
  alfa <- 1-confidence
  k    <- 1-alfa/2
  z    <- qnorm(k)
  sqrtQvector <- sqrt(diag(Q))
  intconfidencebeta  <- z*sqrtQvector        
  intconfidencetheta <- z*sqrt(diag(Finv))  
  
  z <- abs(betaest)/sqrtQvector       
  p <- pnorm(z, lower.tail=FALSE)
  pvalue <- 2*p
  
  result$validtheta   <- TRUE
  result$theta        <- data.frame(estimate=thetakmas1,
                                    std.error=intconfidencetheta)
  result$beta         <- data.frame(coef=betaest,
                                    std.error=intconfidencebeta,
                                    tvalue=betaest/intconfidencebeta,
                                    pvalue=pvalue,
                                    greater.alfa=pvalue>alfa)
  result$goodnessoffit<- c(loglike=loglike, AIC=AIC, BIC=BIC)
  result$estimates    <- mudtest
  
  return (result)
}


