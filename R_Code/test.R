rm(list= ls())

# Daten generieren:

# Installiere Paket
# utils::install.packages(pkgs="../spatioTemporalData/spatioTemporalData_1.0.zip", 
#                         repos = NULL)

require(spatioTemporalData)
require(SAE)
set.seed(1)
setup <- simSetupMarhuenda(nDomains=100, nTime=10, sarCorr=0.5, arCorr=0.5)
output <- simRunMarhuenda(setup)
dat <- slot(output[[1]], "data")[[1]]
sigmaE <- slot(output[[1]], "sigma")

# Prepare Data
dat <- dat[order(dat$Domain, dat$Time), ]
XY <- makeXY(y ~ x, dat)

getNDomains <- function(dat) length(unique(dat$Domain))
getNTime <- function(dat) length(unique(dat$Time))

modelSpecs <- list(y = XY$y,
                   x = XY$x,
                   nDomains = getNDomains(dat),
                   nTime = getNTime(dat),
                   
                   w0 = w0Matrix(getNDomains(dat)),
                   w = wMatrix(getNDomains(dat)),
                   
                   Z = reZ(getNDomains(dat), getNTime(dat)),
                   Z1 = reZ1(getNDomains(dat), getNTime(dat)),
                   beta = c(1, 1),
                   sigma = c(1, 1),
                   rho = c(0.5, 0.5),
                   sigmaE = sigmaE,
                   tol = 1e-3,
                   psiFunction = MASS::psi.huber,
                   K = 0.71)


#' optimizeBeta
#' 
#' @description find optimal beta coefficients for given variance parameters
#' 
#' @param modelSpecs list with all necessary components for estimation
optimizeBeta <- function(modelSpecs) {
  #update necessary components
  Ome1 <- updateOmega1(sarCorr=modelSpecs$rho[1], w0=modelSpecs$w0)
  Ome2 <- updateOmega2(arCorr=modelSpecs$rho[2], nTime=modelSpecs$nTime)
  A <- updateA(sigma2 = modelSpecs$sigma[2], Ome2=Ome2, nDomains = modelSpecs$nDomains, nTime= modelSpecs$nTime)
  V <- updateV(sigma1=modelSpecs$sigma[1], Ome1=Ome1, A=A, Z1=modelSpecs$Z1)
  Vinv <- qr.solve(V)
  #Vinv <- updateSolvedV(sarCorr=, sigma1=, arCorr=, A=, Ome1=, Z1=, )
  sqrtU <- updateSqrtU(V=V)
  sqrtUinv <- diag(1/diag(sqrtU))
  
  # Some precalculations:
  tmp1 <- crossprod(modelSpecs$x, Vinv)
  tmp2 <- tmp1 %*% sqrtU
  
  # Initilize vectors for beta coefficients:
  newBeta <- modelSpecs$beta
  beta <- modelSpecs$beta
  iter <- 1
  
  # Begin NR-Algorithm - see Issue 1 - Paper - Numerical Stability
  while(all((newBeta - beta)^2 > modelSpecs$tol) || iter == 1) {
    beta <- newBeta
    resid <- sqrtUinv %*% (modelSpecs$y - modelSpecs$x %*% beta)
    dOfBeta <- diag(as.numeric(modelSpecs$psiFunction(u = resid, deriv = 1)))
    tmp3 <- tmp1 %*% dOfBeta %*% modelSpecs$x
    #tmp3 <- tmp1 %*% modelSpecs$x
    tmp4 <- solve(tmp3)
    tmp5 <- tmp2 %*% modelSpecs$psiFunction(u = resid)
    
    newBeta <- beta + tmp4 %*% tmp5
    cat(iter)
    iter <- iter + 1 
  }
  modelSpecs$beta <- newBeta
  return(modelSpecs)
}

#' optimizeSigma
#' 
#' @description find optimal sigma coefficients for given variance parameters
#' 
#' @param modelSpecs list with all necessary components for estimation
optimizeSigma <- function(modelSpecs) {
  
  updateDerivatives <- updateDerivativesClosure(modelSpecs$nDomains, modelSpecs$nTime, Z1=modelSpecs$Z1, W=modelSpecs$w)
  
  optimizerClosure <- function(modelSpecs) {
    force(modelSpecs)
    function(sigmas) {
      
      modelSpecs$sigma <<- sigmas
      
      Ome1 <- updateOmega1(sarCorr=modelSpecs$rho[1], w0=modelSpecs$w0)
      Ome2 <- updateOmega2(arCorr=modelSpecs$rho[2], nTime=modelSpecs$nTime)
      A <- updateA(sigma2 = modelSpecs$sigma[2], Ome2=Ome2, nDomains = modelSpecs$nDomains, nTime= modelSpecs$nTime)
      V <- updateV(sigma1=modelSpecs$sigma[1], Ome1=Ome1, A=A, Z1=modelSpecs$Z1)
      Vinv <- qr.solve(V)
      #Vinv <- updateSolvedV(sarCorr=modelSpecs$rho[1], sigma1=, arCorr=, A=, Ome1=, Z1=, )
      sqrtU <- updateSqrtU(V=V)
      sqrtUinv <- diag(1/diag(sqrtU))
      
      resid <- sqrtUinv %*% (modelSpecs$y - modelSpecs$x %*% modelSpecs$beta)
      phiR <- modelSpecs$psiFunction(u = resid)
      
      derivatives <- updateDerivatives(sarCorr=modelSpecs$rho[1], 
                                       sigma1 = modelSpecs$sigma[1],
                                       arCorr = modelSpecs$rho[2],
                                       sigma2 = modelSpecs$sigma[2],
                                       Ome1 = Ome1,
                                       Ome2 = Ome2,
                                       parSet = "sigma")
      tmp1 <- phiR %*% sqrtU %*% Vinv
      tmp2 <- Vinv %*% sqrtU %*% phiR
      
      tmpSig1 <- sum(diag(modelSpecs$K * Vinv %*% derivatives$derVSigma1))
      tmpSig2 <- sum(diag(modelSpecs$K * Vinv %*% derivatives$derVSigma2))
      
      optSig1 <- tmp1 %*% derivatives$derVSigma1 %*% tmp2 - tmpSig1
      optSig2 <- tmp1 %*% derivatives$derVSigma2 %*% tmp2 - tmpSig2
      
      return(optSig1^2 + optSig2^2)
    }
  }
    
  modelSpecs$sigma <- optim(par=modelSpecs$sigma, fn = optimizerClosure(modelSpecs))$par
  return(modelSpecs)
}

#' optimizeRho
#' 
#' @description find optimal rho coefficients for given variance parameters
#' 
#' @param modelSpecs list with all necessary components for estimation
optimizeRho <- function(modelSpecs) {
  
  updateDerivatives <- updateDerivativesClosure(modelSpecs$nDomains, modelSpecs$nTime, Z1=modelSpecs$Z1, W=modelSpecs$w)
  
  optimizerClosure <- function(modelSpecs) {
    force(modelSpecs)
    function(rho) {
      modelSpecs$rho <<- rho
      Ome1 <- updateOmega1(sarCorr=modelSpecs$rho[1], w0=modelSpecs$w0)
      Ome2 <- updateOmega2(arCorr=modelSpecs$rho[2], nTime=modelSpecs$nTime)
      A <- updateA(sigma2 = modelSpecs$sigma[2], Ome2=Ome2, nDomains = modelSpecs$nDomains, nTime= modelSpecs$nTime)
      V <- updateV(sigma1=modelSpecs$sigma[1], Ome1=Ome1, A=A, Z1=modelSpecs$Z1)
      Vinv <- qr.solve(V)
      #Vinv <- updateSolvedV(sarCorr=modelSpecs$rho[1], sigma1=, arCorr=, A=, Ome1=, Z1=, )
      sqrtU <- updateSqrtU(V=V)
      sqrtUinv <- diag(1/diag(sqrtU))
      
      resid <- sqrtUinv %*% (modelSpecs$y - modelSpecs$x %*% modelSpecs$beta)
      phiR <- modelSpecs$psiFunction(u = resid)
      
      derivatives <- updateDerivatives(sarCorr=modelSpecs$rho[1], 
                                       sigma1 = modelSpecs$sigma[1],
                                       arCorr = modelSpecs$rho[2],
                                       sigma2 = modelSpecs$sigma[2],
                                       Ome1 = Ome1,
                                       Ome2 = Ome2,
                                       parSet = "rho")
      tmp1 <- phiR %*% sqrtU %*% Vinv
      tmp2 <- Vinv %*% sqrtU %*% phiR
      
      tmpSig1 <- sum(diag(modelSpecs$K * Vinv %*% derivatives$derVSarCorr))
      tmpSig2 <- sum(diag(modelSpecs$K * Vinv %*% derivatives$derVArCorr))
      
      optSig1 <- tmp1 %*% derivatives$derVSarCorr %*% tmp2 - tmpSig1
      optSig2 <- tmp1 %*% derivatives$derVArCorr %*% tmp2 - tmpSig2
      
      return(optSig1^2 + optSig2^2)
    }
  }
  
  modelSpecs$rho <- optim(par=modelSpecs$rho, fn = optimizerClosure(modelSpecs))$par
  return(modelSpecs)
}

#' optimizeRho
#' 
#' @description find optimal rho coefficients for given variance parameters
#' 
#' @param modelSpecs list with all necessary components for estimation
#' 
optimizeParameters <- function(modelSpecs) {
  for (i in 1:2) {
    cat(paste("beta = ", modelSpecs$beta, "sigma = ", modelSpecs$sigma, "rho = ", modelSpecs$rho, "\n"))
    modelSpecs <- optimizeBeta(modelSpecs)
    modelSpecs <- optimizeRho(modelSpecs)
    modelSpecs <- optimizeSigma(modelSpecs)
    
  }
  return(modelSpecs)
}

tmp <- optimizeBeta(modelSpecs)

tmp <- optimizeBeta(tmp)
summary(lm(y~x, data = dat))
tmp$beta

