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
    #cat(iter)
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
      tmp1 <- crossprod(phiR, sqrtU) %*% Vinv
      tmp2 <- Vinv %*% sqrtU %*% phiR
      
      tmpSig1 <- sum(diag(modelSpecs$K * Vinv %*% derivatives$derVSigma1))
      tmpSig2 <- sum(diag(modelSpecs$K * Vinv %*% derivatives$derVSigma2))
      
      optSig1 <- tmp1 %*% derivatives$derVSigma1 %*% tmp2 - tmpSig1
      optSig2 <- tmp1 %*% derivatives$derVSigma2 %*% tmp2 - tmpSig2
      
      return(optSig1^2 + optSig2^2)
    }
  }
  
  modelSpecs$sigma <- optim(par=modelSpecs$sigma, 
                            fn = optimizerClosure(modelSpecs),
                            method = modelSpecs$method)$par
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
      
      tmp1 <- crossprod(phiR, sqrtU) %*% Vinv
      tmp2 <- Vinv %*% sqrtU %*% phiR
      
      tmpSig1 <- sum(diag(modelSpecs$K * Vinv %*% derivatives$derVSarCorr))
      tmpSig2 <- sum(diag(modelSpecs$K * Vinv %*% derivatives$derVArCorr))
      
      optSig1 <- tmp1 %*% derivatives$derVSarCorr %*% tmp2 - tmpSig1
      optSig2 <- tmp1 %*% derivatives$derVArCorr %*% tmp2 - tmpSig2
      
      return(optSig1^2 + optSig2^2)
    }
  }
  
  modelSpecs$rho <- optim(par=modelSpecs$rho, 
                          fn = optimizerClosure(modelSpecs),
                          method = modelSpecs$method)$par
  return(modelSpecs)
}

#' optimizeParameters
#' 
#' @description find optimal Parameters coefficients for given variance parameters. 
#' This is a an "iterative wraper" of \code{\link{optimizeBeta}}, \code{\link{optimizeRho}}
#' and \code{\link{optimizeSigma}}
#' 
#' @param modelSpecs list with all necessary components for estimation
#' 
optimizeParameters <- function(modelSpecs) {
  
  checkCriterion <- function(modelSpecs, oldParams) 
    !all((c(modelSpecs$beta, modelSpecs$sigma, modelSpecs$rho) - oldParams)^2 < modelSpecs$tol)
  
  oldParams <- rep(100000, length(modelSpecs$beta) + 4)
  
  while (checkCriterion(modelSpecs, oldParams)) {
    #cat(paste("beta = ", modelSpecs$beta, "sigma = ", modelSpecs$sigma, "rho = ", modelSpecs$rho, "\n"))
    oldParams <- c(modelSpecs$beta, modelSpecs$sigma, modelSpecs$rho)
    cat(".")
    modelSpecs <- optimizeBeta(modelSpecs)
    modelSpecs <- optimizeRho(modelSpecs)
    modelSpecs <- optimizeSigma(modelSpecs)  
  }
  return(modelSpecs)
}