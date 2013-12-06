optimizeSigma <- function(modelSpecs) {
  
  updateDerivatives <- updateDerivativesClosure(modelSpecs$nDomains, modelSpecs$nTime, Z1=modelSpecs$Z1, W=modelSpecs$w)
  
  optimizerClosure <- function(modelSpecs) {
    force(modelSpecs)
    function(sigmas) {
      
      modelSpecs$sigma <<- sigmas[c(1,3)]
      modelSpecs$rho[1] <<- sigmas[2]
      
      Ome1 <- updateOmega1(sarCorr=modelSpecs$rho[1], w0=modelSpecs$w0)
      Ome2 <- updateOmega2(arCorr=modelSpecs$rho[2], nTime=modelSpecs$nTime)
      
      listV <- matVinv(W=modelSpecs$w, rho1=modelSpecs$rho[1], sigma1=modelSpecs$sigma[1],
                       rho2 = modelSpecs$rho[2], sigma2 = modelSpecs$sigma[2], Z1=modelSpecs$Z1,
                       modelSpecs$sigmaSamplingError)
      V <- listV$V
      Vinv <- listV$Vinv
      
      sqrtU <- updateSqrtU(V=V)
      sqrtUinv <- diag(1/diag(sqrtU))
      
      resid <- sqrtUinv %*% (modelSpecs$y - modelSpecs$x %*% modelSpecs$beta)
      phiR <- modelSpecs$psiFunction(u = resid)
      
      derivatives <- list(derVSigma1 = matVderS1(Ome1=Ome1, Z1 = modelSpecs$Z1),
                          derVSigma2 = matVderS2(Ome2=Ome2, nDomains=modelSpecs$nDomains),
                          derVRho1 = matVderR1(rho1=modelSpecs$rho[1], sigma1=modelSpecs$sigma[1], 
                                               Z1=modelSpecs$Z1,Ome1=Ome1, W=modelSpecs$w))
      
      ZVuZt <- V 
      diag(ZVuZt) <- diag(ZVuZt) - modelSpecs$sigmaSamplingError
      
      ZVuZtinv <- chol2inv(chol(ZVuZt))
      
      Ome1SigmaBar <- Ome1RhoBar <- Ome2Bar <- matrix(0, nrow = modelSpecs$nDomains * modelSpecs$nTime + modelSpecs$nDomains,
                                                      ncol = modelSpecs$nDomains * modelSpecs$nTime + modelSpecs$nDomains)
      
      Ome1SigmaBar[1:modelSpecs$nDomains, 1:modelSpecs$nDomains] <- Ome1 %*% Ome1 - modelSpecs$rho[1] * Ome1 %*% Ome1 %*% t(modelSpecs$w)
      Ome1RhoBar[1:modelSpecs$nDomains, 1:modelSpecs$nDomains] <- -modelSpecs$sigma[1] * Ome1 %*% Ome1 %*% modelSpecs$w + modelSpecs$sigma[1] * modelSpecs$rho[1] * Ome1 %*% Ome1 %*% crossprod(modelSpecs$w)
      Ome2Bar[(modelSpecs$nDomains + 1):NROW(Ome2Bar), (modelSpecs$nDomains + 1):NROW(Ome2Bar)] <- omega2Diag(Ome2, modelSpecs$nDomains)
      
      Vu <- modelSpecs$sigma[1] * Ome1SigmaBar + modelSpecs$rho[1] * Ome1RhoBar + modelSpecs$sigma[2] * Ome2Bar
      all.equal(ZVuZt, modelSpecs$Z %*% Vu %*% t(modelSpecs$Z))
      
      tmp1 <- crossprod(phiR, sqrtU) %*% Vinv
      
      a <- c(tmp1 %*% tcrossprod(derivatives$derVSigma1, tmp1),
             tmp1 %*% tcrossprod(derivatives$derVRho1, tmp1),
             tmp1 %*% tcrossprod(derivatives$derVSigma2, tmp1))
      
      gamma1 <- modelSpecs$K * Vinv %*% derivatives$derVSigma1 %*% ZVuZtinv
      gamma2 <- modelSpecs$K * Vinv %*% derivatives$derVRho1 %*% ZVuZtinv
      gamma3 <- modelSpecs$K * Vinv %*% derivatives$derVSigma2 %*% ZVuZtinv
      
      tmp1 <- modelSpecs$Z %*% tcrossprod(Ome1SigmaBar, modelSpecs$Z)
      tmp2 <- modelSpecs$Z %*% tcrossprod(Ome1RhoBar, modelSpecs$Z)
      tmp3 <- modelSpecs$Z %*% tcrossprod(Ome2Bar, modelSpecs$Z)
      
      matTrace <- function(x) sum(diag(x))
      A <- diag(nrow = 3)
      A[1,1] <- matTrace(gamma1 %*% tmp1)
      A[1,2] <- matTrace(gamma1 %*% tmp2)
      A[1,3] <- matTrace(gamma1 %*% tmp3)
      A[2,1] <- matTrace(gamma2 %*% tmp1)
      A[2,2] <- matTrace(gamma2 %*% tmp2)
      A[2,3] <- matTrace(gamma2 %*% tmp3)
      A[3,1] <- matTrace(gamma3 %*% tmp1)
      A[3,2] <- matTrace(gamma3 %*% tmp2)
      A[3,3] <- matTrace(gamma3 %*% tmp3)
      
      return(solve(A) %*% a)
    }
  }
  browser()
  est <- fp(optimizerClosure(modelSpecs), 
            c(modelSpecs$sigma[1], modelSpecs$rho[1], modelSpecs$sigma[1]), 
            opts = list(tol = modelSpecs$tol, 
                        maxiter = modelSpecs$maxIter))$x
  
  modelSpecs$sigma <- est[c(1,3)]
  modelSpecs$rho[1] <- est[2]
  
  
  return(modelSpecs)
}