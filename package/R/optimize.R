#' optimizeBeta
#' 
#' @description find optimal beta coefficients for given variance parameters
#' 
#' @param modelSpecs list with all necessary components for estimation
optimizeBeta <- function(modelSpecs) {
  #update necessary components
  
  listV <- matVinv(W=modelSpecs$w, rho1=modelSpecs$rho[1], sigma1=modelSpecs$sigma[1],
                   rho2 = modelSpecs$rho[2], sigma2 = modelSpecs$sigma[2], Z1=modelSpecs$Z1,
                   modelSpecs$sigmaSamplingError)
  
  V <- listV$V
  Vinv <- listV$Vinv
  
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
  while(all((newBeta - beta)^2 > modelSpecs$tol) || iter == 1 & iter < modelSpecs$maxIter) {
    
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
  
  optimizerWrapper <- function(sigma) {
    optimizerSigma(sigma, y = modelSpecs$y, X = modelSpecs$x, Z1 = modelSpecs$Z1, 
                   sigmaSamplingError = modelSpecs$sigmaSamplingError, rho = modelSpecs$rho,
                   W = modelSpecs$w, beta = modelSpecs$beta, K = modelSpecs$K, Z = modelSpecs$Z)
  }
  
  modelSpecs$sigma <- fp(optimizerWrapper, 
                         modelSpecs$sigma, 
                         opts = list(tol = modelSpecs$tol, 
                                     maxiter = modelSpecs$maxIter))$x
  
  return(modelSpecs)
}

#' optimizeRho
#' 
#' @description find optimal rho coefficients for given variance parameters
#' 
#' @param modelSpecs list with all necessary components for estimation
optimizeRho <- function(modelSpecs) {
  modelSpecs$rho <- nloptr(modelSpecs$rho, 
                           optimizerRho,
                           lb = c(-0.99, -0.99), ub = c(0.99, 0.99),
                           opts = list(algorithm = "NLOPT_LN_NELDERMEAD",
                                       xtol_rel = modelSpecs$tol,
                                       maxeval = modelSpecs$maxIter),
                           y = modelSpecs$y, X = modelSpecs$x, sigma1 = modelSpecs$sigma[1], 
                           sigma2 = modelSpecs$sigma[2], Z1 = modelSpecs$Z1, 
                           sigmaSamplingError = modelSpecs$sigmaSamplingError, 
                           W = modelSpecs$w, beta = modelSpecs$beta, K = modelSpecs$K)$solution
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
  
  iter <- 1
  while (checkCriterion(modelSpecs, oldParams) & iter < modelSpecs$maxIter) {
    #cat(paste("beta = ", modelSpecs$beta, "sigma = ", modelSpecs$sigma, "rho = ", modelSpecs$rho, "\n"))
    oldParams <- c(modelSpecs$beta, modelSpecs$sigma, modelSpecs$rho)
    consoleOutput(modelSpecs$consoleOutput)
    modelSpecs <- optimizeBeta(modelSpecs)
    modelSpecs <- optimizeSigma(modelSpecs)
    modelSpecs <- optimizeRho(modelSpecs)
    
    iter <- iter + 1
  }
  return(modelSpecs)
}


#' estimateRE
#' 
#' @description This function estimates the BLUP for the RE-Part, given the results
#' of the estimation of all parameters - \code{\link{optimizeParameters}}
#' 
#' @param modelSpecs list with all necessary components for estimation
#' 
estimateRE <- function(modelSpecs) {
  n <- modelSpecs$nDomains*modelSpecs$nTime
  # Sampling Error Component
  R <- diag(modelSpecs$sigmaSamplingError)
  svdR <- svd(R)
  sqrtRinv <- solve(t(svdR$v%*%(t(svdR$u)*sqrt(svdR$d))))
  
  # RE - Spatio-Temporal
  Ome1 <- updateOmega1(sarCorr=modelSpecs$rho[1], w0=modelSpecs$w0)
  Ome2 <- updateOmega2(arCorr=modelSpecs$rho[2], nTime=modelSpecs$nTime)
  
  G <- matrix(0, ncol = modelSpecs$nDomains + modelSpecs$nDomains * modelSpecs$nTime,
              nrow = modelSpecs$nDomains + modelSpecs$nDomains * modelSpecs$nTime)
  
  G[1:modelSpecs$nDomains, 1:modelSpecs$nDomains] <- modelSpecs$sigma[1] * Ome1
  G[(modelSpecs$nDomains+1):(modelSpecs$nDomains * modelSpecs$nTime + modelSpecs$nDomains), 
    (modelSpecs$nDomains+1):(modelSpecs$nDomains * modelSpecs$nTime + modelSpecs$nDomains)] <- modelSpecs$sigma[2] * omega2Diag(Ome2=Ome2, nDomains=modelSpecs$nDomains)
  
  svdG <- svd(G)
  sqrtGinv <- solve(t(svdG$v%*%(t(svdG$u)*sqrt(svdG$d))))
  
  # Variance-Covariance
  Z <- modelSpecs$Z
  listV <- matVinv(W=modelSpecs$w, rho1=modelSpecs$rho[1], sigma1=modelSpecs$sigma[1],
                   rho2 = modelSpecs$rho[2], sigma2 = modelSpecs$sigma[2], Z1=modelSpecs$Z1,
                   modelSpecs$sigmaSamplingError)
  
  V <- listV$V
  Vinv <- listV$Vinv
  
  #   sqrtU <- updateSqrtU(V=V)
  #   sqrtUinv <- diag(1/diag(sqrtU))
  
  # Starting Values
  y <- modelSpecs$y
  X <- modelSpecs$x
  beta <- modelSpecs$beta
  
  resid <- y - X %*% beta
  vv <- G %*% t(Z) %*% Vinv %*% as.vector(resid)
  
  # Algorithm
  n <- modelSpecs$nDomains*modelSpecs$nTime
  areanumber <- modelSpecs$nDomains
  tol <- modelSpecs$tol
  diff.u<-1
  
  i<-0
  k_v <- 1.345
  maxit <- modelSpecs$maxIter
  while (abs(diff.u)>tol)
  {
    i <- i+1 
    v_robust <- as.vector(vv)
    res1 <- sqrtRinv %*% (resid - Z %*% v_robust)
    res2 <- sqrtGinv %*% v_robust
    w2 <- diag(c(psiOne(res1, k_v)/res1), n, n)
    w3 <- diag(c(psiOne(res2, k_v)/res2), n + areanumber, n + areanumber)
    Atmp1 <- t(Z) %*% (sqrtRinv) %*% w2 %*%(sqrtRinv) %*% Z
    Atmp2 <- sqrtGinv %*% w3 %*% sqrtGinv
    A <- Atmp1 + Atmp2
    B <- t(Z) %*% sqrtRinv %*% w2 %*% sqrtRinv %*% resid
    vv <- solve(A) %*% B
    
    diff.u<-sum((vv-v_robust)^2)
    if (i > maxit) break
  }
  
  modelSpecs$u <- as.numeric(Z %*% as.numeric(vv))
  return(modelSpecs)
}
