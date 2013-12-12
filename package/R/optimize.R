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

   u <- optimizeRE(modelSpecs$sigma, y = modelSpecs$y, X = modelSpecs$x, Z1 = modelSpecs$Z1, 
             sigmaSamplingError = modelSpecs$sigmaSamplingError, rho = modelSpecs$rho,
             W = modelSpecs$w, beta = modelSpecs$beta, K = modelSpecs$K, Z = modelSpecs$Z, tol = modelSpecs$tol, maxit = modelSpecs$maxIter)
   
   modelSpecs$u <- as.numeric(u)
  
  return(modelSpecs)
}
