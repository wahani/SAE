optimizeSigma <- function(modelSpecs) {
  
  modelSpecs$sigma <- fuTools::fp(optimizerSigma, 
                         modelSpecs$sigma, 
                         opts = list(tol = modelSpecs$tol, 
                                     maxiter = modelSpecs$maxIter),
                         y = modelSpecs$y, X = modelSpecs$X, Z1 = modelSpecs$Z1, 
                         sigmaSamplingError = modelSpecs$sigmaSamplingError, rho = modelSpecs$rho,
                         W = modelSpecs$w, beta = modelSpecs$beta, K = modelSpecs$K, Z = modelSpecs$Z)$x
  
  return(modelSpecs)
}

optimizeReVar <- function(modelSpecs) {
  # Generic function: computes varianze Parameters of random effects
  UseMethod("optimizeReVar")
}

optimizeReVar.default <- function(modelSpecs) {
  # Generic function: computes varianze Parameters of random effects
  stop("This type is not supported!")
}

optimizerReVarRFH <- function(reVar, vardir, y, X, beta, k, K, psiFunction) {
  
  # V
  V <- diag(vardir + reVar)
  Vinv <- diag(1 / diag(V))
  
  # G
  G <- diag(reVar, nrow = nrow(V))
  Ginv <- diag(1 / diag(G))
  
  # U
  sqrtU <- updateSqrtU(V=V)
  sqrtUinv <- diag(1 / diag(sqrtU))
  
  # residuals
  resid <- sqrtUinv %*% (y - X %*% beta)
  psiResid <- psiFunction(u = resid, k = k)
  
  # a + A
  a <- t(psiResid) %*% sqrtU %*% Vinv %*% Vinv %*% sqrtU %*% psiResid
  A <- matTrace(K * Vinv %*% Ginv)
  
  est <- as.numeric(a / A)
  
  if(any(is.na(est)) || any(est == Inf)) {
    warning("Varince Parameter is not identified and is set to 0")
    est[is.na(est)] <- 0
    est[est==Inf] <- 0
  }
  est
  
}

optimizeReVar.MSRFH <- function(modelSpecs) {
  startTime <- proc.time()[3]
  # Generic function: computes variance Parameters of random effects
  fit <- fuTools::fp(optimizerReVarRFH, modelSpecs$reVar, opts = list(tol = modelSpecs$tol, 
                                                             maxiter = modelSpecs$maxIter),
            vardir = modelSpecs$vardir, y = modelSpecs$y, X = modelSpecs$X, beta = modelSpecs$beta,
            k = modelSpecs$k, K = modelSpecs$K, modelSpecs$psiFunction)
  
  modelSpecs$reVar <- fit$x
    
  # Reporting for algorithm
  usedTime <- proc.time()[3] - startTime
  n <- NROW(modelSpecs$fitparam)
  modelSpecs$fitparam[n + 1, c("param", "m", "stepIterations", "returnStatus", "timeElapsed")] <- 
    data.frame("reVar", 1, fit$iterations, fit$returnStatus, usedTime,
               stringsAsFactors=FALSE)
  modelSpecs$fitparam$stepParam[n + 1] <- list(fit$x)
  
  modelSpecs
}
