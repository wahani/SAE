#' optimizeSigma
#' 
#' @description find optimal sigma coefficients for given variance parameters
#' 
#' @param modelSpecs list with all necessary components for estimation
optimizeSigma <- function(modelSpecs) {
  
  modelSpecs$sigma <- fp(optimizerSigma, 
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
  Vinv <- diag(1 / (vardir + reVar))
  
  # G
  G <- diag(reVar, nrow = nrow(V))
  Ginv <- diag(1 / reVar, nrow = nrow(V))
  
  # U
  sqrtU <- updateSqrtU(V=V)
  sqrtUinv <- diag(1/diag(sqrtU))
  
  # residuals
  resid <- sqrtUinv %*% (y - X %*% beta)
  psiResid <- psiFunction(u = resid, k = k)
  
  # a + A
  a <- t(psiResid) %*% sqrtU %*% Vinv %*% Vinv %*% sqrtU %*% psiResid
  A <- matTrace(K * Vinv %*% Ginv)
  
  as.numeric(a / A)
  
}

optimizeReVar.MSRFH <- function(modelSpecs) {
  # Generic function: computes varianze Parameters of random effects
  modelSpecs$reVar <- fp(optimizerReVarRFH, modelSpecs$reVar, opts = list(tol = modelSpecs$tol, 
                                                      maxiter = modelSpecs$maxIter),
     vardir = modelSpecs$vardir, y = modelSpecs$y, X = modelSpecs$X, beta = modelSpecs$beta,
     k = modelSpecs$k, K = modelSpecs$K, modelSpecs$psiFunction)$x
  modelSpecs
}







