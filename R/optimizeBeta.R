optimizeBeta <- function(modelSpecs) {
  
  startTime <- proc.time()[3]
  
  #update necessary components
  listV <- compV(modelSpecs)
  
  V <- listV$V
  Vinv <- listV$Vinv
  
  sqrtU <- updateSqrtU(V=V)
  sqrtUinv <- diag(1/diag(sqrtU))
  
  # Some precalculations:
  tmp1 <- crossprod(modelSpecs$X, Vinv)
  tmp2 <- tmp1 %*% sqrtU
  
  # Initilize vectors for beta coefficients:
  newBeta <- modelSpecs$beta
  beta <- modelSpecs$beta
  iter <- 0
  
  if(modelSpecs$progress) cat("\n--Optimizing beta--\n")
  # Begin NR-Algorithm - see Issue 1 - Paper - Numerical Stability
  while(all((abs(newBeta - beta)/newBeta) > modelSpecs$tol) & iter < modelSpecs$maxIter | iter == 0) {
    
    if(modelSpecs$progress) cat("Iteration:", iter, "\n")
    beta <- newBeta
    
    resid <- sqrtUinv %*% (modelSpecs$y - modelSpecs$X %*% beta)
    dOfBeta <- diag(as.numeric(modelSpecs$psiFunction(u = resid, deriv = 1)))
    tmp3 <- tmp1 %*% dOfBeta %*% modelSpecs$X
    tmp4 <- solve(tmp3)
    tmp5 <- tmp2 %*% modelSpecs$psiFunction(u = resid, k = modelSpecs$k)
    
    newBeta <- beta + tmp4 %*% tmp5
    #cat(iter)
    iter <- iter + 1 
  }
  
  modelSpecs$beta <- newBeta
  
  # Reporting for algorithm
  usedTime <- proc.time()[3] - startTime
  convergence <- if(iter == modelSpecs$maxIter & 
                      all((abs(newBeta - beta)/newBeta) > modelSpecs$tol)) 1 else 0
  n <- NROW(modelSpecs$fitparam)
  modelSpecs$fitparam[n + 1, c("param", "m", "stepIterations", "returnStatus", "timeElapsed")] <- 
    data.frame("beta", 1, iter, convergence, usedTime, stringsAsFactors=FALSE)
  modelSpecs$fitparam$stepParam[n + 1] <- list(newBeta)
  
  return(modelSpecs)
}