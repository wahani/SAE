optimizeBeta <- function(modelSpecs) {
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
  iter <- 1
  
  # Begin NR-Algorithm - see Issue 1 - Paper - Numerical Stability
  while(all((newBeta - beta)^2 > modelSpecs$tol) || iter == 1 & iter < modelSpecs$maxIter) {
    
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
  return(modelSpecs)
}