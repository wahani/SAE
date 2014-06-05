optimizeParam <- function(modelSpecs) {
  # Generic function: fit model parameters
  UseMethod("optimizeParam")
}

optimizeParam.default <- function(modelSpecs) {
  # Generic function: fit model parameters
  stop("This type is not supported!")
}

optimizeParam.MSRFH <- function(modelSpecs) {
  # Generic function: fit model parameters
  algorithmExpr <- expression({modelSpecs <- optimizeBeta(modelSpecs)
                               modelSpecs <- optimizeReVar(modelSpecs)})
  modelSpecs <- algorithmFH(algorithmExpr, modelSpecs)
  modelSpecs
  
}

algorithmFH <- function(expr, modelSpecs, paramNames = c("beta", "reVar")) {
  # Helper-Functions
  
  isTolReached <- function(x, y) {
    all(abs((x - y)/x) < modelSpecs$tol)
  }
  
  oldParam <- unlist(lapply(paramNames, get, envir = modelSpecs)) + 1
  i <- 0
  while(i < modelSpecs$maxIter) {
    i <- i + 1
    eval(expr=expr)
    newParam <- unlist(lapply(paramNames, get, envir = modelSpecs))
    if(isTolReached(newParam, oldParam)) break else oldParam <- newParam
  }
  
  modelSpecs$fitparam$m <- sort(rep(1:i, length(unique(modelSpecs$fitparam$param))))
      
  modelSpecs
}