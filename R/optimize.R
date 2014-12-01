#' optimizeRho
#' 
#' @description find optimal rho coefficients for given variance parameters
#' 
#' @param modelSpecs list with all necessary components for estimation
optimizeRho <- function(modelSpecs) {
  modelSpecs$rho <- nloptr::nloptr(modelSpecs$rho, 
                           optimizerRho,
                           lb = c(-0.99, -0.99), ub = c(0.99, 0.99),
                           opts = list(algorithm = "NLOPT_LN_NELDERMEAD",
                                       xtol_rel = modelSpecs$tol,
                                       maxeval = modelSpecs$maxIter),
                           y = modelSpecs$y, X = modelSpecs$X, sigma1 = modelSpecs$sigma[1], 
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
  
  u <- optimizeRE(modelSpecs)
  modelSpecs$u <- as.numeric(u)
  
  return(modelSpecs)
}
