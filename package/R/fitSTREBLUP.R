#' fitSTREBLUP
#' 
#' @description Wrapper function for the estimation of the spatio-temporal-robust-EBLUP
#' 
#' @param formula formula for the fixed effects part
#' @param dat data with dependent and independent variable as they are used in formula. The data should 
#' contain a variable \code{Domain} and \code{Time} indicating the domain and time period for each observation respectively
#' @param beta, sigma, rho start values for the parameters
#' @param sigmaSamplingError True variances for the sampling errors
#' @param nDomains the number of domains. Will be determined by default, see \code{dat}
#' @param nTime the number of time periods. Will be determined by default, see \code{dat}
#' @param w0 proximity matrix
#' @param w at the moment row-standardized proximity matrix
#' @param tol numerical tolerance for algorithm
#' @param method method used by the function \code{\link{optim}}
#' @param maxIter maximal number of iterations - ignored in the function optim. Relevant for "global" algorithm and for the
#' optimization of beta
#' 
#' @export
#'
#' @examples require(spatioTemporalData); require(SAE)
#' set.seed(3)
#' setup <- simSetupMarhuenda(nDomains=30, nTime=5, sarCorr=0.5, arCorr=0.5, n = 200)
#' output <- simRunMarhuenda(setup) 
#' dat <- slot(output[[1]], "data")[[1]]
#' result <- fitSTREBLUP(y~x, dat, c(0,1), c(1,1), c(0.5,0.5))
#' summary(result)
fitSTREBLUP <- function(formula, dat, beta, sigma, rho, 
                        sigmaSamplingError = seSigmaClosure(nDomains, nTime)(),
                        nDomains = getNDomains(dat),
                        nTime = getNTime(dat),
                        w0 = w0Matrix(nDomains), 
                        w = w0/rowSums(w0),
                        tol = 1e-3, method = "Nelder-Mead", maxIter = 500) {
  
  modelSpecs <- prepareData(formula, dat, nDomains, nTime, beta, sigma, rho, sigmaSamplingError, w0, w, tol, method, maxIter)
  
  #try-catch handling for parameter estimation
  modelFit <- try(optimizeParameters(modelSpecs), silent = TRUE)
  modelFit <- if (class(modelFit) != "try-error") 
    estimateRE(modelFit) else {
      modelSpecs$beta <- numeric(length(modelSpecs$beta))
      modelSpecs$sigma <- numeric(length(modelSpecs$sigma))
      modelSpecs$rho <- numeric(length(modelSpecs$rho))
      modelSpecs
    }
      
  dat <- dat[order(dat$Domain, dat$Time), ]
  
  output <- list(estimates = data.frame(yHat = modelFit$x %*% modelFit$beta + modelFit$u),
                 beta = modelFit$beta,
                 sigma = modelFit$sigma,
                 rho = modelFit$rho)
  
  class(output) <- c("fitSTREBLUP", "list")
  output
}