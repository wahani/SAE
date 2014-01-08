prepareData <- function(formula, dat, nDomains, nTime, beta, sigma, rho, sigmaSamplingError, w0, w, tol, method, maxIter, k = 1.345) {
  
  dat <- dat[order(dat$Domain, dat$Time), ]
  XY <- makeXY(formula, dat)
  
  modelSpecs <- list(y = XY$y,
                     X = XY$x,
                     nDomains = nDomains,
                     nTime = nTime,
                     
                     # Proximity Matrix
                     w0 = w0,
                     w = w,
                     
                     # Z Matrix for Random Effects
                     Z = reZ(getNDomains(dat), getNTime(dat)),
                     Z1 = reZ1(getNDomains(dat), getNTime(dat)),
                     
                     # Initial Parameters
                     beta = beta + c(0.01, 0.01),
                     sigma = sigma + c(0.01, 0.01),
                     rho = rho + c(0.01, 0.01),
                     
                     # True Variance Parameters for Sampling Error
                     sigmaSamplingError = sigmaSamplingError,
                     
                     # Specs for estimation
                     tol = tol,
                     psiFunction = psiOne,
                     k = k,
                     K = 2 * pnorm(k) - 1 - 2 * k * dnorm(k) + 2 * (k^2) * (1-pnorm(k)),
                     method = method,
                     maxIter = maxIter,
                     consoleOutput = TRUE,
                     
                     # "Empty" Vector of Random Effects Estimates
                     u = numeric(nDomains * nTime))
  class(modelSpecs) <- c("MSSTRFH")
  return(modelSpecs)
}

#' genOptsRobust
#' 
#' Constructs a list with necessary elements for \code{fitfh} - only a convinience function.
#' 
#' @param k 
#' @param K 
#' 
#' @return A named list with \code{k} and \code{K} as elements.
#' @export
genOptsRobust <- function(k = 1.345, K = 2 * pnorm(k) - 1 - 2 * k * dnorm(k) + 2 * (k^2) * (1-pnorm(k))) {
  list(k = k, K = K)
}

#' genOptsOptim
#' 
#' Constructs a list with necessary elements for \code{fitfh} - only a convinience function.
#' 
#' @param tol tolerance criterion for algorithms
#' @param maxIter maximal number of iterations, used as a parameter in each optimization step (not cummulative)
#' @param progress logical - print progress to the console
#' 
#' @return A named list with \code{tol}, \code{maxIter} and \code{progress} as elements.
#' @export
genOptsOptim <- function(tol = 1e-6, maxIter = 100, progress = FALSE) {
  list(tol = tol, maxIter = maxIter, progress = progress)
}

genModelSpecs <- function(optsRobust = genOptsRobust(), optsOptim = genOptsOptim(), type) {
  out <- new.env()
  mapply("assign", names(optsRobust), optsRobust, MoreArgs = list(envir = out))
  mapply("assign", names(optsOptim), optsOptim, MoreArgs = list(envir = out))
  out$psiFunction <- psiOne
  class(out) <- c(class(out), paste("MS", type, sep = ""))
  out
}
