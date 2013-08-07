prepareData <- function(formula, dat, nDomains, nTime, beta, sigma, rho, sigmaSamplingError, w0, w, tol, method, maxIter) {
  
  dat <- dat[order(dat$Domain, dat$Time), ]
  XY <- makeXY(formula, dat)
  
  modelSpecs <- list(y = XY$y,
                     x = XY$x,
                     nDomains = nDomains,
                     nTime = nTime,
                     
                     w0 = w0,
                     w = w,
                     
                     Z = reZ(getNDomains(dat), getNTime(dat)),
                     Z1 = reZ1(getNDomains(dat), getNTime(dat)),
                     
                     beta = beta,
                     sigma = sigma,
                     rho = rho,
                     
                     sigmaSamplingError = sigmaSamplingError,
                     
                     tol = tol,
                     psiFunction = psiOne,
                     K = 0.71,
                     method = method,
                     maxIter = maxIter)
  return(modelSpecs)
}