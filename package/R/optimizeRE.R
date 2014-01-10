#' optimizeRE - generic funciton
#' 
#' Function for estimation of random effects.
#' 
#' @param modelSpecs modelSpecs object, initialized with \code{mS()}
#' 
#' @rdname optimizeRE
#' @export optimizeRE 
optimizeRE <- function(modelSpecs) {
  UseMethod("optimizeRE")
}

#' @rdname optimizeRE
#' @S3method optimizeRE MSSTRFH
optimizeRE.MSSTRFH <- function(modelSpecs) {
  optimizeRESTR(modelSpecs$sigma, modelSpecs$rho, 
                modelSpecs$y, modelSpecs$X, modelSpecs$Z1, modelSpecs$sigmaSamplingError, 
                modelSpecs$w, modelSpecs$beta, modelSpecs$nDomains, modelSpecs$nTime, 
                modelSpecs$K, modelSpecs$Z, modelSpecs$tol, modelSpecs$maxIter)
}

#' @rdname optimizeRE
#' @S3method optimizeRE default
optimizeRE.default <- function(modelSpecs) {
  stop(paste(sub("MS", "", class(modelSpecs)), "is an unknown model specification.")) 
}

#' @rdname optimizeRE
#' @S3method optimizeRE MSRFH
optimizeRE.MSRFH <- function(modelSpecs) {
  optimizeRER(modelSpecs$reVar, modelSpecs$vardir, modelSpecs$y, modelSpecs$X, 
              modelSpecs$beta, modelSpecs$K, modelSpecs$tol, modelSpecs$maxIter)
}

