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
                modelSpecs$k, modelSpecs$Z, modelSpecs$tol, modelSpecs$maxIter)
}

#' @rdname optimizeRE
#' @S3method optimizeRE default
optimizeRE.default <- function(modelSpecs) {
  stop(paste(sub("MS", "", class(modelSpecs)), "is an unknown model specification.")) 
}

#' @rdname optimizeRE
#' @S3method optimizeRE MSRFH
optimizeRE.MSRFH <- function(modelSpecs) {
    
  # Starting optimization/algorithm
  fitre <- try(optimizeRER(modelSpecs$reVar, modelSpecs$vardir, modelSpecs$y, modelSpecs$X, 
                  modelSpecs$beta, modelSpecs$k, modelSpecs$tol, modelSpecs$maxIter))
  
  # Error handling:
  if (inherits(fitre, "try-error")) {
    fitre <- list(x = NA, returnStatus = 2, errorMessage = fitre[1], errorCall = NA)
  } else {
    if(any(is.na(fitre$x))) fitre$x <- matrix(rep_len(0, length(fitre$x)))
    fitre$returnStatus <- if(fitre$returnStatus) 0 else 1
  }
  
  # Adding fitre to modelSpecs
  fitre$call <- match.call()
  class(fitre) <- "fp"
  modelSpecs$fitre <- fitre
  modelSpecs
}

