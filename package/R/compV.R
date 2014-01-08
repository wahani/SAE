compV <- function(modelSpecs) {
  # Generic function: computes the variance covariance matrix V
  UseMethod("compV")
}

compV.default <- function(modelSpecs) {
  # Default behaviour is not wanted
  stop("This type is not supported!")
}

compV.MSSTRFH <- function(modelSpecs) {
  # Default behaviour is not wanted
  matVinv(W=modelSpecs$w, rho1=modelSpecs$rho[1], sigma1=modelSpecs$sigma[1],
          rho2 = modelSpecs$rho[2], sigma2 = modelSpecs$sigma[2], Z1=modelSpecs$Z1,
          modelSpecs$sigmaSamplingError)
}

compV.MSRFH <- function(modelSpecs) {
  # Default behaviour is not wanted
  V <- diag(modelSpecs$vardir + modelSpecs$reVar)
  Vinv <- diag(1 / (modelSpecs$vardir + modelSpecs$reVar))
  list(V = V, Vinv = Vinv)
}