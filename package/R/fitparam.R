fitparam <- function(modelSpecs) {
  # Generic function: fit model parameters
  UseMethod("fitparam")
}

fitparam.default <- function(modelSpecs) {
  # Generic function: fit model parameters
  stop("This type is not supported!")
}

fitparam.MSRFH <- function(modelSpecs) {
  # Generic function: fit model parameters
  modelSpecs <- optimizeBeta(modelSpecs)
  modelSpecs
}
