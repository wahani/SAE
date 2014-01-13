addPrediction <- function(modelSpecs) {
  # Generic function: computes predicted values (BLUP)
  UseMethod("addPrediction")
}

addPrediction.default <- function(modelSpecs) {
  # Generic function: computes predicted values (BLUP)
  modelSpecs$prediction <- modelSpecs$X %*% modelSpecs$beta + modelSpecs$fitre$x
  modelSpecs
}
