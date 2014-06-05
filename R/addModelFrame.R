addModelFrame <- function(modelSpecs, formula, vardir, idName, data) {
  # Generic function
  UseMethod("addModelFrame")
}

addModelFrame.default <- function(modelSpecs, formula, vardir, idName, data) {
  # Default behaviour is not wanted
  stop("This type is not supported!")
}

addModelFrame.MSRFH <- function(modelSpecs, formula, vardir, idName, data) {  
  # Method for the robust FH
  data <- data[order(data[[idName]]), ]
  XY <- makeXY(formula, data)
  modelSpecs$y <- XY$y
  modelSpecs$X <- XY$x
  modelSpecs$vardir <- data[[vardir]]
  modelSpecs$vardirName <- vardir
  modelSpecs$data <- data
  modelSpecs$formula <- formula
  modelSpecs
}