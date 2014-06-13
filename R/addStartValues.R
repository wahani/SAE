addStartValues <- function(modelSpecs) {
  # Generic function: computes starting values
  UseMethod("addStartValues")
}

addStartValues.default <- function(modelSpecs) {
  # Default behaviour is not wanted
  stop("This type is not supported!")
}

addStartValues.MSRFH <- function(modelSpecs) {
  # Generic function: computes starting values
  resultFH <- eval(parse(text = paste("sae::eblupFH(modelSpecs$formula,", modelSpecs$vardirName, ",data=modelSpecs$data)")))
  if(any(is.na(resultFH$fit$estcoef))) stop("Starting values are NAs!")
  modelSpecs$beta <- resultFH$fit$estcoef$beta
  modelSpecs$reVar <- resultFH$fit$refvar
  modelSpecs
}