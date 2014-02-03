fitfh <- function(formula, vardir, idName, data, optsRobust = genOptsRobust(), 
                  optsOptim = genOptsOptim(), type = "RFH") {
  # Check Input:
  availableTypes <- c("RFH", "tmp")
  if(!(type %in% availableTypes)) 
    stop("The type '", type, "' is not supported. Choose one in: ", 
         paste(availableTypes, collapse = ", "))
  
  # Fit Fay-Herriot model
  modelSpecs <- genModelSpecs(optsRobust, optsOptim, type)
  modelSpecs <- addModelFrame(modelSpecs, formula, vardir, idName, data)
  modelSpecs <- addStartValues(modelSpecs)
  modelSpecs <- optimizeParam(modelSpecs)
  modelSpecs <- optimizeRE(modelSpecs)
  modelSpecs <- addPrediction(modelSpecs)
  out <- list(call = match.call(), prediction = modelSpecs$prediction, fitparam = modelSpecs$fitparam, 
              fitre = modelSpecs$fitre)
  class(out) <- type
  out
}