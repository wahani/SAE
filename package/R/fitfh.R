fitfh <- function(formula, vardir, idName, data, optsRobust = genOptsRobust(), optsOptim = genOptsOptim(), type) {
  modelSpecs <- genModelSpecs(optsRobust, optsOptim, type)
  modelSpecs <- addModelFrame(modelSpecs, formula, vardir, idName, data)
  modelSpecs <- addStartValues(modelSpecs)
  modelSpecs <- optimizeParam(modelSpecs)
  modelSpecs <- optimizeRE(modelSpecs)
  
  out <- list(prediction = modelSpecs$prediction, fitparam = modelSpecs$fitparam, 
              fitre = modelSpecs$fitre)
  
  out
}