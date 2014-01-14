fitfh <- function(formula, vardir, idName, data, optsRobust = genOptsRobust(), 
                  optsOptim = genOptsOptim(), type = "RFH") {
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