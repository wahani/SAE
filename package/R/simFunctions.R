# Functions for simulation
fitSimSetup <- function(fitFunction, simSetup) {
  # function definition
  addPrediction <- function(dat, prediction) data.frame(dat, yHat = prediction)
  
  #
  startBeta <- simSetup@beta
  startSigma <- c(simSetup@sigma1, simSetup@sigma2)
  startRho <- c(simSetup@sarCorr, simSetup@arCorr)
  
  datList <- slot(simSetup, "data")
  
  predictionList <- mclapply(datList,
                           function(dat) {
                             fit <- fitFunction(dat = dat, 
                                                formula = y ~ x, 
                                                beta = startBeta, 
                                                sigma = startSigma, 
                                                rho = startRho)
                             fit$estimates
                           }, 
                             mc.cores = if (grepl("Windows", Sys.getenv("OS"))) 
                               1L else detectCores(),
                             mc.preschedule = FALSE)
  
  predictionList[grepl("try-error", sapply(predictionList, class))] <- 
    list(rep(NA, simSetup@nDomains * simSetup@nTime))
  
  slot(simSetup, "data") <- mapply("addPrediction", datList, predictionList, SIMPLIFY = FALSE)
  simSetup
}

getSimResults <- function(simSetup, fitFunction) {
  funList <- lapply(fitFunction, match.fun)
  simResults <- lapply(funList, fitSimSetup, simSetup = simSetup)
  combineSimResults(simResults, fitFunction)
}
