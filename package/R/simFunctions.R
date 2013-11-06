simScenario <- function(sigma, sigmaCont, ...) {
  require(spatioTemporalData)
  set.seed(1)
  
  svvpp <- simRunContamination(sarCorr=0.5, arCorr=0.5, 
                               spatialCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 3),
                               temporalCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 2),
                               spatioTemporalMessup = TRUE, scenarioName = "(v1, v2, 0.5, 0.5)",
                               ...)
  
  output <- c(svvpp, recursive=TRUE)
  
}

# Functions for simulation
fitSimSetup <- function(fitFunction, simSetup, mc.cores = detectCores()) {
  # function definition
  addPrediction <- function(dat, prediction) {
    if (max(dat$Domain) == nrow(prediction)) {
      dat$yHat <- NA
      dat[dat$Time == max(dat$Time), "yHat"] <- prediction 
      return(dat)
    } else return(data.frame(dat, yHat = prediction))
  }
  
  #
  startBeta <- simSetup@beta
  startSigma <- c(simSetup@sigma1, simSetup@sigma2)
  startRho <- c(simSetup@sarCorr, simSetup@arCorr)
  
  datList <- slot(simSetup, "data")
  
  predictionList <- parallelTools::mclapply(datList,
                                            function(dat) {
                                              fit <- fitFunction(dat = dat, 
                                                                 formula = y ~ x, 
                                                                 beta = startBeta, 
                                                                 sigma = startSigma, 
                                                                 rho = startRho,
                                                                 sigmaSamplingError = simSetup@sigmaSE,
                                                                 w0 = simSetup@neighbourHood)
                                              fit$estimates
                                            }, 
                                            mc.cores = mc.cores,
                                            mc.preschedule = FALSE, 
                                            packageToLoad = c("spatioTemporalData", "SAE"))
  
  predictionList[grepl("try-error", sapply(predictionList, class))] <- 
    list(rep(NA, simSetup@nDomains * simSetup@nTime))
  
  slot(simSetup, "data") <- mapply("addPrediction", datList, predictionList, SIMPLIFY = FALSE)
  simSetup
}

getSimResults <- function(simSetup, fitFunction, mc.cores = detectCores()) {
  funList <- lapply(fitFunction, match.fun)
  simResults <- lapply(funList, fitSimSetup, simSetup = simSetup, mc.cores = mc.cores)
  combineSimResults(simResults, fitFunction)
}
