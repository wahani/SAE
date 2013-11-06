rm(list= ls())

# Daten generieren:

# Installiere Paket
# library(spatioTemporalData)
# library(devtools)
# install_github("spatioTemporalData", username = "wahani", subdir = "package")
# utils::install.packages(pkgs="../spatioTemporalData/spatioTemporalData_1.1.4.tar.gz")

# Helper
"+.simSetup" <- function(x, y) {
  datListX <- slot(x, "data")
  datListY <- slot(y, "data")
  
  datList <- mapply(function(dat1, dat2) data.frame(dat1, dat2$yHat),
                    datListX, datListY, SIMPLIFY = FALSE)
  
  slot(x, "data") <- datList
  x
}

counterClosure <- function() {
  i <- 0
  function() {
    i <<- i + 1
    i
  }
}

# Functions for simulation
fitSimMapReduce <- function(fitFunction, simSetup, ind) {
  # Helper
  addPrediction <- function(dat, prediction) {
    if (max(dat$Domain) == nrow(prediction)) {
      dat$yHat <- NA
      dat[dat$Time == max(dat$Time), "yHat"] <- prediction 
      return(dat)
    } else return(data.frame(dat, yHat = prediction))
  }
  
  ##############################################################################
  startBeta <- simSetup@beta
  startSigma <- c(simSetup@sigma1, simSetup@sigma2)
  startRho <- c(simSetup@sarCorr, simSetup@arCorr)
  
  datList <- slot(simSetup, "data")[ind]
  
  predictionList <- lapply(datList,
                           function(dat) {
                             fit <- fitFunction(dat = dat, 
                                                formula = y ~ x, 
                                                beta = startBeta, 
                                                sigma = startSigma, 
                                                rho = startRho,
                                                sigmaSamplingError = simSetup@sigmaSE,
                                                w0 = simSetup@neighbourHood)
                             fit$estimates
                           })
  
  predictionList[grepl("try-error", sapply(predictionList, class))] <- 
    list(rep(NA, simSetup@nDomains * simSetup@nTime))
  
  slot(simSetup, "data") <- mapply("addPrediction", datList, predictionList, SIMPLIFY = FALSE)
  
  dataForSave <- do.call(rbind, slot(simSetup, "data"))
  dataForSave
}


################################################################################

simSetup <- simScenario(n = 4, nDomains = 50, nTime = 10, beta = c(100, 1), sigma = 1, sigmaCont = 40,
                      xdt = spGenerator, seVar = seSigmaClosure)[[1]]

counter <- counterClosure()

slot(simSetup, "data") <- lapply(simSetup@data, 
       function(dat) {
         dat$ID <- counter()
         dat$scenario <- simSetup@scenarioName
         dat
       })


system.time(simResults <- fitSimMapReduce(fitSTREBLUP, simSetup, ind = 1))

simResults




