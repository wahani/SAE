rm(list= ls())

# Daten generieren:

# Installiere Paket
# library(spatioTemporalData)
# library(devtools)
# install_github("spatioTemporalData", username = "wahani", subdir = "package")
# utils::install.packages(pkgs="../spatioTemporalData/spatioTemporalData_1.1.4.tar.gz")

"+.simSetup" <- function(x, y) {
  datListX <- slot(x, "data")
  datListY <- slot(y, "data")
  
  datList <- mapply(function(dat1, dat2) data.frame(dat1, dat2$yHat),
                    datListX, datListY, SIMPLIFY = FALSE)
  
  slot(x, "data") <- datList
  x
}

simScenario <- function(sigma, sigmaCont, ...) {
  require(spatioTemporalData)
  set.seed(1)
  svv00 <- simRunContamination(sarCorr=0, arCorr=0, 
                               spatialCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 3),
                               temporalCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 2),
                               spatioTemporalMessup = TRUE, scenarioName = "(v1, v2, 0, 0)",
                               ...)
  
  svvpp <- simRunContamination(sarCorr=0.5, arCorr=0.5, 
                               spatialCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 3),
                               temporalCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 2),
                               spatioTemporalMessup = TRUE, scenarioName = "(v1, v2, 0.5, 0.5)",
                               ...)
  
  s0000 <- simRunContamination(sarCorr=0, arCorr=0,
                               spatialCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 0),
                               temporalCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 0),
                               spatioTemporalMessup = FALSE, scenarioName = "(0, 0, 0, 0)",
                               ...)
  
  s00pp <- simRunContamination(sarCorr=0.5, arCorr=0.5,
                               spatialCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 0),
                               temporalCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 0),
                               spatioTemporalMessup = FALSE, scenarioName = "(0, 0, 0.5, 0.5)", 
                               ...)
  
  output <- c(svv00, svvpp, s0000, s00pp, recursive=TRUE)
  
}


output <- simScenario(n = 200, nDomains = 100, nTime = 10, beta = c(10, 5), sigma = 1, sigmaCont = 50,
                      xdt = spGenerator, seVar = seSigmaClosure)
require(SAE)

simResults <- lapply(output, getSimResults, 
                     fitFunction = c("fitEB", "fitSTEBLUP", "fitSTREBLUP"),
                     mc.cores = 27)

save(simResults, file = "Workspaces/simResults11.RData", compress = TRUE)
#load(file = "Workspaces/simResults1.RData")
