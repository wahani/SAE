# Simulation Scenario Skript

rm(list= ls())
utils::install.packages(pkgs="spatioTemporalData_1.0.tar.gz", 
                        repos = NULL)
utils::install.packages(pkgs="SAE_0.1.tar.gz", 
                        repos = NULL)

setTrueY <- function(simSetup) {
  # funciton-definition
  trueY <- function(dat, sigmaE) {
    dat$trueY <- dat$y - sigmaE 
    dat
  }
  
  #
  sigmaE <- as.data.frame(slot(simSetup, "sigma"))
  dat <- slot(simSetup, "data")
  dataList <- mapply("trueY", dat, sigmaE, SIMPLIFY = FALSE)
  slot(simSetup, "data") <- dataList
  simSetup
}

"+.simSetup" <- function(x, y) {
  datListX <- slot(x, "data")
  datListY <- slot(y, "data")
  
  datList <- mapply(function(dat1, dat2) data.frame(dat1, dat2$yHat),
                    datListX, datListY, SIMPLIFY = FALSE)
  
  slot(x, "data") <- datList
  x
}

require(spatioTemporalData)
require(SAE)
set.seed(3)
setup <- simSetupMarhuenda(nDomains=100, nTime=20, sarCorr=0.5, arCorr=0.5, n = 200)
output <- simRunMarhuenda(setup)[[1]]
output <- setTrueY(output)

fitFunction <- c("fitSTEBLUP", "fitSTREBLUP")

simResults <- getSimResults(output, fitFunction)

save(simResults, file = "Workspaces/simResults1.RData")
dataList <- simResults@data

