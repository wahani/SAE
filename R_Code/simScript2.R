rm(list= ls())

"+.simSetup" <- function(x, y) {
  datListX <- slot(x, "data")
  datListY <- slot(y, "data")
  
  datList <- mapply(function(dat1, dat2) data.frame(dat1, dat2$yHat),
                    datListX, datListY, SIMPLIFY = FALSE)
  
  slot(x, "data") <- datList
  x
}

require(SAE)

set.seed(3)
setup <- simSetupMarhuenda(nDomains=100, nTime=10, sarCorr=0.5, arCorr=0.5, n = 200)
output <- simRunMarhuenda(setup)[[1]]
output <- setTrueY(output)

fitFunction <- c("fitEB", "fitSTEBLUP", "fitSTREBLUP")

simResults <- getSimResults(output, fitFunction, 27)

save(simResults, file = "Workspaces/simResults2.RData")
#load(file = "Workspaces/simResults1.RData")
