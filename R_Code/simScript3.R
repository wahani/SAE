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
setup <- simSetupMarhuenda(nDomains=100, nTime=10, sarCorr=c(0, 0.5), arCorr=c(0, 0.5), n = 200)
output <- simRunMarhuenda(setup)
output <- lapply(output, setTrueY)[-4]

simResults <- lapply(output, getSimResults, 
                     fitFunction = c("fitSTREBLUP"), 
                     mc.cores = 25)

save(simResults, file = "Workspaces/simResults3.RData")
#load(file = "Workspaces/simResults1.RData")
