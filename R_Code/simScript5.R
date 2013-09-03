rm(list= ls())

# Daten generieren:

# Installiere Paket
# library(spatioTemporalData)
# library(devtools)
# install_github("spatioTemporalData", username = "wahani", subdir = "package")
# utils::install.packages(pkgs="spatioTemporalData_1.0.tar.gz")

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
output <- simRunContamination(nDomains=100, nTime=10, sarCorr=c(0, 0.5), arCorr=c(0, 0.5), n = 200,
                              spatialCont = list(sigma = 1, sigmaCont = 9, nDomainsCont = 3),
                              temporalCont = list(sigma = 1, sigmaCont = 9, nDomainsCont = 2),
                              spatioTemporalMessup = TRUE)

simResults <- lapply(output, getSimResults, 
                     fitFunction = c("fitEB", "fitSTEBLUP","fitSTREBLUP"),  # "fitEB", "fitSTEBLUP", fitSTREBLUP
                     mc.cores = 27)

save(simResults, file = "Workspaces/simResults5.RData")
#load(file = "Workspaces/simResults1.RData")
