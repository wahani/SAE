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

nDomains <- 100
nTime <- 10
n <- 200

set.seed(11034)
svv00 <- simRunContamination(nDomains=nDomains, nTime=nTime, sarCorr=c(0), arCorr=c(0), n = n,
                              spatialCont = list(sigma = 1, sigmaCont = 9, nDomainsCont = 3),
                              temporalCont = list(sigma = 1, sigmaCont = 9, nDomainsCont = 2),
                              spatioTemporalMessup = TRUE)

svvpp <- simRunContamination(nDomains=nDomains, nTime=nTime, sarCorr=c(0.9), arCorr=c(0.9), n = n,
                             spatialCont = list(sigma = 1, sigmaCont = 9, nDomainsCont = 3),
                             temporalCont = list(sigma = 1, sigmaCont = 9, nDomainsCont = 2),
                             spatioTemporalMessup = TRUE)

s0000 <- simRunContamination(nDomains=nDomains, nTime=nTime, sarCorr=c(0), arCorr=c(0), n = n,
                                   spatialCont = list(sigma = 1, sigmaCont = 1, nDomainsCont = 0),
                                   temporalCont = list(sigma = 1, sigmaCont = 1, nDomainsCont = 0),
                                   spatioTemporalMessup = FALSE)

s00pp <- simRunContamination(nDomains=nDomains, nTime=nTime, sarCorr=c(0.9), arCorr=c(0.9), n = n,
                             spatialCont = list(sigma = 1, sigmaCont = 1, nDomainsCont = 0),
                             temporalCont = list(sigma = 1, sigmaCont = 1, nDomainsCont = 0),
                             spatioTemporalMessup = FALSE)

output <- c(svv00, s00pp, s0000, s00pp, recursive=TRUE)

simResults <- lapply(output, getSimResults, 
                     fitFunction = c("fitEB", "fitSTEBLUP","fitSTREBLUP"),  # "fitEB", "fitSTEBLUP", fitSTREBLUP
                     mc.cores = 27)

save(simResults, file = "Workspaces/simResults6.RData")
#load(file = "Workspaces/simResults1.RData")
