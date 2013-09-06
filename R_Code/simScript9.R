rm(list= ls())

# Daten generieren:

# Installiere Paket
# library(spatioTemporalData)
# library(devtools)
# install_github("spatioTemporalData", username = "wahani", subdir = "package")
#utils::install.packages(pkgs="../spatioTemporalData/spatioTemporalData_1.1.3.tar.gz")

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
beta <- c(20, 1)
sigma <- 1
sigmaCont <- 10

set.seed(1)
svv00 <- simRunContamination(nDomains=nDomains, nTime=nTime, sarCorr=0, arCorr=0, n = n,
                             spatialCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 3),
                             temporalCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 2),
                             beta = beta,
                             spatioTemporalMessup = TRUE)

svvpp <- simRunContamination(nDomains=nDomains, nTime=nTime, sarCorr=0.5, arCorr=0.5, n = n,
                             spatialCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 3),
                             temporalCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 2),
                             spatioTemporalMessup = TRUE)

s0000 <- simRunContamination(nDomains=nDomains, nTime=nTime, sarCorr=0, arCorr=0, n = n,
                             spatialCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 0),
                             temporalCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 0),
                             spatioTemporalMessup = FALSE)

s00pp <- simRunContamination(nDomains=nDomains, nTime=nTime, sarCorr=0.5, arCorr=0.5, n = n,
                             spatialCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 0),
                             temporalCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 0),
                             spatioTemporalMessup = FALSE)

output <- c(svv00, s00pp, s0000, s00pp, recursive=TRUE)

simResults <- lapply(output, getSimResults, 
                     fitFunction = c("fitEB", "fitSTEBLUP", "fitSTREBLUP"),
                     mc.cores = 27)

save(simResults, file = "Workspaces/simResults9.RData", compress = TRUE)
#load(file = "Workspaces/simResults1.RData")
