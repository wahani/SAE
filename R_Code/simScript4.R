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
output <- simRunContamination(nDomains=100, nTime=10, sarCorr=c(0, 0.5), arCorr=c(0, 0.5), n = 200)
output <- lapply(output, setTrueY)

# dat <- output[[1]]@data[[1]]
# fitEB(formula=y~x, dat=dat, beta=c(0,1), sigma=c(1,1), rho=c(0.5, 0.5), sigmaSamplingError=output[[1]]@sigmaSE, w0=output[[1]]@neighbourHood)



simResults <- lapply(output, getSimResults, 
                     fitFunction = c("fitEB", "fitSTEBLUP","fitSTREBLUP"),  # "fitEB", "fitSTEBLUP", fitSTREBLUP
                     mc.cores = 27)

save(simResults, file = "Workspaces/simResults4.RData")
#load(file = "Workspaces/simResults1.RData")
