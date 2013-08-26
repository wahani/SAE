rm(list= ls())

# Daten generieren:

# Installiere Paket
# utils::install.packages(pkgs="../spatioTemporalData/spatioTemporalData_1.0.zip", 
#                         repos = NULL)
# utils::install.packages(pkgs="../spatioTemporalData/spatioTemporalData_1.0.tar.gz", 
#                         repos = NULL)

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
setup <- simSetupMarhuenda(nDomains=30, nTime=20, sarCorr=0.5, arCorr=0.5, n = 50)
output <- simRunMarhuenda(setup)[[1]]
output <- setTrueY(output)

fitFunction <- c("fitSTEBLUP", "fitSTREBLUP")

system.time(simResults <- getSimResults(output, fitFunction))

dataList <- simResults@data

mse <- do.call("rbind", lapply(dataList, 
                        function(dat) {
                          dat <- subset(dat, Time == max(dat$Time))
                          data.frame(mse.ST = mean((dat$trueY - dat$yHat.fitSTEBLUP)^2),
                                     mse.STR = mean((dat$trueY - dat$yHat.fitSTREBLUP)^2),
                                     mse.direct = mean((dat$y - dat$yHat.fitSTREBLUP)^2))
                        }
))

lapply(mse, mean, na.rm = T)


# Prepare Data
# system.time(tmp <- optimizeParameters(modelSpecs))

Rprof(tmp <- tempfile())
out <- fitSTREBLUP(y~x, dat, c(0,1), c(1,1), c(0.5,0.5))
Rprof()
profileSummary<- summaryRprof(tmp)
unlink(tmp)

profileSummary$by.total

fit <- fitSTREBLUP(y~x, dat, c(0,1), c(1,1), c(0.5,0.5))
# modelFit <- fitSTEBLUP(y~x, dat, c(0,1), c(1,1), c(0.5,0.5))

fit <- fitEB(y~x, dat, c(0,1), c(1,1), c(0.5,0.5))
