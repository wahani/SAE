rm(list= ls())

# Daten generieren:

# Installiere Paket
# utils::install.packages(pkgs="../spatioTemporalData/spatioTemporalData_1.0.zip", 
#                         repos = NULL)
# utils::install.packages(pkgs="../spatioTemporalData/spatioTemporalData_1.0.tar.gz", 
#                         repos = NULL)

require(spatioTemporalData)
require(SAE)
set.seed(3)
setup <- simSetupMarhuenda(nDomains=40, nTime=5, sarCorr=0.5, arCorr=0.5, n = 15)
output <- simRunMarhuenda(setup)[[1]]
output <- setTrueY(output)

fitFunction <- c("fitSTEBLUP", "fitSTREBLUP")

simResults <- getSimResults(output, fitFunction)

dataList <- simResults@data

mse <- do.call("rbind", lapply(dataList, 
                        function(dat) {
                          dat <- subset(dat, Time == max(dat$Time))
                          data.frame(mse.ST = calcRRMSE(dat$trueY, dat$yHat.fitSTEBLUP),
                                     mse.STR = calcRRMSE(dat$trueY, dat$yHat.fitSTREBLUP),
                                     mse.direct = calcRRMSE(dat$trueY, dat$y))
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
