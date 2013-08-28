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
setup <- simSetupMarhuenda(nDomains=100, nTime=10, sarCorr=0.5, arCorr=0.5, n = 50)
output <- simRunMarhuenda(setup)[[1]]
output <- setTrueY(output)
datList <- output@data

results <- mclapply(datList,
         function(dat) {
           fit <- fitSTREBLUP(dat = dat, 
                             formula = y ~ x, 
                             beta = c(0,1), 
                             sigma = c(1,1), 
                             rho = c(0.5,0.5))
           fit
         }, 
         mc.cores = if (grepl("Windows", Sys.getenv("OS"))) 
           1L else detectCores(),
         mc.preschedule = FALSE)

results1 <- mclapply(datList,
                    function(dat) {
                      fit <- fitSTEBLUP(dat = dat, 
                                         formula = y ~ x, 
                                         beta = c(0,1), 
                                         sigma = c(1,1), 
                                         rho = c(0.5,0.5))
                      fit
                    }, 
                    mc.cores = if (grepl("Windows", Sys.getenv("OS"))) 
                      1L else detectCores(),
                    mc.preschedule = FALSE)



tmp <- lapply(results1, function(tmp)
  data.frame("value" = c(tmp$beta, tmp$sigma, tmp$rho),
             "parameter" = c("beta1", "beta2", "sigma1", "sigma2", "rho1", "rho2"),
             "model" = "ST"))
tmp <- do.call(rbind, tmp)
evalData$model <- "STR" 
evalData <- rbind(evalData, tmp)
require(ggplot2)
resultPlot <- ggplot(evalData, aes(x = value, fill = model, colour = model)) + 
  geom_density(alpha = 0.1) + 
  facet_grid(parameter~., scales="free_y") + 
  labs(title = "Scenario: 100 Domains, 10 Time, TrueParameters: betas: 0, 1, \n rhos: 0.5, 0.5 sigmas: 1, 1, n = 50")

ggsave(filename = "ParameterRobustNonRobustDensities.pdf", width = 10, height = 5)



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
