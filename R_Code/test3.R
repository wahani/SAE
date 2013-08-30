rm(list= ls())

if(.Platform$OS.type == "windows") {
  require(devtools)
  if(!require(parallelTools)) install_github(repo="parallelTools", username = "wahani", subdir = "package")
}

"+.simSetup" <- function(x, y) {
  datListX <- slot(x, "data")
  datListY <- slot(y, "data")
  
  datList <- mapply(function(dat1, dat2) data.frame(dat1, dat2$yHat),
                    datListX, datListY, SIMPLIFY = FALSE)
  
  slot(x, "data") <- datList
  x
}

require(SAE)
require(parallelTools)

set.seed(4)
setup <- simSetupMarhuenda(nDomains=40, nTime=10, sarCorr=0.5, arCorr=0.5, n = 20)
output <- simRunMarhuenda(setup)[[1]]
output <- setTrueY(output)

fitFunction <- c("fitEB", "fitSTEBLUP", "fitSTREBLUP")

simResults <- getSimResults(output, fitFunction)
save(simResults, file = "Workspaces/tmp1.RData")


getEvalCrit <- function(simResults, critFunctionName = "calcRRMSE") {
  require(reshape2)
  critFunction <- match.fun(critFunctionName)
  data <- subset(do.call("rbind", simResults@data), Time == simResults@nTime)
  dataList <- split(data, as.factor(data$Domain))
  
  result <- do.call("rbind", lapply(dataList, 
                                    function(dat) {
                                      data.frame(STFH = critFunction(dat$trueY, dat$yHat.fitSTEBLUP),
                                                 STRFH = critFunction(dat$trueY, dat$yHat.fitSTREBLUP),
                                                 Direct = critFunction(dat$trueY, dat$y),
                                                 FH = critFunction(dat$trueY, dat$yHat.fitEB))
                                    }
  ))
  
  result <- melt(result, variable.name="model", value.name = critFunctionName)
  result$Domain <- rep(1:simResults@nDomains, 4)
  result
}

calcAABIAS <- function(trueValues, estimates) {
  estimates[estimates == 0] <- NA
  mean(abs(trueValues-estimates), na.rm=T)
}


dat1 <- getEvalCrit(simResults, "calcAABIAS")
dat2 <- getEvalCrit(simResults)
require(ggplot2)
ggplot(dat2, aes(y = calcRRMSE, x = model)) + geom_boxplot() + coord_flip() #  + ylim(c(0, 5))
ggplot(dat, aes(y = calcAABIAS, x = model)) + geom_boxplot() + coord_flip()












