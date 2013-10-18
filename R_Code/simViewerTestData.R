load("H:/R/Projekte/SAE/Workspaces/simResults12.RData")

prepareData <- function(sr) {
  counterClosure <- function() {
    index <- 0
    function() {
      index <<- index + 1
      return(index)
    }
  }
  counter <- counterClosure()
  
  tmp <- do.call(rbind, lapply(sr@data, 
                               function(dat) {
                                 dat$r <- counter()
                                 dat
                               }))
  
  tmp <- tmp[tmp$Time == 10, c("trueY", "Domain", "r", "y", "yHat.fitEB", 
                               "yHat.fitSTEBLUP", "yHat.fitSTREBLUP")]
  names(tmp) <- c("y", "domain", "r", "est.Direct", "est.FH", "est.STEBLUP", 
                  "est.STREBLUP")
  rownames(tmp) <- NULL
  
  tmp$scenario <- sr@scenarioName
  tmp
} 

outData <- do.call(rbind, lapply(simResults, prepareData))
outData$domain <- factor(outData$domain)
outData$scenario <- factor(outData$scenario)

levels(outData$scenario)
str(outData)

simData <- list(estData = outData)
save(simData, file = "Workspaces/simViewerTestData.RData")