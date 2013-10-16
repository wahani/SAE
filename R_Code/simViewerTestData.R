load("H:/R/Projekte/SAE/Workspaces/simResults12.RData")

r1 <- simResults[[1]]

counterClosure <- function() {
  index <- 0
  function() {
    index <<- index + 1
    return(index)
  }
}
counter <- counterClosure()
  
tmp <- do.call(rbind, lapply(r1@data, 
       function(dat) {
         dat$r <- counter()
         dat
       }))

tmp <- tmp[tmp$Time == 10, c("trueY", "Domain", "r", "y", "yHat.fitEB", "yHat.fitSTEBLUP", "yHat.fitSTREBLUP")]
names(tmp) <- c("y", "domain", "r", "est.Direct", "est.FH", "est.STEBLUP", "est.STREBLUP")
rownames(tmp) <- NULL

simData <- list(estData = tmp)
save(simData, file = "Workspaces/simViewerTestData.RData")