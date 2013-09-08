rm(list= ls())
# 
# if(.Platform$OS.type == "windows") {
#   require(devtools)
#   if(!require(parallelTools)) install_github(repo="parallelTools", username = "wahani", subdir = "package")
# }
require(devtools)
# install_github(repo="parallelTools", username = "wahani", subdir = "package")
install_github(repo="spatioTemporalData", username = "wahani", subdir = "package")
# 
#install.packages("../spatioTemporalData/spatioTemporalData_1.1.1.tar.gz")
#

"+.simSetup" <- function(x, y) {
  datListX <- slot(x, "data")
  datListY <- slot(y, "data")
  
  datList <- mapply(function(dat1, dat2) data.frame(dat1, dat2$yHat),
                    datListX, datListY, SIMPLIFY = FALSE)
  
  slot(x, "data") <- datList
  x
}

require(SAE)
# require(parallelTools)
# 
set.seed(3)
output <- simRunContamination(nDomains=50, nTime=10, sarCorr=0.5, arCorr=0.5, n = 50,
                             spatialCont = list(sigma = 1, sigmaCont = 9, nDomainsCont = 2),
                             temporalCont = list(sigma = 1, sigmaCont = 9, nDomainsCont = 2),
                             spatioTemporalMessup = TRUE)

simResults <- lapply(output, getSimResults, 
                     fitFunction = c("fitEB", "fitSTEBLUP", "fitSTREBLUP"), # , "fitSTREBLUP" fitSTEBLUP
                     mc.cores = 27)

plot.simSetup <- function(simSetup) {
  require(ggplot2)
  
  datList <- lapply(c("calcRRMSE", "calcRBIAS"), getEvalCrit, simResults = simSetup)
  
  list("RRMSE" = ggplot(datList[[1]], aes(y = RRMSE, x = model)) + geom_boxplot() + coord_flip(),
       "RBIAS" = ggplot(datList[[2]], aes(y = RBIAS, x = model)) + geom_boxplot() + coord_flip())
}

ggplot(getEvalCrit(simResults[[1]], critFunctionName="calcRBIAS")) + 
  geom_line(aes(x = Domain, y = RBIAS))

plots <- lapply(simResults, plot)

plots[[1]]$RRMSE + coord_flip(ylim=c(0,20))
plots[[1]]$RBIAS + coord_flip(ylim=c(-2,2))

sapply(split(log(plots[[3]]$boxplotRRMSE$data$calcRRMSE), plots[[3]]$boxplotRRMSE$data$model), 
      mean)
