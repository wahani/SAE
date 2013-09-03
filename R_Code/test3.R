# rm(list= ls())
# 
# if(.Platform$OS.type == "windows") {
#   require(devtools)
#   if(!require(parallelTools)) install_github(repo="parallelTools", username = "wahani", subdir = "package")
# }
require(devtools)
# install_github(repo="parallelTools", username = "wahani", subdir = "package")
# install_github(repo="spatioTemporalData", username = "wahani", subdir = "package")
# 
# 
require(SAE)
# require(parallelTools)
# 
# set.seed(4)
setup <- simSetupMarhuenda(nDomains=40, nTime=10, sarCorr=0.5, arCorr=0.5, n = 15)
output <- simRunSetup(setup)[[1]]

# output <- setTrueY(output)
# 
# fitFunction <- c("fitEB", "fitSTEBLUP", "fitSTREBLUP")
# 
# simResults <- getSimResults(output, fitFunction)
# save(simResults, file = "Workspaces/tmp1.RData")
# load(file = "Workspaces/simResults1.RData")

dat <- getEvalCrit(simResults[[2]], "calcAABIAS")
dat$cont <- factor(dat$Domain %in% 97:100, labels = c("normal", "cont"))

aggregate(AABIAS ~ cont + model, dat, mean)



plot.simSetup <- function(simSetup) {
  require(ggplot2)
  
  dat1 <- getEvalCrit(simSetup, "calcAABIAS")
  dat2 <- getEvalCrit(simSetup, "calcRRMSE")
  dat3 <- getEvalCrit(simSetup, "calcRBIAS")
  
  list("boxplotRRMSE" = ggplot(dat2, aes(y = calcRRMSE, x = model)) + geom_boxplot() + coord_flip(),
       "boxplotAABIAS" = ggplot(dat1, aes(y = calcAABIAS, x = model)) + geom_boxplot() + coord_flip(),
       "boxplotRBIAS" = ggplot(dat3, aes(y = calcRBIAS, x = model)) + geom_boxplot() + coord_flip()
  )
}

plots <- lapply(simResults, plot)

plots[[1]]$boxplotAABIAS 
plots[[1]]$boxplotRRMSE + scale_y_log10()
plots[[1]]$boxplotRBIAS + ylim(c(-0.2, 0.2))


sapply(split(log(plots[[3]]$boxplotRRMSE$data$calcRRMSE), plots[[3]]$boxplotRRMSE$data$model), 
      mean)
