rm(list = ls())
require(SAE)

load("Workspaces//simResults11.RData")

l1 <- simResults
load("Workspaces//simResults12.RData")
simResults[4] <- l1

# # Für 6 und 7 
# load("Workspaces//simResults6.RData")
# simResultSTR <- simResults
# load("Workspaces//simResults7.RData")
# simResultsST <- simResults
# load("Workspaces//simResults8.RData")
# simResultsFH <- simResults
# 
# simResultsFH[[1]]@data[[1]][1:100,]
# 
# addPreds <- function(simResult1, simResult2) {
#   
#   addYHAT <- function(dat1, dat2) {
#     yHatPosition <- grepl("yHat", names(dat2))
#     dat <- data.frame(dat1, "d" = dat2[, yHatPosition])
#     names(dat) <- c(names(dat1), names(dat2)[yHatPosition])
#     dat
#   }
#   
#   datList1 <- simResult1@data
#   datList2 <- simResult2@data
#     
#   simResult1@data <- mapply("addYHAT", datList1, datList2, SIMPLIFY=FALSE)
#   simResult1
# }
# 
# tmp <- mapply(addPreds, simResultSTR, simResultsST, SIMPLIFY=FALSE)
# 
# simResults <- mapply(addPreds, tmp, simResultsFH, SIMPLIFY=FALSE)

plot.simSetup <- function(simSetup, scenario = "") {
  require(ggplot2)
  
  if (scenario == "") scenario <- simSetup@scenarioName
  datList <- lapply(c("calcRRMSE", "calcRBIAS"), getEvalCrit, simResults = simSetup, scenario = scenario)
  
  list("RRMSE" = ggplot(datList[[1]], aes(y = RRMSE, x = model)) + geom_boxplot() + coord_flip(),
       "RBIAS" = ggplot(datList[[2]], aes(y = RBIAS, x = model)) + geom_boxplot() + coord_flip())
}

plotSimResultList <- function(simResultsList, scenarioList = as.list(rep("", length(simResultsList))), critFunctionName) {
  require(ggplot2)
  evalDataList <- mapply(getEvalCrit, simResultsList, scenarioList, critFunctionName = critFunctionName, SIMPLIFY = FALSE)
  evalData <- do.call("rbind", evalDataList)
  
  if (critFunctionName == "calcRBIAS")
    return(ggplot(evalData) + geom_boxplot(aes(x = model, y = RBIAS)) +  
             geom_hline(aes(inetercept = 0), colour = "red") + 
             coord_flip() + facet_grid(Scenario~.))
  
  if (critFunctionName == "calcRRMSE")
    return(ggplot(evalData) + geom_boxplot(aes(x = model, y = RRMSE)) +  
             coord_flip() + facet_grid(Scenario~.))
  
  if (critFunctionName == "calcMSE")
    return(ggplot(evalData) + geom_boxplot(aes(x = model, y = MSE)) +  
             coord_flip() + facet_grid(Scenario~.))
  
  if (critFunctionName == "calcABIAS")
    return(ggplot(evalData) + geom_boxplot(aes(x = model, y = ABIAS)) +  
             coord_flip() + facet_grid(Scenario~.))
  
  if (critFunctionName == "calcBIAS")
    return(ggplot(evalData) + geom_boxplot(aes(x = model, y = BIAS)) +  
             geom_hline(aes(inetercept = 0), colour = "red") + 
             coord_flip() + facet_grid(Scenario~.))
}


plotSimResultList(simResults, critFunctionName = "calcABIAS")




# simResults <- output
# simResults <- lapply(simResults, function(sr) {
#   datList <- sr@data
#   sr@data <- lapply(datList, function(dat) {
#     dat$fit.tmp <- 0
#     dat
#   })
#   sr
# })


# simResults6.RData
# ggRBIAS <- plotSimResultList(simResults, critFunctionName = "calcRBIAS", 
#                              scenarioList = list("(v1, v2, 0, 0)", "(v1, v2, 0.5, 0.5)", "(0, 0, 0, 0)", "(0, 0, 0.5, 0.5)"))
# ggRBIAS + coord_flip(ylim=c(-0.1, 0.1))
# 
# ggRRMSE <- plotSimResultList(simResults, critFunctionName = "calcRRMSE", 
#                              scenarioList = list("(v1, v2, 0, 0)", "(v1, v2, 0.5, 0.5)", "(0, 0, 0, 0)", "(0, 0, 0.5, 0.5)"))
# ggRRMSE + coord_flip(ylim=c(0, 0.5))
# 
# ggMSE <- plotSimResultList(simResults, critFunctionName = "calcMSE", 
#                              scenarioList = list("(v1, v2, 0, 0)", "(v1, v2, 0.5, 0.5)", "(0, 0, 0, 0)", "(0, 0, 0.5, 0.5)"))
# ggMSE + coord_flip(ylim=c(0, 5))
# 
# ggAABIAS <- plotSimResultList(simResults, critFunctionName = "calcAABIAS", 
#                            scenarioList = list("(v1, v2, 0, 0)", "(v1, v2, 0.5, 0.5)", "(0, 0, 0, 0)", "(0, 0, 0.5, 0.5)"))
# ggAABIAS + coord_flip(ylim=c(0, 5))
# 
# ggBIAS <- plotSimResultList(simResults, critFunctionName = "calcBIAS", 
#                               scenarioList = list("(v1, v2, 0, 0)", "(v1, v2, 0.5, 0.5)", "(0, 0, 0, 0)", "(0, 0, 0.5, 0.5)"))
# ggBIAS + coord_flip(ylim=c(0, 5))


# # simResults4.RData
# ggRBIAS <- plotSimResultList(simResults, critFunctionName = "calcRBIAS", 
#                              scenarioList = list("(v1, v2, 0, 0)", "(v1, v2, 0, 0.5)", "(v1, v2, 0.5, 0)", "(v1, v2, 0.5, 0.5)"))
# ggRBIAS + coord_flip(ylim=c(-0.1, 0.1))
# 
# ggRRMSE <- plotSimResultList(simResults, critFunctionName = "calcRRMSE", 
#                              scenarioList = list("(v1, v2, 0, 0)", "(v1, v2, 0, 0.5)", "(v1, v2, 0.5, 0)", "(v1, v2, 0.5, 0.5)"))
# ggRRMSE + coord_flip(ylim=c(0, 5))
# 
# ggMSE <- plotSimResultList(simResults, critFunctionName = "calcMSE", 
#                              scenarioList = list("(v1, v2, 0, 0)", "(v1, v2, 0, 0.5)", "(v1, v2, 0.5, 0)", "(v1, v2, 0.5, 0.5)"))
# ggMSE + coord_flip(ylim=c(0, 5))

# # simResults3.RData
# ggRBIAS <- plotSimResultList(simResults, critFunctionName = "calcRBIAS", scenarioList = list("(0, 0, 0, 0)", "(0, 0, 0, 0.5)", "(0, 0, 0.5, 0)"))
# ggRBIAS + coord_flip(ylim=c(-0.1, 0.1))
# 
# ggRRMSE <- plotSimResultList(simResults, critFunctionName = "calcRRMSE", scenarioList = list("(0, 0, 0, 0)", "(0, 0, 0, 0.5)", "(0, 0, 0.5, 0)"))
# ggRRMSE + coord_flip(ylim=c(0, 5))
# 
# ggRRMSE <- plotSimResultList(simResults, critFunctionName = "calcMSE", scenarioList = list("(0, 0, 0, 0)", "(0, 0, 0, 0.5)", "(0, 0, 0.5, 0)"))
# ggRRMSE + coord_flip(ylim=c(0, 5))

# simResults1.RData + 2
# ggRBIAS <- plotSimResultList(list(simResults[[1]]), critFunctionName = "calcRRMSE", scenarioList = list("(0, 0, p1, p2)"))
# ggRBIAS + coord_flip(ylim=c(-0.1, 0.1))
# 
# ggRRMSE <- plotSimResultList(list(simResults), critFunctionName = "calcRRMSE", scenarioList = list("(0, 0, p1, p2)"))
# ggRRMSE + coord_flip(ylim=c(0, 1))
# 
# ggAABIAS <- plotSimResultList(list(simResults), critFunctionName = "calcAABIAS", scenarioList = list("(0, 0, p1, p2)"))
# ggAABIAS + coord_flip(ylim=c(-0.1, 0.1))


# #simResults5.Rdata
# ggRBIAS <- plotSimResultList(simResults[-c(2,3)], critFunctionName = "calcRBIAS")
# ggRBIAS + coord_flip(ylim=c(-2, 2))
# 
# ggRRMSE <- plotSimResultList(simResults[-c(2,3)], critFunctionName = "calcRRMSE")
# ggRRMSE + coord_flip(ylim=c(0, 40))
# 
# 
# critFunctionName <- "calcRBIAS"
# simResultsList <- simResults[-c(2,3)]
# scenarioList <- list("(v1, v2, 0, 0)", "(v1, v2, p1, p2)")
# 
# evalDataList <- mapply(getEvalCrit, simResultsList, scenarioList, critFunctionName = critFunctionName, SIMPLIFY = FALSE)
# evalData <- do.call("rbind", evalDataList)
# 
# ggplot(evalData) + geom_line(aes(x = Domain, y = RBIAS, linetype = model)) + facet_grid(Scenario~.)


