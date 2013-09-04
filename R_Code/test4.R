rm(list = ls())
require(SAE)

load("Workspaces//simResults5.RData")

plot.simSetup <- function(simSetup, scenario = "") {
  require(ggplot2)
  
  datList <- lapply(c("calcRRMSE", "calcRBIAS"), getEvalCrit, simResults = simSetup, scenario = scenario)
  
  list("RRMSE" = ggplot(datList[[1]], aes(y = RRMSE, x = model)) + geom_boxplot() + coord_flip(),
       "RBIAS" = ggplot(datList[[2]], aes(y = RBIAS, x = model)) + geom_boxplot() + coord_flip())
}

plotSimResultList <- function(simResultsList, scenarioList = list("(v1, v2, 0, 0)", "(v1, v2, p1, p2)"), critFunctionName) {
  evalDataList <- mapply(getEvalCrit, simResultsList, scenarioList, critFunctionName = critFunctionName, SIMPLIFY = FALSE)
  evalData <- do.call("rbind", evalDataList)
  
  if (critFunctionName == "calcRBIAS")
    return(ggplot(evalData) + geom_boxplot(aes(x = model, y = RBIAS)) +  
             geom_hline(aes(inetercept = 0), colour = "red") + 
             coord_flip() + facet_grid(Scenario~.))
  
  if (critFunctionName == "calcRRMSE")
    return(ggplot(evalData) + geom_boxplot(aes(x = model, y = RRMSE)) +  
             coord_flip() + facet_grid(Scenario~.))
}

ggRBIAS <- plotSimResultList(simResults[-c(2,3)], critFunctionName = "calcRBIAS")
ggRBIAS + coord_flip(ylim=c(-2, 2))

ggRRMSE <- plotSimResultList(simResults[-c(2,3)], critFunctionName = "calcRRMSE")
ggRRMSE + coord_flip(ylim=c(0, 40))


critFunctionName <- "calcRBIAS"
simResultsList <- simResults[-c(2,3)]
scenarioList <- list("(v1, v2, 0, 0)", "(v1, v2, p1, p2)")

evalDataList <- mapply(getEvalCrit, simResultsList, scenarioList, critFunctionName = critFunctionName, SIMPLIFY = FALSE)
evalData <- do.call("rbind", evalDataList)

ggplot(evalData) + geom_line(aes(x = Domain, y = RBIAS, linetype = model)) + facet_grid(Scenario~.)


