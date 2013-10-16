rm(list = ls())
require(SAE)

load("Workspaces//simResults11.RData")
l1 <- simResults
load("Workspaces//simResults12.RData")
simResults[4] <- l1

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

setting <- "R = 200; D = 100, T = 10"
width = 7
height = 0.75 * width

# Absolute Bias
ggABIAS <- plotSimResultList(simResults, critFunctionName = "calcABIAS")  + 
  labs(title = paste("Absolut BIAS;", setting), x = NULL) + theme_bw()
ggsave(filename="Grafiken/ABIAS.pdf", plot = ggABIAS, width = width, height = height)

ggABIAS_zoom <- ggABIAS + coord_flip(ylim = c(0.5, 1))
ggsave(filename="Grafiken/ABIAS_zoom.pdf", plot = ggABIAS_zoom , width = width, height = height)


# Relative Bias
ggRBIAS <- plotSimResultList(simResults, critFunctionName = "calcRBIAS") + labs(x = NULL) + theme_bw()
ggsave(filename="Grafiken/RBIAS.pdf", plot = ggRBIAS, width = width, height = height)


# Relative Root MSE
ggRRMSE <- plotSimResultList(simResults, critFunctionName = "calcRRMSE") + labs(x = NULL) + theme_bw()
ggsave(filename="Grafiken/RRMSE.pdf", plot = ggRRMSE, width = width, height = height)

ggRRMSE_zoom <- ggRRMSE + coord_flip(y = c(0.02, 0.06))
ggsave(filename="Grafiken/RRMSE_zoom.pdf", plot = ggRRMSE_zoom, width = width, height = height)

# Plot for Data
df1 <- simResults[[1]]@data[[3]]
df1$Scenario <- simResults[[1]]@scenarioName
df2 <- simResults[[2]]@data[[1]]
df2$Scenario <- simResults[[2]]@scenarioName

df <- rbind(df1, df2)
df$Outliers <- factor(ifelse(df$Scenario == simResults[[1]]@scenarioName & df$Domain %in% 96:100, 1, 2), 
                      labels = c("outliers", "no outliers"))

ggData <- ggplot(df) + geom_point(aes(x = x, y = y, colour = Outliers)) + facet_grid(.~Scenario) + 
  scale_colour_manual(values = c("outliers" = "dodgerblue", "no outliers" = "black"), 
                      guide = "none")+ theme_bw()
ggsave(filename="Grafiken/data.pdf", plot = ggData, width = width, height = height)
##
library(xtable)

s <- lapply(simResults, function(sR) {
  dat <- do.call("rbind", sR@data)
  s <- as.matrix(summary(dat$y), ncol = 1)
  colnames(s) <- sR@scenarioName
  s
})

outTable <- do.call("cbind", s)[, c(2, 3, 1, 4)]

print(xtable(outTable, align = "lrrrr"), file = "summarTable.tex")






