library(saedevel)
library(sae)
library(dplyr)
detach("package:plyr", unload=TRUE)
library(reshape2)
library(ggplot2)

source(file="R_Code/Spatial_Data_Generation_Pop.R")

simulationWrapper <- function(i) {
  # Sample:
  sampledPop <- Pop[[i]][sample_id[, i], ]
  
  # Aggregieren der Daten:
  aggSample <- sampledPop %.% 
    group_by(clusterid) %.%
    summarise(n = length(y),
              vardir = var(y),
              y =  mean(y), x = mean(x)) %.%
    mutate(vardirsmooth = sum(vardir * (n - 1)) / (sum(n) - length(n)) / n )
  
  # Wahre Werte:
  trueY <- Pop[[i]] %.% 
    group_by(clusterid) %.% 
    summarise(var = var(y), trueY = mean(y), v = mean(v))
  
  # Zusammenspielen der wahren Werte und aggregiertem Sample
  aggSample <- left_join(aggSample, trueY)
  aggSample$var1 <- sd(residuals(summary(lm(y ~ x, aggSample))))^2/10
  aggSample <- aggSample[order(aggSample$clusterid), ]
  
  # Schätzung des FH-Models
  # Robust:
  if(class(try({
  fitRFH <- fitfh(formula = y ~ x, vardir="vardirsmooth", idName="clusterid", data = aggSample)
  # Nicht-Robust:
  fitFH <- eblupFH(y ~ x, vardirsmooth, method = "REML", data = aggSample)})) == "try-error") return(NULL)
  
  # Zusammenfassung der Ergebniss:
  aggSample$rfhY <- fitRFH$prediction
  aggSample$fhY <- fitFH$eblup
  aggSample$r <- i
  aggSample$scenario <- "(0,0,0,0)"
  simulatedData <- aggSample %.% select(scenario, r, clusterid, trueY, rfhY, fhY, y) %.%
    group_by(scenario, r, clusterid)
  
  # Zusammenfassung der Simulation
  tmp <- summary(fitRFH)
  attr(simulatedData, "simMetaData") <- 
    data.frame(scenario = "(0,0,0,0)",
               r = i,
               proc = c("fh", "rfh"),
               iter = c(fitFH$iterations, max(fitRFH$fitparam$m)),
               refvar = c(fitFH$fit$refvar, tmp$coeffRandom))
  
  attr(simulatedData, "simMetaData")[, rownames(tmp$coeffFixed)] <- 
    rbind(fitFH$fit$estcoef$beta, as.numeric(tmp$coeffFixed))
  
  # Rückgabe
  class(simulatedData) <- c("simDataFH", class(simulatedData))
  simulatedData
}

rbind.simDataFH <- function(..., deparse.level = 1) {
  args <- list(...)
  attrList <- lapply(args, function(x) attr(x, "simMetaData"))
  args <- do.call("rbind.data.frame", args)
  attr(args, "simMetaData") <- do.call("rbind", attrList)
  args
}

tmp <- do.call(rbind, lapply(1:length(Pop), simulationWrapper))


dat <- melt(tmp, id.vars = c("scenario", "r", "clusterid", "trueY")) %.% 
  group_by(scenario, clusterid, variable) %.%
  summarise(BIAS = calcBIAS(trueY, value),
            ABIAS = calcABIAS(trueY, value),
            RBIAS = calcRBIAS(trueY, value),
            MSE = calcMSE(trueY, value),
            RRMSE = calcRRMSE(trueY, value))

ggplot(dat) + geom_boxplot(aes_string(x = "variable", y = "RBIAS"))
ggplot(dat) + geom_boxplot(aes_string(x = "variable", y = "RRMSE"))







