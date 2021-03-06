library(saedevel)
library(sae)
library(dplyr)
library(reshape2)
library(ggplot2)
i <- 1
simulationWrapper <- function(i) {
  # Sample:
  sampledPop <- Pop[[i]][sample_id[, i], ]
  
  # Aggregieren der Daten:
  tmp <- var(Pop[[i]]$e)
  
  aggSample <- sampledPop %.% 
    group_by(clusterid) %.%
    summarise(n = length(y),
              vardir = var(y),
              y =  mean(y), x = mean(x)) %.%
    mutate(vardirsmooth = sum(vardir * (n - 1)) / (sum(n) - length(n)) / n ,
           varTrue =  tmp / n)
  
  # Wahre Werte:
  trueY <- Pop[[i]] %.% 
    group_by(clusterid) %.% 
    summarise(var = var(y), trueY = mean(y), v = mean(v))
  
  # Zusammenspielen der wahren Werte und aggregiertem Sample
  aggSample <- left_join(aggSample, trueY, by = "clusterid")
  aggSample$var1 <- sd(residuals(summary(lm(y ~ x, aggSample))))^2/10
  aggSample <- aggSample[order(aggSample$clusterid), ]
  
  tmp <- lm(y~x, data = aggSample)
  var(tmp$residuals)
  
  # Schätzung des FH-Models
  # Robust:
  if(class(try({
    fitRFH <- fitfh(formula = y ~ x, vardir="varTrue", idName="clusterid", data = aggSample, optsRobust = genOptsRobust(k = 1000))
    # Nicht-Robust:
    fitFH <- eblupFH(y ~ x, varTrue, method = "REML", data = aggSample)})) == "try-error") return(NULL)
  
  # Zusammenfassung der Ergebniss:
  aggSample$rfhY <- as.numeric(fitRFH$prediction)
  aggSample$fhY <- as.numeric(fitFH$eblup)
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
