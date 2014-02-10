library(saedevel)
library(sae)
library(dplyr)
detach("package:plyr", unload=TRUE)
library(reshape2)
library(ggplot2)

i <- 1

sampledPop <- Pop[[i]][sample_id[, i], ]
aggSample <- sampledPop %.% 
  group_by(clusterid) %.%
  summarise(n = length(y),
            vardir = sd(y)^2/n*(n-1)/nrow(sampledPop),
            y =  mean(y), x = mean(x))

trueY <- Pop[[i]] %.% group_by(clusterid) %.% summarise(var = sd(y)^2/length(y), trueY = mean(y))

aggSample <- left_join(aggSample, trueY)
aggSample$var1 <- sd(residuals(summary(lm(y ~ x, aggSample))))^2/10
aggSample <- aggSample[order(aggSample$clusterid), ]

fitRFH <- fitfh(formula = y ~ x, vardir="var1", idName="clusterid", data = aggSample)
fitFH <- eblupFH(y ~ x, var1, method = "REML", data = aggSample)

ggplot(melt(data.frame(
  DIRECTBIAS = (aggSample$trueY - aggSample$y)/aggSample$trueY,
  RFHBIAS = (aggSample$trueY - predict(fitRFH, type = "EBLUP"))/aggSample$trueY, 
  FHBIAS = (aggSample$trueY - fitFH$eblup)/aggSample$trueY))) +
  geom_boxplot(aes(y = value, x = variable))


summary(fitRFH)
plot(resid(fit, "sampling"), resid(fit, "RE"))

alt <- eblupFH(yi ~ MajorArea, SD, method = "REML", MAXITER = 100, PRECISION = 0.0001, data = milk)

plot(predict(fit))
