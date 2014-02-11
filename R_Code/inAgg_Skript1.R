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

fitRFH <- fitfh(formula = y ~ x, vardir="bootSD", idName="clusterid", data = aggSample)
fitFH <- eblupFH(y ~ x, bootSD, method = "REML", data = aggSample)

ggplot(melt(data.frame(
  DIRECTBIAS = (aggSample$trueY - aggSample$y)/aggSample$trueY,
  RFHBIAS = (aggSample$trueY - predict(fitRFH, type = "EBLUP"))/aggSample$trueY, 
  FHBIAS = (aggSample$trueY - fitFH$eblup)/aggSample$trueY))) +
  geom_boxplot(aes(y = value, x = variable))


summary(fitRFH)
plot(resid(fit, "sampling"), resid(fit, "RE"))

alt <- eblupFH(yi ~ MajorArea, SD, method = "REML", MAXITER = 100, PRECISION = 0.0001, data = milk)

plot(predict(fit))



bootstrapMean <- do.call("rbind",
                         lapply(split(Pop[[1]], Pop[[1]]$clusterid),
                                function(dat) {
                                  data.frame(k = 1:500, clusterid = unique(dat$clusterid),
                                             y = replicate(500, mean(dat[sample(nrow(dat), 5), "y"]), simplify = "numeric"))
                                }))

bootSD <- bootstrapMean %.% group_by(clusterid) %.% summarise(sd = sd(y))
aggSample$bootSD <- bootSD[order(bootSD$clusterid), "sd"]












