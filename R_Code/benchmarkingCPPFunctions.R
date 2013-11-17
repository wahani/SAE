library(SAE)
library(microbenchmark)
# rho1 <- 0.5
# rho2 <- 0.5
# sigma1 <- 2
# sigma2 <- 1
# 
# nDomains <- 100
# nTime <- 20
# W <- wMatrix(nDomains)
# 
# ome2 <- updateOmega2(rho, nTime)
# ome1 <- reSAR1(W=W, rho=rho)
# Z1 <- reZ1(nDomains, nTime)
# sigmaSamplingError <- 1:(nDomains * nTime)




# Schnellste Lösung ist matV2
# tmp <- microbenchmark(tmp1 <- matA(sigma2=1, Ome2=ome2, nDomains=nDomains, sigmaSamplingError=1:(nDomains * nTime)),
#                       tmp2 <- updateA(sigma2=1, Ome2=ome2, nDomains=nDomains, nTime = nTime, sigmaSamplingError=1:(nDomains * nTime)))
# 
# tmp <- microbenchmark({A <- matA(sigma2=1, Ome2=ome2, nDomains=nDomains, sigmaSamplingError=1:(nDomains * nTime))
#                 V1 <- matV1(1, ome1, A, Z1)},
#                V2 <- matV2(W=W, rho1=rho, sigma1=1, sigma2=1, Ome2=ome2, Z1=Z1, nDomains=nDomains, sigmaSamplingError=1:(nDomains * nTime)),
# {A <- updateA(sigma2=1, Ome2=ome2, nDomains=nDomains, nTime = nTime, sigmaSamplingError=1:(nDomains * nTime))
#  V3 <- updateV(1, ome1, A, Z1)}
#               )



# R-Implementation ist schneller?
# tmp <- microbenchmark(matVSolved1(rho1=rho, sigma1=1, rho2=rho, sigma2=1, A=A, Ome1=ome1, Z1=Z1),
#                       updateSolvedV(sarCorr=rho, sigma1=1, arCorr=rho, sigma2=1, A=A, Ome1=ome1, Z1=Z1))

# Inverse 'normal bilden' ist schneller als R-Implementierung
# tmp <- microbenchmark(tmp1 <- matV1(rho, 1, 1, ome2, W, 1:(nDomains * nTime), Z1),
#                       tmp2 <- matV(W=W, rho1=rho, sigma1=1, sigma2=1, Ome2=ome2, 
#                                           Z1=Z1, nDomains=nDomains, sigmaSamplingError=1:(nDomains * nTime)))

# library(stats4)
# X <- cbind(1, rnorm(nDomains * nTime, 2, 5))
# y <- 2*X[,2] + rnorm(nDomains * nTime, 5)
# 
# # y, X, rho1, sigma1, rho2, sigma2, Z1, sigmaSamplingError, W
# minusloglClosure <- function(y, X, sigmaSamplingError, W) {
#   function(rho1, sigma1, rho2, sigma2) -llr(y, X, rho1, sigma1, rho2, sigma2, Z1, sigmaSamplingError, W)
# }
# 
# 
# start <- 
# 
# mle(minusloglClosure(y = y, X = X, sigmaSamplingError = sigmaSamplingError, W = W), 
#     list(rho1 = 0.1, rho2 = 0.1, sigma1 = 1, sigma2 = 1))
# 

rm(list= ls())

# Daten generieren:

# Installiere Paket
# library(spatioTemporalData)
# library(devtools)
# install_github("spatioTemporalData", username = "wahani", subdir = "package")
# utils::install.packages(pkgs="../spatioTemporalData/spatioTemporalData_1.1.4.tar.gz")

"+.simSetup" <- function(x, y) {
  datListX <- slot(x, "data")
  datListY <- slot(y, "data")
  
  datList <- mapply(function(dat1, dat2) data.frame(dat1, dat2$yHat),
                    datListX, datListY, SIMPLIFY = FALSE)
  
  slot(x, "data") <- datList
  x
}

simScenario <- function(sigma, sigmaCont, ...) {
  require(spatioTemporalData)
  set.seed(1)
  svv00 <- simRunContamination(sarCorr=0, arCorr=0, 
                               spatialCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 3),
                               temporalCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 2),
                               spatioTemporalMessup = TRUE, scenarioName = "(v1, v2, 0, 0)",
                               ...)
  
  svvpp <- simRunContamination(sarCorr=0.5, arCorr=0.5, 
                               spatialCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 3),
                               temporalCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 2),
                               spatioTemporalMessup = TRUE, scenarioName = "(v1, v2, 0.5, 0.5)",
                               ...)
  
  s0000 <- simRunContamination(sarCorr=0, arCorr=0,
                               spatialCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 0),
                               temporalCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 0),
                               spatioTemporalMessup = FALSE, scenarioName = "(0, 0, 0, 0)",
                               ...)
  
  s00pp <- simRunContamination(sarCorr=0.5, arCorr=0.5,
                               spatialCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 0),
                               temporalCont = list(sigma = sigma, sigmaCont = sigmaCont, nDomainsCont = 0),
                               spatioTemporalMessup = FALSE, scenarioName = "(0, 0, 0.5, 0.5)", 
                               ...)
  
  output <- c(s00pp, recursive=TRUE)
  
}

output <- simScenario(n = 4, nDomains = 100, nTime = 10, beta = c(100, 1), sigma = 1, sigmaCont = 40,
                      xdt = spGenerator, seVar = seSigmaClosure)[-2]



# simResults <- lapply(output, getSimResults, 
#                      fitFunction = c("fitSTREBLUP"),
#                      mc.cores = 1)


simSetup <- output[[1]]
dat <- simSetup@data[[1]]
y <- dat$y
X <- cbind(1, dat$x)
sigmaSamplingError <- simSetup@sigmaSE
W <- wMatrix(nDomains=simSetup@nDomains)
Z1 <- reZ1(nDomains=simSetup@nDomains, nTime=simSetup@nTime)


system.time(
tmp1 <- fitSTREBLUP(y~x, dat, c(100,1), c(1,1), c(0.5, 0.5))
)

minusloglClosure <- function(y, X, sigmaSamplingError, W, Z1) {
  function(rho1, sigma1, rho2, sigma2) -llr(y, X, rho1, sigma1, rho2, sigma2, Z1, sigmaSamplingError, W)
}

library(stats4)
system.time(
  tmp2 <- mle(minusloglClosure(y = y, X = X, sigmaSamplingError = sigmaSamplingError, W = W, Z1), 
                list(rho1 = 0.5, rho2 = 0.5, sigma1 = 1, sigma2 = 1), lower = c(-0.99, 0.01, -0.99, 0.01), method = "SANN",
                upper = c(0.99, 100000, 0.99, 10000))
)

system.time(
tmp3 <- mle(
  minusloglClosure(y = y, X = X, sigmaSamplingError = sigmaSamplingError, W = W, Z1), 
    list(rho1 = 0.5, rho2 = 0.5, sigma1 = 1, sigma2 = 1), method = "SANN")
)

summary(tmp3)

system.time(
  tmp4 <- fitSTEBLUP(y~x, dat, c(100,1), c(1,1), c(0.5, 0.5))
)




Vinv <- matVinv(W, rho1 = 0.5, sigma1 = 1, rho2 = 0.5, sigma2 = 1, Z1, sigmaSamplingError)


ome2 <- updateOmega2(0.5, 10)
ome1 <- reSAR1(W=W, rho=0.5)

A <- updateA(sigma2=1, Ome2=ome2, nDomains=simSetup@nDomains, nTime = simSetup@nTime, sigmaSamplingError=sigmaSamplingError)
V <- updateV(1, ome1, A, Z1)

tmp1 <- cppChol(V)
tmp2 <- chol(V)

all.equal(chol2inv(tmp1), tcrossprod(cppChol2Inv(tmp1)) )


all.equal(tmp1, tmp2)

all.equal(Vinv$V, V)

all.equal(Vinv$Vinv, solve(V))

blue(y, X, Vinv$Vinv)




