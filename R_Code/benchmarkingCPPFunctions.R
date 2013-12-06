library(SAE)
library(microbenchmark)
library(stats4)
library(mnormt)
library(nloptr)
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
  
  output <- c(svvpp, recursive=TRUE)
  
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
tmp1$beta
tmp1$rho
tmp1$sigma

# 
# 
# minusloglClosure <- function(y, X, sigmaSamplingError, W, Z1) {
#   function(rho1, sigma1, rho2, sigma2) -llr(y, X, rho1, sigma1, rho2, sigma2, Z1, sigmaSamplingError, W)
# }
# # 
# library(stats4)
# system.time(
#   tmp2 <- mle(minusloglClosure(y = y, X = X, sigmaSamplingError = sigmaSamplingError, W = W, Z1), 
#                 list(rho1 = 0.5, rho2 = 0.5, sigma1 = 1, sigma2 = 1), lower = c(-0.99, 0.01, -0.99, 0.01), method = "L-BFGS-B",
#                 upper = c(0.99, 100000, 0.99, 10000))
# )
# # 
# system.time(
# tmp3 <- mle(
#   minusloglClosure(y = y, X = X, sigmaSamplingError = sigmaSamplingError, W = W, Z1), 
#     list(rho1 = 0.5, rho2 = 0.5, sigma1 = 1, sigma2 = 1), method = "Nelder-Mead")
# )
# nobs(tmp3)
# # 
# # summary(tmp3)
# # 
# system.time(
#   tmp4 <- fitSTEBLUP(y~x, dat, c(100,1), c(10,10), c(0, 0))
# )
# 
# 
# minusloglClosure1 <- function(y, X, sigmaSamplingError, W, Z1) {
#   function(x) llr(y, X, x[1], x[2], x[3], x[4], Z1, sigmaSamplingError, W)
# }

lltClosure <- function(y, X, sigmaSamplingError, W, Z1) {
  function(theta) {
    df <- 2
    lV <- matVinv(W=W, rho1=rev(theta)[4], sigma1 = rev(theta)[3], rho2=rev(theta)[2],
                     sigma2 = rev(theta)[1], Z1, sigmaSamplingError)
    beta <- blue(y=y, X=X, Vinv=lV$Vinv)
    mPred <- as.numeric(X %*% beta)
    pt <- -0.5 * log(1 + df + 4) *  determinant(crossprod(X, lV$Vinv) %*% X)$modulus
    -sum(dmt(as.numeric(y), mean = mPred, 
                           S = lV$V, df = df, log = TRUE)) + pt
  }
}

system.time(
  tmp <- nloptr(c(0.5, 1, 0.5, 1), 
                lltClosure(y = y, X = X, sigmaSamplingError = sigmaSamplingError, W = W, Z1),
                lb = c(-0.99, 0.01, -0.99, 0.01), ub = c(0.99, 10000, 0.99, 10000),
                opts = list(algorithm = "NLOPT_LN_NELDERMEAD",
                            xtol_rel = 10^(-3),
                            maxeval = 500)))



# 
# system.time(
#   tmp <- nloptr(c(0, 10, 0, 10), 
#                 minusloglClosure1(y = y, X = X, sigmaSamplingError = sigmaSamplingError, W = W, Z1),
#                 lb = c(-0.99, 0.01, -0.99, 0.01), ub = c(0.99, 10000, 0.99, 10000),
#                 opts = list(algorithm = "NLOPT_LN_NELDERMEAD",
#                             xtol_rel = 10^(-3),
#                             maxeval = 200)))
# 
# algorithms <- c("NLOPT_GN_DIRECT", "NLOPT_GN_DIRECT_L",
#   "NLOPT_GN_DIRECT_L_RAND", "NLOPT_GN_DIRECT_NOSCAL",
#   "NLOPT_GN_DIRECT_L_NOSCAL",
#   "NLOPT_GN_DIRECT_L_RAND_NOSCAL",
#   "NLOPT_GN_ORIG_DIRECT", "NLOPT_GN_ORIG_DIRECT_L",
#   "NLOPT_GD_STOGO", "NLOPT_GD_STOGO_RAND",
#   "NLOPT_LD_SLSQP", "NLOPT_LD_LBFGS_NOCEDAL",
#   "NLOPT_LD_LBFGS", "NLOPT_LN_PRAXIS", "NLOPT_LD_VAR1",
#   "NLOPT_LD_VAR2", "NLOPT_LD_TNEWTON",
#   "NLOPT_LD_TNEWTON_RESTART",
#   "NLOPT_LD_TNEWTON_PRECOND",
#   "NLOPT_LD_TNEWTON_PRECOND_RESTART",
#   "NLOPT_GN_CRS2_LM", "NLOPT_GN_MLSL", "NLOPT_GD_MLSL",
#   "NLOPT_GN_MLSL_LDS", "NLOPT_GD_MLSL_LDS",
#   "NLOPT_LD_MMA", "NLOPT_LN_COBYLA", "NLOPT_LN_NEWUOA",
#   "NLOPT_LN_NEWUOA_BOUND", "NLOPT_LN_NELDERMEAD",
#   "NLOPT_LN_SBPLX", "NLOPT_LN_AUGLAG", "NLOPT_LD_AUGLAG",
#   "NLOPT_LN_AUGLAG_EQ", "NLOPT_LD_AUGLAG_EQ",
#   "NLOPT_LN_BOBYQA", "NLOPT_GN_ISRES")
# 
# nloptList <- list()
# for(algorithm in algorithms[grep("LD", algorithms)]) {
#   try(nloptList[algorithm] <- list(nloptr(c(0.5, 1, 0.5, 1), 
#                 minusloglClosure1(y = y, X = X, sigmaSamplingError = sigmaSamplingError, W = W, Z1),
#                 lb = c(-0.99, 0.01, -0.99, 0.01), ub = c(0.99, 10000, 0.99, 10000),
#                 opts = list(algorithm = algorithm,
#                             xtol_rel = 10^(-3)))))
# }
# 
# nloptList[sapply(nloptList, function(alg) alg$iterations < 100)]





# Vinv <- matVinv(W, rho1 = 0.5, sigma1 = 1, rho2 = 0.5, sigma2 = 1, Z1, sigmaSamplingError)
# 
# 
ome2 <- updateOmega2(0.5, 10)
ome1 <- reSAR1(W=W, rho=0.5)
# 

# V <- updateV(1, ome1, A, Z1)
# 
# tmp1 <- cppChol(V)
# tmp2 <- chol(V)
# 
# all.equal(chol2inv(tmp1), tcrossprod(cppChol2Inv(tmp1)) )
# 
# 
# all.equal(tmp1, tmp2)
# 
# all.equal(Vinv$V, V)
# 
# all.equal(Vinv$Vinv, solve(V))
# 
# blue(y, X, Vinv$Vinv)

# microbenchmark(
# {A <- updateA(sigma2=1, Ome2=ome2, nDomains=simSetup@nDomains, nTime = simSetup@nTime, sigmaSamplingError=sigmaSamplingError)
#  Ainv1 <- solve(A)},
# Ainv2 <- matAinv(sigma2=1, Ome2=ome2, nDomains=simSetup@nDomains, sigmaSamplingError=sigmaSamplingError))
# 
# all.equal(Ainv1, Ainv2)
# 
# microbenchmark(
# {
# ome2 <- updateOmega2(0.5, 10)
# ome1 <- reSAR1(W=W, rho=0.5)
# A <- updateA(sigma2=1, Ome2=ome2, nDomains=simSetup@nDomains, nTime = simSetup@nTime, sigmaSamplingError=sigmaSamplingError)
# V <- updateV(1, ome1, A, Z1)
# Vinv1 <- solve(V)
# },
# 
# )
# 
# microbenchmark(
# {
#   ome2 <- updateOmega2(0.5, 10)
#   ome1 <- reSAR1(W=W, rho=0.5)
#   A <- updateA(sigma2=1, Ome2=ome2, nDomains=simSetup@nDomains, nTime = simSetup@nTime, sigmaSamplingError=sigmaSamplingError)
#   V <- updateV(1, ome1, A, Z1)
#   Vinv1 <- solve(V)
# },  
# Vinv1 <- matVinv1(W, rho1 = 0.5, sigma1 = 1, rho2 = 0.5, sigma2 = 1, Z1, sigmaSamplingError),
# Vinv2 <- matVinv(W, rho1 = 0.5, sigma1 = 1, rho2 = 0.5, sigma2 = 1, Z1, sigmaSamplingError))
# 
# 
# 
# all.equal(Vinv1, Vinv2$Vinv)



Rprof(tmp <- tempfile())
tmp1 <- fitSTREBLUP(y~x, dat, c(100,1), c(1,1), c(0.5, 0.5))
Rprof()
profileSummary <- summaryRprof(tmp)
unlink(tmp)


profileSummary$by.total









