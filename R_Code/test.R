rm(list= ls())

# Daten generieren:

# Installiere Paket
# utils::install.packages(pkgs="../spatioTemporalData/spatioTemporalData_1.0.zip", 
#                         repos = NULL)
# utils::install.packages(pkgs="../spatioTemporalData/spatioTemporalData_1.0.tar.gz", 
#                         repos = NULL)


require(spatioTemporalData)
require(SAE)
set.seed(1)
setup <- simSetupMarhuenda(nDomains=50, nTime=10, sarCorr=0.5, arCorr=0.5)
output <- simRunMarhuenda(setup)
dat <- slot(output[[1]], "data")[[1]]
sigmaE <- slot(output[[1]], "sigma")

# Prepare Data
prepareData <- function(dat, sigmaE, beta, sigma, rho) {
  
  dat <- dat[order(dat$Domain, dat$Time), ]
  XY <- makeXY(y ~ x, dat)
  
  modelSpecs <- list(y = XY$y,
                     x = XY$x,
                     nDomains = getNDomains(dat),
                     nTime = getNTime(dat),
                     
                     w0 = w0Matrix(getNDomains(dat)),
                     w = wMatrix(getNDomains(dat)),
                     
                     Z = reZ(getNDomains(dat), getNTime(dat)),
                     Z1 = reZ1(getNDomains(dat), getNTime(dat)),
                     beta = beta,
                     sigma = sigma,
                     rho = rho,
                     sigmaE = sigmaE,
                     tol = 1e-3,
                     psiFunction = psiOne,
                     K = 0.71,
                     method = "Nelder-Mead")
  return(modelSpecs)
}

modelSpecs <- prepareData(dat, sigmaE, c(0,1), c(2,1), c(0.5,0.5))
out <- optimizeParameters(modelSpecs)




# system.time(tmp <- optimizeParameters(modelSpecs))

# Rprof(tmp <- tempfile())
# out <- optimizeParameters(modelSpecs)
# Rprof()
# profileSummary <- summaryRprof(tmp)
# unlink(tmp)
# 
# profileSummary$by.total



