rm(list= ls())

# Daten generieren:

# Installiere Paket
# utils::install.packages(pkgs="../spatioTemporalData/spatioTemporalData_1.0.zip", 
#                         repos = NULL)

require(spatioTemporalData)
require(SAE)
set.seed(2)
setup <- simSetupMarhuenda(nDomains=50, nTime=10, sarCorr=0.5, arCorr=0.5)
output <- simRunMarhuenda(setup)
dat <- slot(output[[1]], "data")[[1]]
sigmaE <- slot(output[[1]], "sigma")

# Prepare Data
dat <- dat[order(dat$Domain, dat$Time), ]
XY <- makeXY(y ~ x, dat)

getNDomains <- function(dat) length(unique(dat$Domain))
getNTime <- function(dat) length(unique(dat$Time))

modelSpecs <- list(y = XY$y,
                   x = XY$x,
                   nDomains = getNDomains(dat),
                   nTime = getNTime(dat),
                   
                   w0 = w0Matrix(getNDomains(dat)),
                   w = wMatrix(getNDomains(dat)),
                   
                   Z = reZ(getNDomains(dat), getNTime(dat)),
                   Z1 = reZ1(getNDomains(dat), getNTime(dat)),
                   beta = c(1, 1),
                   sigma = c(1, 1),
                   rho = c(0.5, 0.5),
                   sigmaE = sigmaE,
                   tol = 1e-3,
                   psiFunction = psiOne,
                   K = 0.71,
                   method = "SANN")

tmp <- optimizeParameters(modelSpecs)


