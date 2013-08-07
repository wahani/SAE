rm(list= ls())

# Daten generieren:

# Installiere Paket
# utils::install.packages(pkgs="../spatioTemporalData/spatioTemporalData_1.0.zip", 
#                         repos = NULL)
# utils::install.packages(pkgs="../spatioTemporalData/spatioTemporalData_1.0.tar.gz", 
#                         repos = NULL)


require(spatioTemporalData)
require(SAE)
set.seed(3)
setup <- simSetupMarhuenda(nDomains=30, nTime=5, sarCorr=0.5, arCorr=0.5, n = 200)
output <- simRunMarhuenda(setup)
dat <- slot(output[[1]], "data")[[1]]
sigmaE <- slot(output[[1]], "sigma")

# Prepare Data
# system.time(tmp <- optimizeParameters(modelSpecs))

# Rprof(tmp <- tempfile())
# out <- optimizeParameters(modelSpecs)
# Rprof()
# profileSummary <- summaryRprof(tmp)
# unlink(tmp)
# 
# profileSummary$by.total

fit <- fitSTREBLUP(y~x, dat, c(0,1), c(1,1), c(0.5,0.5))

summary.fitSTREBLUP <- function(fit) {
  out <- matrix(c(fit$beta, fit$sigma, fit$rho), ncol = 1)
  rownames(out) <- c(rownames(fit$beta), "sigmaSquaredSAR", "sigmaSquaredAR", "rhoSAR", "rhoAR")
  colnames(out) <- "estimated coefficient"
  as.table(out)
}
class(fit)
summary(fit)