require(spatioTemporalData)
require(SAE)
set.seed(3)
setup <- simSetupMarhuenda(nDomains=100, nTime=10, sarCorr=0.5, arCorr=0.5, n = 50)
output <- simRunMarhuenda(setup)[[1]]
datList <- output@data

results <- mclapply(datList,
                    function(dat) {
                      fit <- fitSTREBLUP(dat = dat, 
                                         formula = y ~ x, 
                                         beta = c(0,1), 
                                         sigma = c(1,1), 
                                         rho = c(0.5,0.5))
                      fit
                    }, 
                    mc.cores = if (grepl("Windows", Sys.getenv("OS"))) 
                      1L else detectCores(),
                    mc.preschedule = FALSE)

results1 <- mclapply(datList,
                     function(dat) {
                       fit <- fitSTEBLUP(dat = dat, 
                                         formula = y ~ x, 
                                         beta = c(0,1), 
                                         sigma = c(1,1), 
                                         rho = c(0.5,0.5))
                       fit
                     }, 
                     mc.cores = if (grepl("Windows", Sys.getenv("OS"))) 
                       1L else detectCores(),
                     mc.preschedule = FALSE)



tmp <- lapply(results1, function(tmp)
  data.frame("value" = c(tmp$beta, tmp$sigma, tmp$rho),
             "parameter" = c("beta1", "beta2", "sigma1", "sigma2", "rho1", "rho2"),
             "model" = "ST"))
tmp <- do.call(rbind, tmp)
evalData$model <- "STR"
evalData <- rbind(evalData, tmp)
require(ggplot2)
resultPlot <- ggplot(evalData, aes(x = value, fill = model, colour = model)) + 
  geom_density(alpha = 0.1) + 
  facet_grid(parameter~., scales="free_y") + 
  labs(title = "Scenario: 100 Domains, 10 Time, TrueParameters: betas: 0, 1, \n rhos: 0.5, 0.5 sigmas: 1, 1, n = 50")

ggsave(filename = "ParameterRobustNonRobustDensities.pdf", width = 10, height = 5)


results1
results[[1]]$estimates - results1[[1]]$estimates

tmp1 <- prepareData(formula=y~x, dat=datList[[1]], nDomains=100, nTime=10, 
                    beta=as.numeric(results1[[1]]$beta), 
                    sigma=as.numeric(results1[[1]]$sigma), 
                    rho = as.numeric(results1[[1]]$rho), 
                    sigmaSamplingError = unlist(lapply(1:100, seSigmaClosure(100, 10), t = 1:10)), 
                    w0=w0Matrix(nDomains=100), w=wMatrix(100), 
                    tol = 0.0003, method = "", maxIter = 500)

tmp1 <- estimateRE(tmp1)