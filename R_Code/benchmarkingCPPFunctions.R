library(SAE)
library(microbenchmark)
rho <- 0.5

nDomains <- 100
nTime <- 20
W <- wMatrix(nDomains)

ome2 <- updateOmega2(rho, nTime)
ome1 <- reSAR1(W=W, rho=rho)
Z1 <- reZ1(nDomains, nTime)


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
tmp <- microbenchmark({A <- updateA(sigma2=1, Ome2=ome2, nDomains=nDomains, nTime = nTime, sigmaSamplingError=1:(nDomains * nTime))
                       tmp1 <- updateSolvedV(sarCorr=rho, sigma1=1, arCorr=rho, sigma2=1, A=A, Ome1=ome1, Z1=Z1)},
                      tmp2 <- matVSolved2(W=W, rho1=rho, sigma1=1, sigma2=1, Ome2=ome2, 
                                          Z1=Z1, nDomains=nDomains, sigmaSamplingError=1:(nDomains * nTime)))



mean(tmp1 - tmp2)

plot(tmp)

