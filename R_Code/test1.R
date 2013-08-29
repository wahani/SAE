# Siehe Funktion SAE::estimateRE
modelSpecs <- tmp1
n <- modelSpecs$nDomains*modelSpecs$nTime
# Sampling Error Component
R.tmp=modelSpecs$sigmaSamplingError*diag(1, modelSpecs$nDomains*modelSpecs$nTime)
svd.R.tmp=svd(R.tmp)
sqrt.R.tmp.inv=solve(t(svd.R.tmp$v%*%(t(svd.R.tmp$u)*sqrt(svd.R.tmp$d))))

# RE - Spatio-Temporal
ome1 <-  updateOmega1(sarCorr=modelSpecs$rho[1], w0=modelSpecs$w0)
ome2Tmp <- updateOmega2(arCorr=modelSpecs$rho[2], nTime=modelSpecs$nTime)
ome2 <- omega2Diag(Ome2=ome2Tmp, nDomains=modelSpecs$nDomains)
G.tmp <- matrix(0, ncol = modelSpecs$nDomains + modelSpecs$nDomains * modelSpecs$nTime,
            nrow = modelSpecs$nDomains + modelSpecs$nDomains * modelSpecs$nTime)
G.tmp[1:modelSpecs$nDomains, 1:modelSpecs$nDomains] <- modelSpecs$sigma[1] * ome1
G.tmp[(modelSpecs$nDomains+1):(modelSpecs$nDomains * modelSpecs$nTime + modelSpecs$nDomains), 
  (modelSpecs$nDomains+1):(modelSpecs$nDomains * modelSpecs$nTime + modelSpecs$nDomains)] <- modelSpecs$sigma[2] * ome2

svd.G.tmp=svd(G.tmp)
sqrt.G.tmp.inv=solve(t(svd.G.tmp$v%*%(t(svd.G.tmp$u)*sqrt(svd.G.tmp$d))))

# Variance-Covariance
z <- modelSpecs$Z
V <- z%*%G.tmp%*%t(z) + R.tmp

A <- updateA(sigma2 = modelSpecs$sigma[2], Ome2=ome2Tmp, nDomains = modelSpecs$nDomains, nTime= modelSpecs$nTime,
             modelSpecs$sigmaSamplingError)
#V <- updateV(sigma1=modelSpecs$sigma[1], Ome1=Ome1, A=A, Z1=modelSpecs$Z1)
Vinv <- updateSolvedV(sarCorr=modelSpecs$rho[1], sigma1=modelSpecs$sigma[1], 
                      arCorr=modelSpecs$rho[2], A=A, Ome1=ome1, Z1=modelSpecs$Z1)

sqrt.u <- updateSqrtU(V=V)
sqrt.u.inv <- diag(1/diag(sqrt.u))

# Starting Values
y <- modelSpecs$y
XS <- modelSpecs$x
beta.q <- modelSpecs$beta
vv.tmp<-G.tmp%*%t(z)%*%Vinv%*%as.vector(y-XS%*%beta.q)

# Algorithm
n <- modelSpecs$nDomains*modelSpecs$nTime
areanumber <- modelSpecs$nDomains
my.psi <- psiOne
tol <- modelSpecs$tol
diff.u<-1
iter2<-0
k_v <- 1.345
maxit <- modelSpecs$maxIter
while (abs(diff.u)>tol)
{
  cat(".")
  iter2<-iter2+1 
  v_robust=as.vector(vv.tmp)
  res1<-sqrt.R.tmp.inv%*%(y-XS%*%beta.q-z%*%v_robust)
  res2<-sqrt.G.tmp.inv%*%v_robust
  w2<-diag(c(my.psi(res1,k_v)/res1),n,n)
  w3<-diag(c(my.psi(res2,k_v)/res2),n + areanumber,n + areanumber)
  Atmp1 <- t(z)%*%(sqrt.R.tmp.inv)%*%w2%*%(sqrt.R.tmp.inv)%*%z
  Atmp2 <- sqrt.G.tmp.inv%*%w3%*%sqrt.G.tmp.inv
  A=Atmp1+Atmp2
  B<-t(z)%*%(sqrt.R.tmp.inv)%*%w2%*%(sqrt.R.tmp.inv)%*%(y-XS%*%beta.q)
  vv.tmp<-solve(A)%*%B
  
  diff.u<-sum(c((vv.tmp-v_robust)^2))
  if (iter2>maxit)
  {warning(paste("failed to converge in", maxit, "steps"))
   break}
}

rand.eff.robust=as.numeric(vv.tmp)
summary(XS%*%modelSpecs$beta + z %*% rand.eff.robust - results1[[1]]$estimates)


