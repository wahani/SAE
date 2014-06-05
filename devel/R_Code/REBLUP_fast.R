
#Robust EBLUP
library(MASS)
library(nlme)
library(lme4)
library(rsae)

REBLUP_fast<-function(y, x, Xmean, saind, x.out=NULL,tol, start_true = FALSE,maxit,k,k_ccst,k_v,MSE=FALSE)
{
  
  ################################################################################
  ############# Hubert psi function ##############################################
  ################################################################################
  
  hub.psi <- function(x, k)
  {
    psi <- ifelse(abs(x) <= k, x, sign(x) * k)
    der.psi <- ifelse(abs(x) <= k, 1, 0)
    list(psi = psi, der.psi = der.psi)
  }

  my.psi<-function(u,c){
    
    sm<-median(abs(u))/0.6745
    w <- psi.huber(u/sm,c)
    w*u
  }
  n=length(y)
  nii<-as.numeric(table(saind))
  ni=table(saind)
  samplesizez<-length(saind)
  areanumber<- length(unique(saind))
  Z=matrix(0,n,areanumber)
  kk=0
  for (j in 1:areanumber){for(i in 1:ni[j]){
    kk=kk+1
    Z[kk,j]=1}}
  
  
  if(start_true==TRUE){
    sigmasq0.v<-10
    sigmasq0<-40
    mod   <- lmer(y ~ x + (1|saind))
  }else{
    # Estimation of starting values  
    w=cbind(cbind(1,x),Z)
    fit.H=lm(y~w)
    e.1.H=residuals(fit.H)
    sigma.e.hat.H=sum(e.1.H^2)/(samplesizez-(qr(w)$rank))
    
    xx=cbind(rep(1,samplesizez),x)
    fit2.H=lm(y~xx)
    e.2.H=residuals(fit2.H)
    A.H= sum(e.2.H^2)
    xx=as.matrix(xx)
    B.H=sigma.e.hat.H*(samplesizez-qr(xx)$rank)
    o.H=diag(samplesizez)-(xx%*%(ginv(t(xx)%*%xx)%*%t(xx)))
    C.H=sum(diag(t(Z)%*%o.H%*%Z))
    sigma.v.hat.H=(A.H-B.H)/C.H
    sigmasq0= sigma.e.hat.H #initial values of sigma_sq_e#
    sigmasq0.v=sigma.v.hat.H #initial values sigma_sq_V#
    
    
    
  }
  ls <- lm(y ~ x)
  beta1<-c(ls$coefficients[1], ls$coefficients[2])
  Xpop<-cbind(1,Xmean)
  XS             <- cbind(1, x)
 
  
  const          <- 2 * pnorm(k) - 1 - 2 * k * dnorm(k) + 2 * (k^2) * (1-pnorm(k))
  K              <- diag(rep(const, samplesizez))

  sigma1<-c(sigmasq0,sigmasq0.v)
  rm("sigmasq0","sigmasq0.v")
  ZZ<-Z%*%t(Z)
  
  iter=1
  # Iterations  
  while(TRUE) 
  {
     xx=cbind(rep(1,samplesizez),x)
     xbeta <- c(xx %*% beta1)
     R <- Diagonal(x = rep(sigma1[1], samplesizez))
     G <- Diagonal(x = rep(sigma1[2], areanumber))
     V <- R + Z %*% G %*% t(Z)
     V.inv <- chol2inv(chol(V))
     U <- Diagonal(x=diag(V))
     U.inv <- chol2inv(chol(U))
     r <- c(sqrt(U.inv) %*% (y - xbeta))
     psi.fn <- hub.psi(as.vector(r[[1]]), k = k)
     psi.r <- psi.fn$psi
     der.psi.r <- psi.fn$der.psi
     qq <- V.inv %*% sqrt(U) %*% psi.r
     q1 <- t(xx) %*% qq
     KK5<-sqrt(U) %*% V.inv
     M1 <- t(xx) %*% sqrt(U.inv)%*%Diagonal(x=der.psi.r)%*%KK5%*%xx
     beta <- beta1 + c(solve(M1) %*% q1)[[1]]
     KK0<-t(psi.r)%*% KK5
     A1=KK0%*%qq
     A2=KK0%*%ZZ%*%qq
     A <- matrix(c(A1[1],A2[1]), nrow = 2, ncol=1)
     KK1<-K%*%V.inv
     KK2<-KK1%*%diag(samplesizez)%*%V.inv
     KK3<-KK1%*%ZZ%*%V.inv
     t1=sum(diag(KK2%*%diag(samplesizez)))
     t2=sum(diag(KK2%*%ZZ))
     t3=sum(diag(KK3%*%diag(samplesizez)))
     t4=sum(diag(KK3%*%ZZ))
     T1<- matrix(c(t1,t3,t2,t4), nrow = 2, ncol=2)
     sigma=chol2inv(chol(T1))%*%A
    
    fehler         <- sum(abs(beta[1] - beta1[1]) + abs(beta[2] -beta1[2]) + abs(sigma[2] -sigma1[2]) +  abs(sigma[1] -sigma1[1]))
    ifelse (iter < maxit, ifelse(fehler > tol, {iter=iter+1}, {break}), {break})
    beta1 <- as.vector(beta) 
    sigma1<-sigma

  }
  
  conv           <- ifelse(iter < maxit, 1, 0)
  n.iter         <- iter
  
  # Prediction Random Effects
  
  it<-1
  I.k            <- diag(rep(1,areanumber))
  
  R <- diag(rep(sigma[1], samplesizez))
  G <- diag(rep(sigma[2], areanumber))
  V              <- R + Z %*% G %*% t(Z)
  
  G.inv          <- chol2inv(chol(G))
  V1             <- diag(diag(V))
  R.inv          <- chol2inv(chol(R))
  V.inv          <- chol2inv(chol(V))
  
  V1             <- diag(diag(V))
  V1.inv         <- chol2inv(chol(V1))
  
  xbeta          <- c(XS %*% beta1)
  eblup.v        <- c(G %*% t(Z) %*% V.inv %*% (y - xbeta))
  
  r              <- c(sqrt(V1.inv) %*% (y - xbeta))
  psi.r          <- hub.psi(r, k_v)$psi
  
  sum.v          <- NULL
  vv             <- NULL
  v0             <- eblup.v
  
  repeat
  {
    J            <- c(rep(0,it-1), -1, 1)
    J            <- J[2:(it+1)]
    
    zv0            <- rep(v0,times=c(nii))
    r1             <- c(sqrt(R.inv) %*% (y - xbeta - zv0))
    psi.rr         <- hub.psi(c(r1), k_v)
    psi.r1         <- psi.rr$psi
    der.psi.r1     <- psi.rr$der.psi
    
    v1             <-  c(sqrt(G.inv) %*% v0)
    psi.v          <- hub.psi(v1, k_v)
    psi.v1         <- psi.v$psi
    der.psi.v1     <- psi.v$der.psi
    
    E.psi.der      <- pnorm(k_v) - pnorm(-k_v)
    
    H              <- t(Z) %*% sqrt(R.inv) %*%  psi.r1 - sqrt(G.inv) %*% psi.v1
    
    #	HH             <- t(Z) %*% sqrt(solve(R)) %*% diag(der.psi.r1) %*% sqrt(solve(R)) %*% Z + 
    #            				sqrt(solve(G)) %*% diag(der.psi.v1) %*% sqrt(solve(G))
    
    HH             <- t(Z) %*% sqrt(R.inv) %*% diag(rep(E.psi.der, samplesizez)) %*% sqrt(R.inv) %*% Z + 
      sqrt(G.inv) %*% diag(rep(E.psi.der, areanumber)) %*% sqrt(G.inv)
    
    v0             <- v0 + c(solve(HH) %*% H)			
    vv             <- cbind(vv, v0)
    
    sum.v0         <- sum(abs(v0))
    sum.v          <- c(sum.v, sum.v0)
    
    erreur         <- abs(sum(J * sum.v))
    ifelse (it < maxit, ifelse(erreur > tol, {it=it+1}, {break}), {break})  
  }
  
  conv           <- ifelse(it < maxit, 1, 0)
  n.iter         <- it
  
  rob.eblup.v       <- v0
  
  
  #REBLUP estimator: G is the variance-cov matrix for the random effects
  y.hat.sample<-NULL
  nii<-as.numeric(table(saind))
  if(length(x.out)==0){
    
    thetaREBLUP<-(Xpop%*%beta1) + rob.eblup.v
    
    
  }else{
    sample.sizer<-as.numeric(table(x.out$clusterid))
    Ni<-sample.sizer+nii
    Xmean_out<-cbind(1,tapply(x.out$x,x.out$clusterid,mean))
    y.sample<-tapply(y,saind,sum)
    thetaREBLUP<-as.vector((Xpop%*%beta1) + rob.eblup.v)
    thetaREBLUP_out<-as.vector((Xmean_out%*%beta1) + rob.eblup.v)
    y.hat.sample<-as.vector((1/Ni)*((y.sample)+sample.sizer*(thetaREBLUP_out)))
  }
  
 
  
  ###################################################
  # MSE
  ###################################################

  REBLUP.Mean_SR_m<-NULL
  REBLUP.MSE.Robust_SR_m<-NULL
  REBLUP.MSE.Robust_SR_m_conditional<-NULL
  if(MSE==TRUE){
  z<-Z
  tmp<-V
  tmp.inv<-V.inv
  tmp.U<-diag(diag(tmp),n,n)
  svd.tmp.U<-svd(tmp.U)
  sqrt.ui<-solve(t(svd.tmp.U$v%*%(t(svd.tmp.U$u)*sqrt(svd.tmp.U$d))))
  sqrt.u<-t(svd.tmp.U$v%*%(t(svd.tmp.U$u)*sqrt(svd.tmp.U$d)))
  res1<-(sqrt.ui%*%y-sqrt.ui%*%XS%*%beta1)
  w1<-diag(c(my.psi(res1,k)/res1),n,n) # w1 aus Chambers et al. 2013 Seite 7
  R.tmp<-R
  svd.R.tmp<-svd(R.tmp)
  sqrt.R.tmp.inv=solve(t(svd.R.tmp$v%*%(t(svd.R.tmp$u)*sqrt(svd.R.tmp$d))))
  
  G.tmp<-G
  svd.G.tmp=svd(G.tmp)
  sqrt.G.tmp.inv=solve(t(svd.G.tmp$v%*%(t(svd.G.tmp$u)*sqrt(svd.G.tmp$d))))
  
  v_robust<-as.vector(rob.eblup.v)
  res1<-sqrt.R.tmp.inv%*%(y-XS%*%beta1-z%*%v_robust)
  res2<-sqrt.G.tmp.inv%*%v_robust
  w2<-diag(c(my.psi(res1,k_v)/res1),n,n)
  w3<-diag(c(my.psi(res2,k_v)/res2),areanumber,areanumber)
  
  huber <- function(x){x_rob  <-rep(0,times=length(x))
                       for(i in 1:length(x)){x_rob[i] <-   x[i]*min(1,k/abs(x[i]))}
                       return(x_rob)}
  K2<-mean(huber(rnorm(10000,0,1))^2)
  AA<-solve(t(XS)%*%tmp.inv%*%sqrt.u%*%w1%*%sqrt.ui%*%XS)%*%t(XS)%*%tmp.inv%*%sqrt.u%*%w1%*%sqrt.ui # Matrix A_s aus Chambers et al. 2013 Seite 7
  BB<-solve(t(z)%*%(sqrt.R.tmp.inv)%*%w2%*%(sqrt.R.tmp.inv)%*%z+sqrt.G.tmp.inv%*%w3%*%sqrt.G.tmp.inv)%*%t(z)%*%(sqrt.R.tmp.inv)%*%w2%*%(sqrt.R.tmp.inv) # Matrix B_s aus Chambers et al. 2013 Seite 7
  Vi<-V.inv 
  U<-diag(c(diag(V)),n,n)
  svd.U.tmp=svd(U)  # Matrix U 
  sqrt.U.tmp=t(svd.U.tmp$v%*%(t(svd.U.tmp$u)*sqrt(svd.U.tmp$d)))  # U^{1/2}
  sqrt.U.tmp.inv=solve(t(svd.U.tmp$v%*%(t(svd.U.tmp$u)*sqrt(svd.U.tmp$d))))  # U^{-1/2}
  
  u.hat<-c(z%*%solve(t(z)%*%sqrt.R.tmp.inv%*%w2%*%sqrt.R.tmp.inv%*%z)
           %*%(t(z)%*%sqrt.R.tmp.inv%*%w2%*%sqrt.R.tmp.inv%*%z+sqrt.G.tmp.inv
               %*%w3%*%sqrt.G.tmp.inv)%*%rob.eblup.v) # Was ist u.hat
  
  
  r.cond<- NULL
  for(i in 1:areanumber){
    yi<-y[saind==i]
    xi<-x[saind==i]
    ui<-u.hat[saind==i]
    r.cond<- c(r.cond,(yi-beta1[1]-xi*beta1[2]-ui))
  }
  
  ri<-sqrt.U.tmp.inv%*%(y-beta1[1]-x*beta1[2]) # U^{-1/2}(y-X\beta1)
  ti<-sqrt.R.tmp.inv%*%(y-beta1[1]-x*beta1[2]-z%*%rob.eblup.v) # R^{-1/2}(y-X\beta1-Zv)
  qi<-sqrt.G.tmp.inv%*%rob.eblup.v # G^{-1/2}v^rob
  
  Iri<-diag((-k< c(ri) & c(ri) <k),n,n)
  Iti<-diag((-k< c(ti) & c(ti) <k),n,n)
  Iqi<-diag((-k< c(qi) & c(qi) <k),areanumber,areanumber)
  
  E11<-(-t(cbind(1,x))%*%Vi%*%sqrt.U.tmp%*%Iri%*%sqrt.U.tmp.inv%*%cbind(1,x))
  E12<-matrix(0,2,areanumber)
  
  E21<-(-t(z)%*%sqrt.R.tmp.inv%*%Iti%*%sqrt.R.tmp.inv%*%cbind(1,x))
  E22<-(-t(z)%*%sqrt.R.tmp.inv%*%Iti%*%sqrt.R.tmp.inv%*%z-sqrt.G.tmp.inv%*%Iqi%*%sqrt.G.tmp.inv)
  
  EE<-rbind(cbind(E11,E12),cbind(E21,E22))
  
  p=2
  
  lambda1<-(1+(p/n)*(var(diag(Iri))/mean(diag(Iri))))
  lambda2<-(1+(p/n)*(var(diag(Iti))/mean(diag(Iti))))
  
  tmp.est<-lambda1*sum(my.psi(ri,k)^2)/(n-p)
  tmp1.est<-lambda2*sum(my.psi(ti,k)^2)/(n-p)
  tmp2.est<-sum(my.psi(ri,1.345)*my.psi(ti,k))/(n-p)
  
  V11<-tmp.est*t(cbind(1,x))%*%Vi%*%sqrt.U.tmp%*%sqrt.U.tmp%*%Vi%*%cbind(1,x)
  V22<-tmp1.est*t(z)%*%sqrt.R.tmp.inv%*%sqrt.R.tmp.inv%*%z
  V12<-tmp2.est*t(cbind(1,x))%*%Vi%*%sqrt.U.tmp%*%sqrt.R.tmp.inv%*%z
  V21<-t(V12)
  
  VV<-rbind(cbind(V11,V12),cbind(V21,V22))
  
  Var.Cov.BV<-solve(EE)%*%VV%*%t(solve(EE))
  Var.Cov.B<-Var.Cov.BV[1:2,1:2]
  Var.Cov.V<-Var.Cov.BV[3:(areanumber+2),3:(areanumber+2)]
  
  D<-t(z)%*%sqrt.R.tmp.inv%*%w2%*%sqrt.R.tmp.inv%*%z+sqrt.G.tmp.inv%*%w3%*%sqrt.G.tmp.inv
  
  C1<-solve(D)%*%(0.5*solve(G.tmp)%*%sqrt.G.tmp.inv%*%w3%*%sqrt.G.tmp.inv+0.5*solve(G.tmp)%*%sqrt.G.tmp.inv%*%Iqi%*%sqrt.G.tmp.inv)%*%solve(D)%*%t(z)%*%sqrt.R.tmp.inv%*%w2%*%sqrt.R.tmp.inv
  
  C2<-solve(D)%*%(0.5*t(z)%*%solve(R.tmp)%*%sqrt.R.tmp.inv%*%w2%*%sqrt.R.tmp.inv%*%z+0.5*t(z)%*%solve(R.tmp)%*%sqrt.R.tmp.inv%*%Iti%*%sqrt.R.tmp.inv%*%z)%*%solve(D)%*%t(z)%*%sqrt.R.tmp.inv%*%w2%*%sqrt.R.tmp.inv-solve(D)%*%(0.5*t(z)%*%solve(R.tmp)%*%sqrt.R.tmp.inv%*%w2%*%sqrt.R.tmp.inv+0.5*t(z)%*%solve(R.tmp)%*%sqrt.R.tmp.inv%*%Iti%*%sqrt.R.tmp.inv)
  
  Cov.tmp<-matrix(0,n,n)
  
  v.i<-NULL
  
  for (i in 1:areanumber)
  {
    v.i<-c(v.i,rep(rob.eblup.v[i],ni[i]))
  }
  
  for (i in 1:n)
  {
    for (j in 1:n)
    {if (i==j)
    {
      Cov.tmp[i,j]<-v.i[i]*v.i[j]+sigma1[1]
    }
     
     if (i!=j)
     {
       Cov.tmp[i,j]<-v.i[i]*v.i[j]
     }
     
    }
  }
  
  
  #Cov delta
  
  score.i=function(sigma2u,sigma2e,size,p,yi,xi)
  {
    
    zi=matrix(rep(1,size),size,1)
    Hessian=matrix(0,2,2)
    s<-matrix(0,2,1) 
    Vh<-sigma2e*diag(1,size,size)+sigma2u*zi%*%t(zi)
    Uh<-diag(c(diag(Vh)),size,size)
    Vih<-solve(Vh)
    svd.tmp2=svd(Uh)
    sqrt.v2=solve(t(svd.tmp2$v%*%(t(svd.tmp2$u)*sqrt(svd.tmp2$d))))
    sqrt.v3=(t(svd.tmp2$v%*%(t(svd.tmp2$u)*sqrt(svd.tmp2$d))))
    y.star2=sqrt.v2%*%yi
    x.star2=sqrt.v2%*%xi
    res.new<-sqrt.v2%*%(yi-xi%*%beta1)
    res.robust=my.psi(res.new,k)
    
    
    s[1,1]<-0.5*t(matrix(c(res.robust),size,1))%*%sqrt.v3%*%Vih%*%zi%*%t(zi)%*%Vih%*%sqrt.v3%*%matrix(c(res.robust),size,1)-0.5*K2*sum(diag(Vih%*%zi%*%t(zi)))
    s[2,1]<-0.5*t(matrix(c(res.robust),size,1))%*%sqrt.v3%*%Vih%*%diag(1,size,size)%*%Vih%*%sqrt.v3%*%matrix(c(res.robust),size,1)-0.5*K2*sum(diag(Vih%*%diag(1,size,size)))
    
    Hessian[1,1]=((-1/2)*t(sqrt.v2%*%sqrt.v2%*%(res.new*(-k<=res.new & res.new<k)))%*%sqrt.v3%*%Vih+(1/2)*t(matrix(c(res.robust),size,1))%*%sqrt.v2%*%Vih-t(matrix(c(res.robust),size,1))%*%sqrt.v3%*%Vih%*%zi%*%t(zi)%*%Vih)%*%zi%*%t(zi)%*%Vih%*%sqrt.v3%*%matrix(c(res.robust),size,1)+0.5*K2*sum(diag(Vih%*%zi%*%t(zi)%*%Vih%*%zi%*%t(zi)))
    
    Hessian[1,2]=((-1/2)*t(sqrt.v2%*%sqrt.v2%*%(res.new*(-k<=res.new & res.new<k)))%*%sqrt.v3%*%Vih+(1/2)*t(matrix(c(res.robust),size,1))%*%sqrt.v2%*%Vih-t(matrix(c(res.robust),size,1))%*%sqrt.v3%*%Vih%*%Vih)%*%zi%*%t(zi)%*%Vih%*%sqrt.v3%*%matrix(c(res.robust),size,1)+0.5*K2*sum(diag(Vih%*%zi%*%t(zi)%*%Vih))
    
    
    
    Hessian[2,1]=Hessian[1,2]
    
    Hessian[2,2]=((-1/2)*t(sqrt.v2%*%sqrt.v2%*%(res.new*(-k<=res.new & res.new<k)))%*%sqrt.v3%*%Vih+(1/2)*t(matrix(c(res.robust),size,1))%*%sqrt.v2%*%Vih-t(matrix(c(res.robust),size,1))%*%sqrt.v3%*%Vih%*%Vih)%*%Vih%*%sqrt.v3%*%matrix(c(res.robust),size,1)+0.5*K2*sum(diag(Vih%*%Vih))
    
    
    list(score=s,He=Hessian)
    
  }
  
  si=array(c(rep(0,(areanumber*2))),dim=c(2,1,areanumber))
  Hi=array(c(rep(0,(areanumber*2))),dim=c(2,2,areanumber))
  
  for (i in 1:areanumber)
  {
    tmp=score.i(sigma2u=sigma1[2],sigma2e=sigma1[1],size=ni[i],p=2,yi=y[saind==i],xi=matrix(cbind(1,x[saind==i]),ni[i],2))
    si[,1,i]=tmp$score
    Hi[,,i]=tmp$He
  }
  
  C=matrix(0,2,2)
  
  C[1,1]=(1/areanumber)*sum(Hi[1,1,])
  C[1,2]=(1/areanumber)*sum(Hi[1,2,])
  C[2,1]=(1/areanumber)*sum(Hi[2,1,])
  C[2,2]=(1/areanumber)*sum(Hi[2,2,])
  
  
  B=matrix(0,2,2)
  
  B[1,1]=(1/areanumber)*sum(si[1,1,]*si[1,1,])
  B[1,2]=(1/areanumber)*sum(si[1,1,]*si[2,1,])
  B[2,1]=(1/areanumber)*sum(si[2,1,]*si[1,1,])
  B[2,2]=(1/areanumber)*sum(si[2,1,]*si[2,1,])
  
  
  I.matrix=solve(C)%*%(B)%*%solve(C)
  
  V.sigma=C1%*%Cov.tmp%*%t(C1)*(I.matrix[1,1]/areanumber)+C2%*%Cov.tmp%*%t(C2)*(I.matrix[2,2]/areanumber)+(C1%*%Cov.tmp%*%t(C2)+C2%*%Cov.tmp%*%t(C1))*(I.matrix[1,2]/areanumber)
  
  

  M1<-NULL
  M2<-NULL
  M3<-NULL
  M3Cond<-NULL
  M4<-NULL
  M5<-NULL
  
  for (i in 1:areanumber)
  {
    yi<-y[saind==i]
    xi<-x[saind==i]
    Xr<-x.out[x.out$clusterid==i,3]
    Xr<-cbind(1,Xr)
    
    tr<-apply(Xr,2,mean)
    tr1<-apply(Xr,2,sum)
    ir<-rep(0,n)
    ir[saind==i]<-1
    b<-rep(0,areanumber)
    b[i]<-1
    Is<-diag(1,n,n)
    wi<-(Ni[i]^(-1))*(ir+(Ni[i]-ni[i])*(tr%*%AA+b%*%BB%*%(Is-cbind(1,x)%*%AA))) # Gewichte für REBLUP
    
    REBLUP.Mean_SR_m[i]=sum(wi*y)
    
    wi.REBLUP=wi[saind==i] # Gewichte Sample
    wi.REBLUP.c=wi[saind!=i] # Gewichte Non-Sample
    
    ri.cond<-r.cond[saind==i]
    ci1=(Ni[i]*wi.REBLUP-1)^2
    ci2=(Ni[i]-ni[i])/n
    ci=(ci1+ci2)/(Ni[i]^2)
    ci1.c=(Ni[i]*wi.REBLUP.c)^2
    ci.c=(ci1.c+ci2)/(Ni[i]^2)
    
    yi.c<-y[saind!=i]
    xi.c<-cbind(1,x[saind!=i])
    ri.c.cond<-r.cond[saind!=i]    # Conditional residual
    areaeff<- NULL
    areaeff<-u.hat
    
    areaeffi<-areaeff[saind==i]   	
    areaeffi.c<-areaeff[saind!=i]  
    
    #bias expression		
    Est.REBLUP.MSEBias.Robust<-(t(wi.REBLUP)%*%areaeffi +t(wi.REBLUP.c)%*%areaeffi.c)-unique(areaeffi) # verzerrung des REBLUP vgl. Chambers et al. 2013 Seite 8 oben
    
    cost<-(1-ni[i]/Ni[i])^2
    tmp.matrix=matrix(c(tr,b),ncol=areanumber+2)
    
    # estimation c3
    
    z.i<-matrix(0,1,areanumber)
    z.i[1,i]<-(Ni[i]-ni[i])/Ni[i] # Sampling Anteil 1-ni/N
    
    c3<-as.numeric(z.i%*%V.sigma%*%t(z.i))
    
    M1[i]<-as.numeric(cost*matrix(c(tr),1,2)%*%Var.Cov.B%*%matrix(c(tr),2,1))
    M2[i]<-as.numeric(cost*matrix(c(b),1,areanumber)%*%Var.Cov.V%*%matrix(c(b),areanumber,1))
    M3[i]<-cost*sum(r.cond^2)/((n-1)*(Ni[i]-ni[i]))
    M3Cond[i]<-cost*sum(r.cond[saind==unique(saind)[i]]^2)/((ni[i]-1)*(Ni[i]-ni[i]))
    M4[i]<-c3
    M5[i]<-(Est.REBLUP.MSEBias.Robust)^2 # Verzerrung im Quadrat
    
    
    REBLUP.MSE.Robust_SR_m[i]<-M1[i]+M2[i]+M3[i]+M4[i]+M5[i]
    
    
    REBLUP.MSE.Robust_SR_m_conditional[i]<-M1[i]+M2[i]+M3Cond[i]+M4[i]+M5[i]
    
    
  } 
  }else{}
  
  
  list(est_syn = thetaREBLUP, est_sample = y.hat.sample, 
       mse_conditional = REBLUP.MSE.Robust_SR_m_conditional, mse_linearized = REBLUP.MSE.Robust_SR_m, 
       randeff = rob.eblup.v, beta = beta1, sigma_v = sigma1[2], sigma_e = sigma1[1],
       areasize = Ni, samplesize = nii, iterations = n.iter)
}




