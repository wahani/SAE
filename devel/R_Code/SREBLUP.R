#Robust EBLUP
library(MASS)
library(nlme)
library(lme4)
library(rsae)
library(spdep)

#Functions
SREBLUP_SR<-function(y, x, Xmean, saind, x.out=NULL,tol, start_true = FALSE, maxit, k, k_ccst, k_v, rob_starting=FALSE,W_matrix)
{
  my.psi<-function(u,c){
    
    sm<-median(abs(u))/0.6745
    w <- psi.huber(u/sm,c)
    w*u
  }
    assign("y",y,pos=1)
    assign("x",x,pos=1)
    assign("k",k,pos=1)
    assign("saind",saind,pos=1)
    assign("tol",tol,pos=1)
    assign("k_v",k_v,pos=1)
    assign("maxit",maxit,pos=1)
  
  if(start_true==TRUE){
    startsigma_v<-1
    startsigma_e<-4
    startwert_p<-DataSetup$p
    mod   <- lmer(y ~ x + (1|saind))
  }else{
    mod   <- lmer(y ~ x + (1|saind))
    ymean<-tapply(y,saind, mean)
    xmean<-tapply(x,saind, mean)
    Y.mod <- errorsarlm(as.vector(ymean) ~ as.vector(xmean), listw = W_list)
    
    startwert_p<-as.numeric(Y.mod$lambda) 
    if(rob_starting==TRUE){
      df <- data.frame(y = y, x = x, areaid = saind)
      unitmodel <- saemodel(y ~ x,area = ~areaid, data = df)
      huberfit<-try(fitsaemodel("huberm",unitmodel, k=k, init="s"))
      startsigma_v<-huberfit$theta[2]
      startsigma_e<-huberfit$theta[1]
    }else{
      startsigma_v <- (as.numeric(attr(VarCorr(mod)$saind, "stddev")))^2  # starting value for simga2_v
      startsigma_e <- (as.numeric(attr(VarCorr(mod), "sc")))^2            # starting value for simga2_e
    }
  }
  
  startbeta<-c(as.numeric(fixef(mod)[1]),as.numeric(fixef(mod)[2]))
  
    XS<-cbind(1,x)
    n=length(y)
    ni=table(saind)
    samplesizez<-length(saind)
    areanumber<- length(unique(saind))
    pp=ncol(XS)
    Xpop<-cbind(1,Xmean)
    
    z=matrix(0,n,areanumber)
    kk=0
    for (j in 1:areanumber){for(i in 1:ni[j]){
                kk=kk+1
                z[kk,j]=1}}
  
    huber <- function(x){x_rob  <-rep(0,times=length(x))
                       for(i in 1:length(x)){x_rob[i] <-   x[i]*min(1,k/abs(x[i]))}
                       return(x_rob)}
    K2<-mean(huber(rnorm(10000,0,1))^2)
    
    
    gstable1=function(est){
                sigma.v<-est[1]
                sigma.e<-est[2]
                V<-matrix(0,n,n)
                V<-sigma.e*diag(1,n)+sigma.v*z%*%SpW%*%t(z)
                U<-diag(diag(V),n,n)
                svd.tmp.U=svd(U)
                s<-matrix(0,2,1)
                sqrt.U=t(svd.tmp.U$v%*%(t(svd.tmp.U$u)*sqrt(svd.tmp.U$d)))
                sqrt.U.inv<-solve(t(svd.tmp.U$v%*%(t(svd.tmp.U$u)*sqrt(svd.tmp.U$d))))
                Vi<-solve(V)
                y.star1=sqrt.U.inv%*%y
                x.star1=sqrt.U.inv%*%XS
                res.new<-y.star1-x.star1%*%betastim
                res.robust=my.psi(res.new,k)
                
                s[1,1]<-t(matrix(c(res.robust),n,1))%*%sqrt.U%*%Vi%*%z%*%t(z)%*%Vi%*%sqrt.U%*%matrix(c(res.robust),n,1)-K2*sum(diag(Vi%*%z%*%t(z)))
                s[2,1]<-t(matrix(c(res.robust),n,1))%*%sqrt.U%*%Vi%*%diag(1,n,n)%*%Vi%*%sqrt.U%*%matrix(c(res.robust),n,1)-K2*sum(diag(Vi%*%diag(1,n,n)))
                
                
                (s[1,1]^2)+(s[2,1]^2)
                
                }
                
    pstable=function(p_est){
                WW<-solve((diag(1,areanumber)-diag(p_est,areanumber)%*%W_matrix)%*%(diag(1,areanumber)-diag(p_est,areanumber)%*%t(W_matrix)))
                V<-matrix(0,n,n)
                V<-estsigma2e[iter+1]*diag(1,n)+estsigma2u[iter+1]*z%*%WW%*%t(z)
                U<-diag(diag(V),n,n)
                svd.tmp.U=svd(U)
                s<-matrix(0,2,1)
                sqrt.U=t(svd.tmp.U$v%*%(t(svd.tmp.U$u)*sqrt(svd.tmp.U$d)))
                sqrt.U.inv<-solve(t(svd.tmp.U$v%*%(t(svd.tmp.U$u)*sqrt(svd.tmp.U$d))))
                Vi<-solve(V)
                y.star1=sqrt.U.inv%*%y
                x.star1=sqrt.U.inv%*%XS
                res.new<-y.star1-x.star1%*%betastim
                res.robust=my.psi(res.new,k)
                diffV_p=estsigma2u[iter+1]*z%*%(-WW%*%(2*p_est*W_matrix%*%t(W_matrix)-2*W_matrix)%*%WW)%*%t(z)
                
                ss<--t(matrix(c(res.robust),n,1))%*%sqrt.U%*%(-Vi%*%diffV_p%*%Vi)%*%sqrt.U%*%matrix(c(res.robust),n,1)-K2*sum(diag(Vi%*%diffV_p))
                #ss<-t(matrix(c(res.robust),n,1))%*%sqrt.U%*%Vi%*%z%*%t(z)%*%Vi%*%sqrt.U%*%matrix(c(res.robust),n,1)-K2*sum(diag(Vi%*%z%*%t(z)))
                
                
                ((ss)^2)
                
                }
    
    #STEP 1    
          
    
          beta.q<-matrix(startbeta,pp,1)
          estsigma2u<-startsigma_v
          estsigma2e<-startsigma_e
          estp<-startwert_p
          
           
          diff.s<-1
          iter<-0
          while (abs(diff.s)>tol)
           {
                       iter<-iter+1
                       SpW<-solve((diag(1,areanumber)-diag(estp[iter],areanumber)%*%W_matrix)
                                  %*%(diag(1,areanumber)-diag(estp[iter],areanumber)%*%t(W_matrix)))
                       
                #STEP 1 Computation of beta
                       diff.b<-1
                       iter1<-0
                
                       tmp=estsigma2e[iter]*diag(1,n)+estsigma2u[iter]*z%*%SpW%*%t(z)
                       tmp.inv=solve(tmp)
                       tmp.U=diag(diag(tmp),n,n)
                       svd.tmp.U=svd(tmp.U)
                       sqrt.ui=solve(t(svd.tmp.U$v%*%(t(svd.tmp.U$u)*sqrt(svd.tmp.U$d))))
                       sqrt.u=t(svd.tmp.U$v%*%(t(svd.tmp.U$u)*sqrt(svd.tmp.U$d)))
                
                       while (abs(diff.b)>tol)
                           {iter1<-iter1+1
                           
                          
                           y.star=sqrt.ui%*%y
                           x.star=sqrt.ui%*%XS
                           
                           res1<-(y.star-x.star%*%beta.q)
                           w1<-diag(c(my.psi(res1,k)/res1),n,n)
                    
                           betastim=solve(t(XS)%*%tmp.inv%*%sqrt.u%*%w1%*%sqrt.ui%*%XS)%*%t(XS)%*%tmp.inv%*%sqrt.u%*%w1%*%sqrt.ui%*%y
                    
                           
                    
                           diff.b<-sum((betastim-beta.q)^2)
                    
                           beta.q<-betastim
                    
                           if (iter1>maxit)
                               {warning(paste("failed to converge in", maxit, "steps for beta"))
                               break}
                    
                            }
                      
                #STEP 2 Computation of sigma
                        
                     
                opt.stable=optim(c(estsigma2u[iter],estsigma2e[iter]),gstable1,method="Nelder-Mead")
                
                estsigma2u[iter+1]<-opt.stable$par[1]
                estsigma2e[iter+1]<-opt.stable$par[2]
                
                diff.s<-sum((estsigma2u[iter+1]-estsigma2u[iter])^2+(estsigma2e[iter+1]-estsigma2e[iter])^2)
                
                           
                #STEP 3 Computation of p
     
                #opt.p=optim(estp[iter],pstable,method="Nelder-Mead")
                #estp[iter+1]<-opt.p$par
  
                opt.p=optimize(pstable,lower = 0, upper = 1)
                #estp[iter+1]<-opt.p$objective
                estp[iter+1]<-opt.p$minimum
                      
                       
                #diff.s<-sum((estp[iter+1]-estp[iter])^2)
                
                       if(sum(estsigma2u[iter+1]+estsigma2e[iter+1])>var(y)){warning(paste("stopped after", maxit, "steps"))
                                                                             break}      
                       if (iter>maxit){warning(paste("failed to converge in", maxit, "steps for theta"))
                                       break}      
                       if (iter>maxit)
                       {warning(paste("failed to converge in", maxit, "steps for p"))
                        break}
                 
           }
           
    SpW_neu<-solve((diag(1,areanumber)-diag(estp[iter+1],areanumber)%*%W_matrix)%*%(diag(1,areanumber)-diag(estp[iter+1],areanumber)%*%t(W_matrix)))
    
    R.tmp=estsigma2e[iter+1]*diag(1,n)
    svd.R.tmp=svd(R.tmp)
    sqrt.R.tmp.inv=solve(t(svd.R.tmp$v%*%(t(svd.R.tmp$u)*sqrt(svd.R.tmp$d))))
    
    G.tmp=estsigma2u[iter+1]*SpW_neu
    svd.G.tmp=svd(G.tmp)
    sqrt.G.tmp.inv=solve(t(svd.G.tmp$v%*%(t(svd.G.tmp$u)*sqrt(svd.G.tmp$d))))
    
    tmp=estsigma2e[iter+1]*diag(1,n)+estsigma2u[iter+1]*z%*%SpW_neu%*%t(z)
    tmp.inv=solve(tmp)
    tmp.U=diag(diag(tmp),n,n)
    svd.tmp.U=svd(tmp.U)
    sqrt.ui=solve(t(svd.tmp.U$v%*%(t(svd.tmp.U$u)*sqrt(svd.tmp.U$d))))
    sqrt.u=t(svd.tmp.U$v%*%(t(svd.tmp.U$u)*sqrt(svd.tmp.U$d)))
    
    
    vv.tmp<-G.tmp%*%t(z)%*%solve(R.tmp+z%*%G.tmp%*%t(z))%*%as.vector(y-XS%*%beta.q)
    
    
    diff.u<-1
    iter2<-0
              while (abs(diff.u)>tol)
               {
               iter2<-iter2+1 
                v_robust=as.vector(vv.tmp)
                res1<-sqrt.R.tmp.inv%*%(y-XS%*%beta.q-z%*%v_robust)
                res2<-sqrt.G.tmp.inv%*%v_robust
                w2<-diag(c(my.psi(res1,k_v)/res1),n,n)
                w3<-diag(c(my.psi(res2,k_v)/res2),areanumber,areanumber)
                A=t(z)%*%(sqrt.R.tmp.inv)%*%w2%*%(sqrt.R.tmp.inv)%*%z+sqrt.G.tmp.inv%*%w3%*%sqrt.G.tmp.inv
                B<-t(z)%*%(sqrt.R.tmp.inv)%*%w2%*%(sqrt.R.tmp.inv)%*%(y-XS%*%beta.q)
                vv.tmp<-solve(A)%*%B
                diff.u<-sum(c((vv.tmp-v_robust)^2))
                if (iter2>maxit)
                    {warning(paste("failed to converge in", maxit, "steps"))
                    break}
                }
         
        rand.eff.robust=as.numeric(vv.tmp)
    
  #SREBLUP estimator: G is the variance-cov matrix for the random effects
  y.hat.sample<-NULL
  nii<-as.numeric(table(saind))
  if(length(x.out)==0){
    
    thetaSREBLUP<-(Xpop%*%beta.q) + rand.eff.robust
    
    
  }else{
    sample.sizer<-as.numeric(table(x.out$clusterid))
    Ni<-sample.sizer+nii
    Xmean_out<-cbind(1,tapply(x.out$x,x.out$clusterid,mean))
    y.sample<-tapply(y,saind,sum)
    thetaSREBLUP<-as.vector((Xpop%*%beta.q) + rand.eff.robust)
    thetaSREBLUP_out<-as.vector((Xmean_out%*%beta.q) + rand.eff.robust)
    y.hat.sample<-as.vector((1/Ni)*((y.sample)+sample.sizer*(thetaSREBLUP_out)))
  }
  
  #Bias-correction of the SREBLUP (locally, Schmid et al. (2013))
  random.effects<-rep(rand.eff.robust,times=c(nii))
  resid<-y-XS%*%beta.q-random.effects
  MAD<-as.numeric(tapply(X=resid,INDEX=saind,FUN=mad))
  resid_MAD<-resid/(rep(MAD,times=c(nii)))
  #bias_corr_local<-((1-ni/Ni)/ni)*tapply(X=rep(MAD,times=c(nii))*my.psi(resid_MAD,k_ccst),INDEX=saind,FUN=sum)
  bias_corr_local<-((1-ni/Ni)/ni)*tapply(X=rep(MAD,times=c(nii))*(resid_MAD)*psi.huber(resid_MAD,k_ccst),INDEX=saind,FUN=sum)   
  thetaSREBLUP_SCCST<-as.vector(thetaSREBLUP+bias_corr_local)
  y.hat.sample_SCCST<-as.vector(y.hat.sample+bias_corr_local) 
  
  
###################################################
  list(est_syn = thetaSREBLUP, est_sample = y.hat.sample, est_syn_SCCST = thetaSREBLUP_SCCST,
       est_sample_SCCST = y.hat.sample_SCCST, randeff = rand.eff.robust, beta = beta.q, 
       sigma_v = estsigma2u[iter+1], sigma_e = estsigma2e[iter+1], spatial_p = estp[iter+1],
       areasize = Ni, samplesize = nii, iterations = iter)
       
       
}
