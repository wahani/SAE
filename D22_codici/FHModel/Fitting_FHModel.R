###############################################################################
###
### This function fits a Fay-Herriot model 
### Fitting method can be chosen between REML and FH methods
###
### Work for European project SAMPLE
###
### Author: Isabel Molina Peralta
### File name: Fitting_FHModel.R
### Updated: Example
###
###############################################################################

fitFH<-function(X,y,Dvec,method="REML",MAXITER=500) {

  m<-length(y) # Sample size or number of areas
  p<-dim(X)[2] # Num. of X columns of num. of auxiliary variables (including intercept)
  Xt<-t(X)
  
  # Fisher-scoring algorithm for REML estimator of variance A starts
  
  if (method=="REML") {    

    # Initial value of variance A is fixed to the median of sampling variances Dvec
    Aest.REML<-0
    Aest.REML[1]<-median(Dvec)
    
    k<-0
    diff<-1
    while ((diff>0.0001)&(k<MAXITER))
    {
      k<-k+1
      Vi<-1/(Aest.REML[k]+Dvec)
      XtVi<-t(Vi*X)
      Q<-solve(XtVi%*%X)
      P<-diag(Vi)-t(XtVi)%*%Q%*%XtVi
      Py<-P%*%y
      # Score function obtained from restricted log-likelihood
      s<-(-0.5)*sum(diag(P))+0.5*(t(Py)%*%Py) 
      # Fisher information obtained from restricted log-likelihood 
      F<-0.5*sum(diag(P%*%P))                  
      # Updating equation
      Aest.REML[k+1]<-Aest.REML[k]+s/F
      # Relative difference of estimators in 2 iterations for stopping condition
      diff<-abs((Aest.REML[k+1]-Aest.REML[k])/Aest.REML[k])
    } # End of while

    # Final estimator of variance A
    A.REML<-max(Aest.REML[k+1],0)
    #print(Aest.REML)
    
    # Indicator of convergence
  
    if(k<MAXITER) {conv<-TRUE} else {conv<-FALSE}
  
    # Computation of the coefficients'estimator beta
  
    Vi<-1/(A.REML+Dvec)
    XtVi<-t(Vi*X)
    Q<-solve(XtVi%*%X)
    beta.REML<-Q%*%XtVi%*%y
  
    # Significance of the regression coefficients
  
    varA<-1/F
  
    std.errorbeta<-sqrt(diag(Q))
    tvalue<-beta.REML/std.errorbeta
    pvalue<-2*pnorm(abs(tvalue),lower.tail=FALSE)
  
    # Goodness of fit measures: loglikelihood, AIC, BIC
    
    Xbeta.REML<-X%*%beta.REML
    resid<-y-Xbeta.REML
    
    loglike<-(-0.5)*(sum(log(2*pi*(A.REML+Dvec))+(resid^2)/(A.REML+Dvec)))
    AIC<-(-2)*loglike+2*(p+1)
    BIC<-(-2)*loglike+(p+1)*log(m)
  
    goodness<-c(loglike=loglike,AIC=AIC,BIC=BIC)
    
    # Computation of the empirical best (EB) predictor
  
    thetaEB.REML<-Xbeta.REML+A.REML*Vi*resid
    coef<-data.frame(beta.REML,std.errorbeta,tvalue,pvalue)
  
    return(list(convergence=conv,modelcoefficients=coef,variance=A.REML,goodnessoffit=goodness,EBpredictor=thetaEB.REML))
 
    # Fisher-scoring algorithm for REML estimator of variance A starts
 
  } else if (method=="FH") {
  
    # Initial value of variance A is fixed to the median of sampling variances Dvec
    Aest.FH<-NULL
    Aest.FH[1]<-median(Dvec)
  
    k<-0
    diff<-1
    while ((diff>0.0001)&(k<MAXITER)){
      k<-k+1
      Vi<-1/(Aest.FH[k]+Dvec)
      XtVi<-t(Vi*X)
      Q<-solve(XtVi%*%X)
      betaaux<-Q%*%XtVi%*%y
      resaux<-y-X%*%betaaux   
      # Left-hand side of equation for FH estimator
      s<-sum((resaux^2)*Vi)-(m-p)
      # Expectation of negative derivative of s 
      F<-sum(Vi)
      # Updating equation
      Aest.FH[k+1]<-Aest.FH[k]+s/F
      # Relative difference of estimators in 2 iterations for stopping condition
      diff<-abs((Aest.FH[k+1]-Aest.FH[k])/Aest.FH[k])
  
    } # End of while
  
    A.FH<-max(Aest.FH[k+1],0)
    print(Aest.FH)
  
    # Indicator of convergence
  
    if(k<MAXITER) {conv<-TRUE} else {conv<-FALSE}
    
    # Computation of the coefficients'estimator beta
  
    Vi<-1/(A.FH+Dvec)
    XtVi<-t(Vi*X)
    Q<-solve(XtVi%*%X)
    beta.FH<-Q%*%XtVi%*%y
  
    # Significance of the regression coefficients
  
    varA<-1/F
    varbeta<-diag(Q)
    std.errorbeta<-sqrt(varbeta)
    zvalue<-beta.FH/std.errorbeta
    pvalue<-2*pnorm(abs(zvalue),lower.tail=FALSE)
  
    # Goodness of fit measures: loglikelihood, AIC, BIC
    
    Xbeta.FH<-X%*%beta.FH
    resid<-y-Xbeta.FH
  
    loglike<-(-0.5)*(sum(log(2*pi*(A.FH+Dvec))+(resid^2)/(A.FH+Dvec)))
    AIC<-(-2)*loglike+2*(p+1)
    BIC<-(-2)*loglike+(p+1)*log(m)
  
    goodness<-c(loglike=loglike,AIC=AIC,BIC=BIC)
  
    # Computation of the empirical best (EB) predictor
  
    thetaEB.FH<-Xbeta.FH+A.FH*Vi*resid
    coef<-data.frame(beta.FH,std.errorbeta,zvalue,pvalue)
  
    return(list(convergence=conv,modelcoefficients=coef,variance=A.FH,goodnessoffit=goodness,EBpredictor=thetaEB.FH))
    
    # Error printing when method is different from REML of FH.
  }  else { print("Error: Unknown method") }
  
}
