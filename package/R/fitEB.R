#' fitSTEBLUP
#' 
#' @description Wrapper function for the estimation of a FH model. Uses the
#' function \code{fitFH} provided by Sample Deliverable 22: Software on Small Area Estimation
#' 
#' @param formula formula for the fixed effects part
#' @param dat data with dependent and independent variable as they are used in formula. The data should 
#' contain a variable \code{Domain} and \code{Time} indicating the domain and time period for each observation respectively - only
#' the last time period will be used for estimation and prediction.
#' @param beta, sigma, rho start values for the parameters - not needed in the call to \code{fitFH}
#' @param sigmaSamplingError True variances for the sampling errors
#' @param nDomains the number of domains. Will be determined by default, see \code{dat} - not needed in the call to \code{fitFH}
#' @param nTime the number of time periods. Will be determined by default, see \code{dat} - not needed in the call to \code{fitFH}
#' @param w0 proximity matrix - not needed in the call to \code{fitFH}
#' @param w at the moment row-standardized proximity matrix - not needed in the call to \code{fitFH}
#' @param tol numerical tolerance for algorithm - not needed in the call to \code{fitFH}
#' @param method method used by the function \code{\link{optim}} - not needed in the call to \code{fitFH}
#' @param maxIter maximal number of iterations
#'
#'@examples require(spatioTemporalData); require(SAE)
#' set.seed(3)
#' setup <- simSetupMarhuenda(nDomains=30, nTime=5, sarCorr=0.5, arCorr=0.5, n = 200)
#' output <- simRunMarhuenda(setup) 
#' dat <- slot(output[[1]], "data")[[1]]
#' result <- fitEB(y~x, dat, c(0,1), c(1,1), c(0.5,0.5))
#' #summary(result)
fitEB <- function(formula, dat, beta, sigma, rho, 
                  sigmaSamplingError = seSigmaClosure(nDomains, nTime)(),
                  nDomains = getNDomains(dat),
                  nTime = getNTime(dat),
                  w0 = w0Matrix(nDomains), 
                  w = wMatrix(nDomains),
                  tol = 1e-3, method = "Nelder-Mead", maxIter = 500) {
  
  dat <- subset(dat, Time == max(Time))
  sigmaSamplingError <- sigmaSamplingError[dat$Time == max(dat$Time)]
  
  modelSpecs <- prepareData(formula, dat, nDomains, nTime, beta, sigma, rho, sigmaSamplingError, w0, w, tol, method, maxIter)
  
  modelFit <- fitFH(modelSpecs$x, modelSpecs$y, modelSpecs$sigmaSamplingError, method = "REML", MAXITER = modelSpecs$maxIter)
  
  output <- list(estimates = data.frame(yHat = modelFit$EBpredictor),
                 beta = modelFit$modelcoefficients,
                 sigma = modelSpecs$sigma,
                 rho = modelSpecs$rho)
  class(output) <- c("fitEB", "fitSTREBLUP", "list")
  output
}


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
