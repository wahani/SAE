###############################################################################
###
### This function gives the MSE estimator of the EB estimator under a FH model
### The EB estimator is obtained either by REML or by FH fitting methods
###
### Work for European project SAMPLE
###
### Author: Isabel Molina Peralta
### File name: MSE_FHModel.R
### Updated: March 15th, 2010
###
###############################################################################

MSE.FHmodel<-function(X,Dvec,A,method="REML"){

  m<-dim(X)[1] # Sample size or number of areas
  p<-dim(X)[2] # Num. of X columns of num. of auxiliary variables (including intercept)

  # Initialize vectors containing the values of g1-g3 for each area along with the final mse 
  g1d<-rep(0,m)
  g2d<-rep(0,m)
  g3d<-rep(0,m)
  mse2d<-rep(0,m)
  
  # Elements of the inverse covariance matrix in a vector
  Vi<-1/(A+Dvec)
  # Auxiliary calculations
  Bd<-Dvec/(A+Dvec)
  SumAD2<-sum(Vi^2)
  XtVi<-t(Vi*X)
  Q<-solve(XtVi%*%X)

  # Calculation of g1-g3 and final MSE when fitting method is REML

  if (method=="REML"){
  
    # Asymptotic variance of REML estimator of variance A
    VarA<-2/SumAD2
  
    for (d in 1:m){
      g1d[d]<-Dvec[d]*(1-Bd[d])
      xd<-matrix(X[d,],nr=1,nc=p)
      g2d[d]<-(Bd[d]^2)*xd%*%Q%*%t(xd)
      g3d[d]<-(Bd[d]^2)*VarA/(A+Dvec[d])
      mse2d[d]<-g1d[d]+g2d[d]+2*g3d[d]
    }

  return(mse=mse2d)  

  # Calculation of g1-g3 and final MSE when fitting method is FH
    
  } else if (method=="FH") {
  
    SumAD<-sum(Vi)
    # Asymptotic variance of FH estimator of variance A
    VarA<-2*m/(SumAD^2)

    # Asymptotic bias of FH estimator of A
    b<-2*(m*SumAD2-SumAD^2)/(SumAD^3)
    
    for (d in 1:m){  
      g1d[d]<-Dvec[d]*(1-Bd[d])
      xd<-matrix(X[d,],nr=1,nc=p)                                              
      g2d[d]<-(Bd[d]^2)*xd%*%Q%*%t(xd)
      g3d[d]<-(Bd[d]^2)*VarA/(A+Dvec[d])
      mse2d[d]<-g1d[d]+g2d[d]+2*g3d[d]-b*(Bd[d]^2)
    }
  
  return(mse=mse2d)  
  
  # Error printing when fitting method is different from REML of FH.
  }  else { print("Error: Unknown method") }

}
