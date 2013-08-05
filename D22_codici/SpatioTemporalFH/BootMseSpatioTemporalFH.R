#####################################################################
###           Spatio-temporal Fay Herriot Models
###                   SAMPLE Project
###
### Author:    Yolanda Marhuenda (y.marhuenda@umh.es)
### File name: BootMseSpatioTemporalFH.R 
### Updated:   January 20th, 2011
###       
#####################################################################

BootMSE.SpatioTemporalFH <- function(model,nB,Xdt,nD,nT,sigma2dt,
                                     beta,theta,rho1_0_b,rho2_0_b,
                                     W,MAXITER,PRECISION,confidence)
{
   if (model!="A" && model!="B")
   {
      print("Error: Model must be A or B.")
      return;    
   }
   msedt <- 0
   Id    <- diag(1,nrow=nD, ncol=nD)

   sigma21 <- theta[1]
   rho1    <- theta[2]
   sigma22 <- theta[3]

   Omega1rho1 <- solve(crossprod(Id-rho1*W))  
   sigma21Omega1rho1 <- sigma21*Omega1rho1
   M <- nD*nT
 
   if (model=="B")
   {
      rho2   <- theta[4]
      Unomenrho22_05 <- (1-rho2^2)^(-0.5)
      u2dt_b <- matrix(0, nrow=M, ncol=1)
   }

   b<-1
   while (b<=nB)
   {
      #cat("B:", b ,"Start:", date(), "\n")
      u1_b   <- matrix(data=mvrnorm(n=1, mu=rep(0,nD),
                       Sigma=sigma21Omega1rho1), nrow=nD, ncol=1)
      u1dt_b <- matrix(data=rep(u1_b, each=nT),nrow=M,ncol=1)   

      if (model=="A")       
         u2dt_b <- matrix(data=rnorm(M, mean=0, sd=sqrt(sigma22)), 
                          nrow=M, ncol=1) 
      else
      {
         epsilondt_b <- matrix(data=rnorm(M, mean=0, sd=sqrt(sigma22)),
                               nrow=M, ncol=1)
         i <- 1
         for (d in 1:nD)
         {
            u2dt_b[i] <- Unomenrho22_05*epsilondt_b[i]
            for (t in 2:nT)
            {
               i <- i+1
               u2dt_b[i] <- rho2*u2dt_b[i-1]+epsilondt_b[i]
            }
            i<-i+1
         }
      }

      edt_b  <- matrix(data=rnorm(M, mean=0, sd=sqrt(sigma2dt)), 
                       nrow=M,ncol=1) 
      ydt_b  <- Xdt%*%beta + u1dt_b + u2dt_b + edt_b
      mudt_b <- ydt_b - edt_b			

                                               ### fitting the model  
      seedsigma_b <- Henderson(Xdt,ydt_b,sigma2dt)
      if (seedsigma_b<0)
      {
         cat("Warning: Assigning Henderson seed for sigma in 
              bootstrap sample b=", b,".\n")
         seedsigma_b<-min(sigma2dt)
         cat("sigma is established to the minimum value of sigma2dt:",
              seedsigma_b,".\n")
      }
      sigma21_0_b <- sigma22_0_b <- 0.5*seedsigma_b 

      if (model=="A")  
         theta0_b <- rbind(sigma21_0_b, rho1_0_b, sigma22_0_b)
      else
         theta0_b <- rbind(sigma21_0_b, rho1_0_b, sigma22_0_b, rho2_0_b)

      result <- FitSpatioTemporalFH(model,X,ydt_b,nD,nT,sigma2dt,theta0_b,
                                    W,MAXITER,PRECISION,confidence) 

      if (result$convergence==FALSE || result$validtheta==FALSE)
      {
         ##cat("Bootstrap sample:(",b,")\n")
         ##MessageErrorFitting(model,b,result$convergence,result$iterations,
         ##                            result$validtheta,result$theta)
         ##cat("Generate other sample\n")
         next
      } 
                                                ### calculate msedt
      difference <- result$estimates - mudt_b
      msedt      <- msedt + difference^2
      b <- b+1     
   } #while (b<=nB)
      
   msedt <- msedt/nB
   return (msedt)
}












