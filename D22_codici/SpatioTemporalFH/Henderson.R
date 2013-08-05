#####################################################################
###           Spatio-temporal Fay Herriot Models
###                   SAMPLE Project
###
### Author:    Yolanda Marhuenda (y.marhuenda@umh.es)
### File name: Henderson.R 
### Updated:   January 20th, 2011
###       
#####################################################################

Henderson <- function(X,y,sigma2dt)
{
   invVe    <- 1/sigma2dt       
   invVeX   <- invVe*X
   tXinvVeX <- t(X)%*%invVeX
   P        <- diag(invVe)-invVeX%*%solve(tXinvVeX)%*%t(invVeX)
   trtZPZ   <- sum(diag(P)) 
   tyPy     <- t(y)%*%P%*%y    # for vectors t(y)%*%P%*%y=sum(y*P%*%y)
   sigma2   <- as.numeric((tyPy-(nrow(X)-ncol(X)))/trtZPZ)
  
   return (sigma2)
}



