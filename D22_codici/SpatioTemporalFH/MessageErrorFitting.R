#####################################################################
###           Spatio-temporal Fay Herriot Models
###                   SAMPLE Project
###
### Author:    Yolanda Marhuenda (y.marhuenda@umh.es)
### File name: MessageErrorFitting.R 
### Updated:   January 20th, 2011
###       
#####################################################################

MessageErrorFitting <- function(model,nsample,convergence,niter,
                                validtheta,theta)
{
   linea <- paste("Warning: Sample=",nsample,"\n")
   cat(linea)

   if (convergence==FALSE)
   {
      linea <- paste("Maximum number of iterations is reached.\n")
      cat(linea)
   }
   else 
   {
      sigma21_REML <- theta[1]
      rho1_REML    <- theta[2]
      sigma22_REML <- theta[3]

      if (sigma21_REML<0)
      {
         linea <- paste(" sigma21<0 (",sigma21_REML,")" )
         cat(linea)
      }
      if (rho1_REML<(-1) || rho1_REML>1)
      {
         linea <- paste(" rho1 must be in [-1,1] (",rho1_REML,")" )
         cat(linea)
      }
      if (sigma22_REML<0)
      {
         linea <- paste(" sigma22<0 (",sigma22_REML,")" )
         cat(linea)
      }

      if (model=="B")
      {
         rho2_REML <- theta[4]
         if (rho2_REML<(-1) || rho2_REML>1)
         {
            linea <- paste(" rho2 must be in [-1,1] (",rho2_REML,")" )
            cat(linea)
         }
      }   
      linea <- paste("\nTotal number of iterations:",niter,".\n")
      cat (linea)
   }
} 

