#################################################################
###
### Example to calculate de SEBA estimator and the MSE 
### 
### Author:    Yolanda Marhuenda (y.marhuenda@umh.es)
### File name: Example_ModelA.R 
### Updated:   January 20th, 2011
###       
#################################################################


rm(list=ls(all.names=TRUE))                  # delete previous objects 

setwd("C:/PROYECTOS/SAMPLE/04_softwareSPFH") # set path where data sets 
                                             # and functions are

# Load library "MASS" and files where functions are
library("MASS")                                 
source("Henderson.R")
source("MessageErrorFitting.R")
source("FitSpatioTemporalFH.R")
source("BootMseSpatioTemporalFH.R") 


# Set the seed for mvrnorm function used in BootMSE.SpatioTemporalFH 
# function. This step is not necessary
set.seed(301)                                
                  
# Define data input parameters
model <- "A"                            # model: "A" or "B"
file_data <- "DataExample.txt"          # data input filename 
file_W    <- "WExample.txt"             # proximity matrix filename

                                        # Set the column numbers or labels
                                        # of variables included in 
                                        # file_data used in the model

c_dom       <- 1                        # area column number
c_period    <- 2                        # time instant column number
c_varaux    <- c("Intercept","X1","X2") # auxiliary variables labels
c_ydt       <- 6                        # direct estimator column number
c_sigma2dt  <- 7                        # direct estimator variance column 
                                        # number

nT          <- 3                        # total number of time instants
intconfidence <- 0.95                   # confidence interval to calculate
                                        # theta, beta and pvalue
MAXITER     <- 200                      # maximum number of iterations of 
                                        # Fisher-scoring algorithm
PRECISION   <- 0.0001                   # convergence tolerance limit for
                                        # the Fisher-scoring algorithm
nB          <- 100                      # number of bootstrap replications


# Read input data
datos    <- read.table(file=file_data, header=TRUE)
X        <- as.matrix(datos[,c_varaux]) # auxiliary variables
ydt      <- datos[,c_ydt]               # y
sigma2dt <- datos[,c_sigma2dt]          # variance
nD       <- nrow(datos)/nT              # number of areas
W        <- as.matrix(read.table(file=file_W,header=FALSE)) # proximity 
                                                            # matrix


#Obtain initial values of theta 
seedsigma <- Henderson(X,ydt,sigma2dt)
if (seedsigma<0)
{
   cat("Henderson sigma:", seedsigma,"<0 \n")
   seedsigma<-min(sigma2dt)
   cat("Warning: Henderson sigma is set to the minimum variance.",seedsigma)
   #quit()
}
sigma21 <- sigma22 <- 0.5*seedsigma 
rho1 <- 0.3   
rho2 <- 0.3                             # model B

if (model=="A") 
   theta0 <- rbind(sigma21, rho1, sigma22)
if (model=="B")   
   theta0 <- rbind(sigma21, rho1, sigma22, rho2)

# Fit the spatio-temporal Fay Herriot A or B
result <- FitSpatioTemporalFH(model,X,ydt,nD,nT,sigma2dt,theta0,W,MAXITER,
                                                   PRECISION,intconfidence)

# Validate output
if (result$convergence==FALSE || result$validtheta==FALSE)
{
   MessageErrorFitting(model,0,result$convergence,result$iterations,
                                            result$validtheta,result$theta)
   quit();
} 

# Print results
result$model         
result$convergence
result$iterations
result$validtheta
result$theta
result$beta
result$goodnessoffit

# Calculate mean squared error
mse <- BootMSE.SpatioTemporalFH(model,nB,X,nD,nT,sigma2dt,
                                result$beta[,"coef"],
                                result$theta[,"estimate"],rho1,rho2,
                                W,MAXITER,PRECISION,intconfidence)  

# Print direct estimates, variance and SEBA or SEBB estimates, mse and 
# residuals of the last period nT.
dom    <- datos[,c_dom]
period <- datos[,c_period]
output <- data.frame(Area=dom, Time=period, Direct=ydt, 
                     Model=result$estimates, VarDirect=sigma2dt, 
                     MSEModel=mse, residuals=ydt-result$estimates) 
print (output[output[,"Time"]==nT,], row.names=FALSE)



